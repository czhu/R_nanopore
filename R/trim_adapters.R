get_adapter_hit = function(x, minMatchLength=6,minPid=60){
    ## we favor adapter that's close to the begining of the read -> we penalise on hits that's
    ## far from the start to resolve borderline case
    x %>% filter(match_length>minMatchLength, pid>minPid) %>% mutate(pos_weighted_score=score-start%/%100) %>%
    group_by(group) %>% summarise(
            hit=which.max(pos_weighted_score),
            start=start[hit],
            end=end[hit],
            score=score[hit],
            match_length=match_length[hit],
            pid=pid[hit],
            match_seq_read=match_seq_read[hit],
            adapter_id=adapter_id[hit]
        )
}

reduce_column = function(x,columanToRemove=c("group","hit")){
    ## remove extra column from the table
    x[,setdiff(colnames(x),columanToRemove)]
}

## current pipeline to trim reads
trim_adapters = function(infile, adapters, outfolder="workspace",ncpu=20,readsChunkSize = 1e5,
    minMedianQualityScore = 1, minReadLength = 100, middleAdaptersTolerance=0.15, windowSize = 300) {
    ## infile -> path to fastq file
    ## adapter -> named DNAStringSet
    polyADistanceToAdapter = 50

    ## output files
    outFileFilteredReads = file.path(outfolder,"reads_qc_failed.fastq.gz")
    outFileProcessedReads = file.path(outfolder,"reads_processed.gz")
    processlogFile = file.path(outfolder,"pipeline.log")
    summaryStatFile = file.path(outfolder,"summary_stat.txt")

    ## initialise
    nRead = 0
    nFiltered = 0
    nReadWithAdapters = 0

    ### adapter sequences
    #adapters = readDNAStringSet("/g/steinmetz/czhu/pore/dev_adapter_trimming/adapters/adapter_Set2.fasta")
    #revAdapters = reverseComplement(adapters)
    #names(revAdapters) = paste0(names(revAdapters),"-rev")

    ####################
    ## read in chunks with FastqStream
    fs= FastqStreamer(infile,n=readsChunkSize)
    sink(file=file(processlogFile, open = "wt"),type="message")
    while (length(fq <- yield(fs))) {
        message("Processing ", length(fq)," Reads starting at read ", nRead+1)
        ##################
        ### median quality and short reads
        ##################
        seqStats = tibble(
            seqnames = sapply(strsplit(as.character(ShortRead::id(fq))," "),"[",1),
            median_qscore = median(get_quality(x)),
            read_length = width(fq),
        ) %>% mutate(qc_pass=median_qscore >= minMedianQualityScore & read_length >= minReadLength)

        ## define trimming window
        trimLeftWindow = cbind(1,ifelse(windowSize>=width(fq), width(fq),windowSize))
        trimRightWindow = cbind(ifelse(width(fq)-windowSize < 1, 1,width(fq)-windowSize), width(fq))

        ## split data
        readIndex = 1:length(fq)
        readIndexChuncks = split(readIndex,cut(readIndex, ncpu))

        ####### parallel here
        res = mclapply(1:ncpu, function(i){
            thisIndexChunck = readIndexChuncks[[i]]

            ## find polyA/T
            polyAT = find_polyA(fq[thisIndexChunck])
            polyAT$seqnames = seqStats$seqnames[thisIndexChunck][polyAT$group]

            ## find adapter on read start
            leftAdapters = match_adapter(
                narrow(fq[thisIndexChunck],trimLeftWindow[thisIndexChunck,1],trimLeftWindow[thisIndexChunck,2]),
                adapters, use.names=TRUE,avoidPartialMappingOnEnd=TRUE
            )
            leftAdapters = get_adapter_hit(leftAdapters)
            leftAdapters$seqnames = seqStats$seqnames[thisIndexChunck][leftAdapters$group]

            ## find adapter on read ends
            rightAdapters = match_adapter(
                ## we reverseComplement the read not the adapter, becase we penalise on the adapters being
                ## far from the end
                reverseComplement(narrow(fq[thisIndexChunck],trimRightWindow[thisIndexChunck,1],trimRightWindow[thisIndexChunck,2])),
                adapters, use.names=TRUE,avoidPartialMappingOnEnd=TRUE
            )
            rightAdapters = get_adapter_hit(rightAdapters)
            rightAdapters$seqnames = seqStats$seqnames[thisIndexChunck][rightAdapters$group]
            ## convert back to read coordinates
            readRightAdapterStart = rightAdapters$start
            readRightAdapterEnd = rightAdapters$end
            rightAdapters$start =  seqStats$read_length[thisIndexChunck][rightAdapters$group] - readRightAdapterEnd + 1
            rightAdapters$end =  seqStats$read_length[thisIndexChunck][rightAdapters$group] - readRightAdapterStart + 1

            ## find all middle adapters
            ## only for those that we have enough space between the two ends + 50 bp
            # hasMiddleWindow = trimRightWindow[thisIndexChunck,1] > (trimLeftWindow[thisIndexChunck,2] + 50)
            # middleAdapters = find_middle_adapter(
            #     narrow( fq[thisIndexChunck][hasMiddleWindow],
            #             trimLeftWindow[thisIndexChunck,2][hasMiddleWindow],
            #             trimRightWindow[thisIndexChunck,1][hasMiddleWindow]),
            #     c(adapters,revAdapters),
            #     tol=middleAdaptersTolerance,posShift=windowSize-1
            #     )
            # if(!is.null(middleAdapters)) {
            #     middleAdapters$seqnames = seqStats$seqnames[thisIndexChunck][hasMiddleWindow][middleAdapters$group]
            #     middleAdapters = middleAdapters[,setdiff(colnames(middleAdapters),c("group","hit","start",
            #     "end","adapter_id"))]
            #     colnames(middleAdapters) = paste0(colnames(middleAdapters),".middle")
            # }

            ### output the table -> make sure always the same column

            colnames(polyAT) = paste0(colnames(polyAT),".polyA")
            leftAdapters = reduce_column(leftAdapters)
            rightAdapters = reduce_column(rightAdapters)
            ans = full_join(leftAdapters,rightAdapters,by="seqnames",suffix=c(".left",".right")) %>%
                full_join(polyAT,by=c("seqnames"="seqnames.polyA"))
            # if(!is.null(middleAdapters)) {
            #     ans = ans %>% full_join(middleAdapters,by=c("seqnames"="seqnames.middle"))
            # } else{
            #     ans[,c("allStart.middle", "allEnd.middle", "allWidth.middle","allAdapter_id.middle")] = as.character(NA)
            # }
            ans
            }, mc.cores=ncpu)
        res= do.call(rbind,res)

        finalRes = left_join(seqStats,res,by="seqnames") %>% mutate(
            read_config_pass =  start.right > end.left & (
                ## theoretically this should be big than zero ...
                (start.polyA - end.left)>polyADistanceToAdapter |
                (start.right - end.polyA)>polyADistanceToAdapter),
            read_config=paste(
                ifelse(is.na(start.left),"N","P"),
                ifelse(is.na(adapter_id.polyA),"N", ifelse(adapter_id.polyA=="polyA","A","T")),
                ifelse(is.na(start.right),"N","P"),
                sep="-"
                ),
            trimStart = ifelse(is.na(end.left),1,end.left),
            trimEnd = ifelse(is.na(start.right),read_length,start.right)
                )
        ## if the right adapter is after the left adapter the matching has wrong -> no trimming
        isWrongConfig = which(finalRes$end.left > finalRes$start.right)
        finalRes$read_config[isWrongConfig] = "N-N-N"
        finalRes$trimStart[isWrongConfig] = 1
        finalRes$trimEnd[isWrongConfig] = finalRes$read_length[isWrongConfig]
        #all(finalRes$start.right> finalRes$end.left,na.rm=TRUE)

        #####################output
        ### write out the stats
        write_tsv(finalRes,summaryStatFile, append=file.exists(summaryStatFile))

        ### write out filtered reads
        if(sum(!finalRes$qc_pass)>0){
            writeFastq(fq[!finalRes$qc_pass],file=outFileFilteredReads,
                mode=ifelse(file.exists(outFileFilteredReads),"a","w"), compress=TRUE)
        }

        ###################### write out processed reads
        ## trimming
        trimmedFq = narrow(fq,finalRes$trimStart, finalRes$trimEnd)
        outSeq = ShortReadQ(
            sread=sread(trimmedFq), quality = quality(trimmedFq),
            id = BStringSet(paste(finalRes$seqnames,finalRes$read_config,sep=":")))

        writeFastq(outSeq[finalRes$qc_pass],
            file=outFileProcessedReads, mode=ifelse(file.exists(outFileProcessedReads),"a","w"), compress=TRUE)

        ################
        nRead = nRead + length(fq)
        ## make sure all left is left than right adapter
        message("Filtered ",sum(!finalRes$qc_pass)," due to QC")
        nFiltered = nFiltered + sum(!finalRes$qc_pass)
        message(sum(finalRes$qc_pass & !finalRes$read_config_pass,na.rm=TRUE),
            " passed QC wihtout the correct adapter configuration")

        thisNReadWithAdatpers = sum(grepl("P-.-P",finalRes$read_config) & finalRes$qc_pass)
        message(thisNReadWithAdatpers," having adapters on both ends")
        nReadWithAdapters = nReadWithAdapters + thisNReadWithAdatpers
    } ## loop end here
    ####
    message("Processed ", nRead, " in total\n",
        nFiltered," filtered due to QC\n",
        nRead - nFiltered," retained after QC\n",
        nReadWithAdapters," having adapters on both ends\n",
        nReadWithAdapters/nRead," raw yield\n",
        nReadWithAdapters/(nRead - nFiltered)," net yield\n"
    )
    close(fs)
    message(done)
}
