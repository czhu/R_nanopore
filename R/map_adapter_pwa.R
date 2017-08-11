map_adapter_pwa = function(fq,adapters, windowSize = 200,mapEnd = FALSE, minLength = 10, ncpu=20,
    outfile) {
    ### caveat only find the first hit!!!
    require(ShortRead)
    require(Biostrings)

    readNames = sapply(strsplit(as.character(ShortRead::id(fq))," "),function(x)x[1])

    ## get the ends

    if(mapEnd){
        ## mapping the three end of the read
        mypos = ifelse(width(fq) < windowSize+1, 1, width(fq) - windowSize)
        mysubseq = narrow(fq, mypos ,width(fq))
    } else {
        ## five end
        mysubseq = narrow(fq,1,ifelse(width(fq) < windowSize + 1, width(fq), windowSize))
    }


    hits = mclapply(1:length(mysubseq), function(i){
        if(!(i %% 1000)) message("Process ",i," samples\n")
        #cat(i,"\n")
        pwa = pairwiseAlignment(adapters,
            sread(mysubseq)[i],type="overlap",
            gapOpening=4,gapExtension=2)

        if(any( Biostrings::nchar(pwa) > minLength)){
            mys = score(pwa)
            ## chooe the one with highest score if that's also the one with the higest score
            ## otherwise the longest alignment
            ## problem, short alignment from end better than full alignment
            wh = ifelse(which.max(nchar(pwa)) == which.max(mys),which.max(mys),which.max(nchar(pwa)))
            rv = list(
                readName = readNames[i],
                adapter = names(adapters)[wh],
                start = start(subject(pwa[wh])),
                end = end(subject(pwa[wh])),
                score = mys[wh],
                matchLength = nchar(pwa)[wh],
                pid = pid(pwa)[wh]
                )
            return(rv)
        } else {
            return(list(
                readName = readNames[i],
                adapter = as.character(NA),
                start = as.numeric(NA),
                end = as.numeric(NA),
                score = as.numeric(NA),
                matchLength = as.numeric(NA),
                pid = as.numeric(NA)
                ))
        }
    },mc.cores=ncpu)
    rv = do.call(rbind.data.frame,hits)
    rv$readName = as.character(rv$readName)
    rv$adapter = as.character(rv$adapter)

    if(mapEnd){
        rv$start = rv$start + mypos
        rv$end = rv$end + mypos
    }
    if(!missing(outfile)){
        toAppend = ifelse(file.exists(outfile),TRUE,FALSE)
        write.table(rv,outfile,append=toAppend,quote=FALSE,row.names=FALSE,
            col.names=!toAppend, sep="\t")
    } else {
        return(rv)
    }
}

summarise_map_adapter_pwa = function(fwdFile,revFile){
    ## summarise the results from identyfing adapters using pwa
    require(tidyverse)

    myfwd = read_tsv(fwdFile)
    myrev = read_tsv(revFile)

    message(mean(is.na(myfwd$adapter))," fwd adapter")
    message(mean(is.na(myrev$adapter))," rev adapter")

    inner_join(myfwd,myrev,by=c("readName"="readName"))
}
