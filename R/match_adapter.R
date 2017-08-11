match_adapter = function(myseqs, adapters,use.names=TRUE,avoidPartialMappingOnEnd=FALSE){
    ## use.names -> use name of adapter instead of index
    ## loop over adapters
    pwa_2_tibble = function(x) {
        tibble(
            start = start(pattern(x)),
            end = end(pattern(x)),
            group = 1:length(x),
            score = score(x),
            match_length = nchar(x),
            pid = pid(x),
            match_seq_read = as.character(pattern(x)),
            match_seq_adapter = as.character(subject(x)))
    }

    ### avoid partial at the end by appending long enough N
    ## -> only allow partial match at the begining
    ## it could happen that one or two N is considered as match, small artefact, but fine for us
    if(avoidPartialMappingOnEnd) {
        myseqs = ShortReadQ(
                sread=DNAStringSet(paste0(as.character(sread(myseqs)),paste(rep("N",100),collapse=""))),
                quality = BStringSet(paste0(as.character(PhredQuality(quality(myseqs))),paste(rep("I",100),collapse=""))),
                id = ShortRead::id(myseqs))
    }

    res = lapply(1:length(adapters), function(i) {
        pwa = pairwiseAlignment(sread(myseqs), adapters[[i]],patternQuality=PhredQuality(quality(myseqs)),
              type="overlap", gapOpening=4,gapExtension=2)
        ans = pwa_2_tibble(pwa)
        if(nrow(ans)>0) {ans$adapter_id = i}
        ans
    })
    res = do.call(rbind,res)
    res = res %>% arrange(group,adapter_id)
    if(use.names) res$adapter_id = names(adapters)[res$adapter_id]

    res
}
