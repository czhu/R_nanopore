find_middle_adapter = function(myseqs, adapters, tol=0.2, use.names=TRUE,posShift=0) {
    ## myseqs ShortRead class
    ## adapters DNAStringSet
    mypaste = function(y) paste(y,collapse=",")

    mymatch_func = function(subject,pattern, max.mismatch) {
        vmatchPattern_direct(pattern,subject, max.mismatch=max.mismatch,
            min.mismatch=0L, with.indels=TRUE, fixed=c(TRUE,TRUE),algo="auto")
    }

    res = lapply(1:length(adapters), function(i) {
        max.mismatch = as.integer(floor(nchar(adapters[[i]])*tol))
        ans = as_tibble(mymatch_func(sread(myseqs), adapters[[i]],max.mismatch))
        if(nrow(ans)>0) {ans$adapter_id = i}
        ans
    })
    res = do.call(rbind,res)
    if(nrow(res)==0) {return(NULL)}
    res = res %>% arrange(group,start)
    if(use.names) res$adapter_id = names(adapters)[res$adapter_id]
    ## shift position in given bp, to convert from coordinates in truncated reads back to
    ## full reads
    res$start = res$start + posShift
    res$end = res$end + posShift

    res %>% group_by(group) %>% summarise(
        hit = which.max(width),
        allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width),
        allAdapter_id = mypaste(adapter_id),
        start = start[hit],
        end = end[hit],
        adapter_id = adapter_id[hit]
        )
}
