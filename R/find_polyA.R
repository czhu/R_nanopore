find_polyA = function(myseqs, n=12, max.mismatch=1L, use.names=TRUE) {
    ## myseqs shortread class
    ## n number of polyA/T
    mypaste = function(y) paste(y,collapse=",")

    adapters = DNAStringSet(c(
        polyA=paste(rep("A",12),collapse=""),
        polyT=paste(rep("T",12),collapse=""))
        )

    res = lapply(1:length(adapters), function(i) {
        ans = as_tibble(find_polyN(sread(myseqs), adapters[[i]]),max.mismatch=max.mismatch)
        if(nrow(ans)>0) ans$adapter_id = i
        ans
    })
    res = do.call(rbind,res)
    if(nrow(res)==0) {return(NULL)}
    res = res %>% arrange(group,start)
    if(use.names) res$adapter_id = names(adapters)[res$adapter_id]

    res %>% group_by(group) %>% summarise(
        hit = which.max(width),
        allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width),
        allAdapter_id = mypaste(adapter_id),
        start = start[hit],
        end = end[hit],
        adapter_id = adapter_id[hit]
        )
}
