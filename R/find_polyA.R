find_polyA = function(myseqs, n=12, max.mismatch=1L,hitBy="length",use.names=TRUE) {
    ## hitBy either by length (longest polyA/T is the hit)
    ## or distance to the end (shortest to the end is the hit)

    ## myseqs shortread class
    ## n number of polyA/T
    mypaste = function(y) paste(y,collapse=",")

    adapters = DNAStringSet(c(
        polyA=paste(rep("A",n),collapse=""),
        polyT=paste(rep("T",n),collapse=""))
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

    if(hitBy=="length"){
        rv = res %>% group_by(group) %>% summarise(
            hit = which.max(width),
            start = start[hit],
            end = end[hit],
            allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width),
            allAdapter_id = mypaste(adapter_id),
            adapter_id = adapter_id[hit]
            )
    } else if (hitBy =="distanceToEnd"){
        res$read_length = width(myseqs)[res$group]
        rv = res %>% group_by(group) %>% summarise(
            hit = which.min(pmin(start,read_length-end)),
            start = start[hit],
            end = end[hit],
            allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width),
            allAdapter_id = mypaste(adapter_id),
            adapter_id = adapter_id[hit]
            )
    }
    rv
}
