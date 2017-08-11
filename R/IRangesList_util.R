IRangesList_2_tibble_psl = function(x) {
     mypaste = function(y) paste(y,collapse=",")
     as_tibble(x) %>% group_by(group) %>%
     summarise(
         start = start(which.max(width)),
         end = start(which.max(width)),
         allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width))
}

tibble_psl_2_IRangesList = function(x, listLength = max(x$group)) {
    myunlist = function(y) as.numeric(unlist(strsplit(y,",")))

    rv = vector("list", listLength)
    rv[x$group] = strsplit(x$start,",")

    relist(IRanges(start=myunlist(x$start), end=myunlist(x$end)),rv)
}
