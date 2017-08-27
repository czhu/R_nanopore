IRangesList_2_tibble_psl = function(x) {
     mypaste = function(y) paste(y,collapse=",")
     as_tibble(x) %>% group_by(group) %>%
     summarise(
         start = start(which.max(width)),
         end = start(which.max(width)),
         allStart=mypaste(start),allEnd=mypaste(end),allWidth=mypaste(width))
}

tibble_psl_2_IRangesList = function(x) {
    # x is a tibble or data frame containing start and end column, both are concatenated by ","
    myunlist = function(y) as.numeric(unlist(strsplit(y,",")))
    #rv = vector("list", nrow(x))
    rv = strsplit(x$start,",")
    relist(IRanges(start=myunlist(x$start), end=myunlist(x$end)),rv)
}
