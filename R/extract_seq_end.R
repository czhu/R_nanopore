extract_seq_end = function(x,windowSize=200,fromEnd=FALSE){
    ## x is ShortRead class
    mystart = ifelse(fromEnd,
            ifelse(width(x) - windowSize +1  <1, 1, width(x) - windowSize + 1),
            1)
    myend = ifelse(fromEnd,
        width(x), ifelse(windowSize>width(x),width(x),windowSize))
    narrow(x,mystart,myend)
}
