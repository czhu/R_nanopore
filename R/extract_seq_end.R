extract_seq_end = function(x,windowSize=200,fromEnd=FALSE){
    ## x is ShortRead class
    if(fromEnd) {
        mystart = ifelse(width(x) - windowSize +1  <1, 1, width(x) - windowSize + 1)
        myend = width(x)
    } else {
        mystart = 1
        myend = ifelse(windowSize>width(x),width(x),windowSize)
    }
    narrow(x,mystart,myend)
}
