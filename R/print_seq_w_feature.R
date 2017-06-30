print_seq_w_feature = function(x,featureRange,featureBase="X",
    lineWidth=80, outfile,header,append=TRUE) {
    ## featurebase can be as long as featureRange to define base for each range
    mystringwrap = function(s, width=80) {
        substring(s,
            seq(1,max(c(nchar(s),width)),width),
            seq(min(c(nchar(s),width)),max(c(nchar(s),width)),width))
    }

    if(!is.character(x)) x= as.character(x)

    if(length(featureBase)==1) {
        featureBase = rep(featureBase,length(featureRange))
    }

    ## in line replacement
    rv = rep(" ",nchar(x))
    rv[as.integer(featureRange)] = rep(featureBase,width(featureRange))

    outPretty = mystringwrap(x,lineWidth)
    rvPretty = rep("\n",3*length(outPretty))
    rvPretty[seq(1,length(rvPretty),3)] = outPretty
    rvPretty[seq(2,length(rvPretty),3)] = mystringwrap(paste(rv, collapse=""),lineWidth)

    if(!missing(header)){
        rvPretty = c(header,rvPretty)
    }

    if(!missing(outfile)){
        writeLines(rvPretty,file(outfile, open = ifelse(append,"a","w")))
    }

    rvPretty
}
