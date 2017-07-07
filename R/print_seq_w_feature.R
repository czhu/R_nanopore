print_seq_w_feature = function(x,featureRange,featureBase="X",
    lineWidth=80, outfile,header,append=TRUE) {
    ## featurebase can be as long as featureRange to define base for each range
    mystringwrap = function(s, width=80) {
        substring(s,
            seq(1,max(c(nchar(s),width)),width),
            seq(min(c(nchar(s),width)),max(c(nchar(s),width)),width))
    }
    myrange = IRanges(start(featureRange), end(featureRange))
    bins <- disjointBins(myrange)

    if(!is.character(x)) x= as.character(x)

    if(length(featureBase)==1) {
        featureBase = rep(featureBase,length(featureRange))
    }

    outPretty = mystringwrap(x,lineWidth)
    nLinePerEntry = 1+max(bins)
    rvPretty = rep("\n", nLinePerEntry * length(outPretty))

    rvPretty[seq(1,length(rvPretty),nLinePerEntry)] = outPretty
    for(thisBin in unique(bins)){
        ## in line replacement
        wh = bins==thisBin
        rv = rep(" ",nchar(x))
        rv[as.integer(myrange[wh])] = rep(featureBase[wh],width(myrange[wh]))

        rvPretty[seq(thisBin+1,length(rvPretty),nLinePerEntry)] = mystringwrap(paste(rv, collapse=""),lineWidth)
    }
    #rvPretty = rvPretty[rvPretty!=paste(rep(" ",lineWidth),collapse="")]

    if(!missing(header)){
        rvPretty = c(header,rvPretty)
    }

    if(!missing(outfile)){
        mycon = file(outfile, open = ifelse(append,"a","w"))
        writeLines(rvPretty,mycon)
        close(mycon)
    }

    rvPretty
}
