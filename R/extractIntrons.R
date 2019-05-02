## given multi exonic tx to extract the ranges for their introns
extractIntrons = function(x){
    ## x is Granges with blokcs to exons
    ## this is a vectoriesd implementation
    ## this should give the same results as clname_to_iv
    if(!all(lengths(blocks(x))))
        stop("extractIntrons only works with multi exonic tx")
    tmpVal = blocks(x)
    nexons = lengths(tmpVal)
    allstarts = unlist(start(tmpVal))
    idxShift = cumsum(nexons)
    allstarts = allstarts[-c(1, (idxShift+1)[-length(idxShift)] )]
    allends = unlist(end(tmpVal))
    allends = allends[ -idxShift]
    myfac = rep(1:length(x), nexons-1)
    ans = GRanges(
        rep(as.character(seqnames(x)), nexons-1),
        IRanges(start=allends+1, end=allstarts-1),
        strand=rep(as.character(strand(x)), nexons-1)
    )
    ans = asBED(split(ans, factor(myfac,1:length(x))))
    if( !is.null(names(blocks(x))) ){
        names(ans) = names(blocks(x))
    }
    return(ans)
}
