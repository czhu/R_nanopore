find_polyA = function(thisRead,featureType="polyA",
    feature=list(
        "fwd"="AAAAAAAAAA",
        "rev"="TTTTTTTTTT"
        )){
    fwd = find_polyN(thisRead$read, feature$fwd)
    rev = find_polyN(thisRead$read, feature$rev)

    # rv = c(
    #     GRanges(rep(thisRead$qname,length(fwd)),fwd,rep("+",length(fwd))),
    #     GRanges(rep(thisRead$qname,length(rev)),rev,rep("-",length(rev)))
    # )
    rv = GRanges(
        c(rep(thisRead$qname,length(fwd)),rep(thisRead$qname,length(rev))),
        c(fwd,rev),
        c(rep("+",length(fwd)),rep("-",length(rev)))
        )
        #     GRanges(rep(thisRead$qname,length(rev)),rev,rep("-",length(rev)))

    rv$featureType=rep(featureType,length(rv))
    rv
}
