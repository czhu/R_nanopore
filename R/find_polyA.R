find_polyA = function(thisRead,featureType="polyA",
    feature=list(
        "fwd"="AAAAAAAAAAAA",
        "rev"="TTTTTTTTTTTT"
        )){
    rv = c(
        GRanges(thisRead$qname,get_polyN(thisRead$read,feature$fwd),"+"),
        GRanges(thisRead$qname,get_polyN(thisRead$read,feature$rev),"-")
    )
    rv$featureType=featureType
    rv
}
