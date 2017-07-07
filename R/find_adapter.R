## a convience function to map adapter
## feature alternative is acually a partial string of feature
## thanks to TSO
find_adapter = function(thisRead,featureType="adapter",
    feature=list(
        "fwd"="TATAGCGGCCGCAAGCAGTGGTATCAACGCAGAGTAC",
        "rev"="GTACTCTGCGTTGATACCACTGCTTGCGGCCGCTATA"
        ),
    featureAlternative=list(
            "fwd"="AAGCAGTGGTATCAACGCAGAGTAC",
            "rev"="GTACTCTGCGTTGATACCACTGCTT"
        ),tol=0.2){

    rv = matchPatternStranded(feature$fwd, feature$rev, thisRead$read, thisRead$qname, tol=tol)
    if(!is.null(featureAlternative)){
        rv2 = matchPatternStranded(featureAlternative$fwd, featureAlternative$rev, thisRead$read,
            thisRead$qname, tol=tol)
        rv = c(rv2[!rv2 %within% rv],rv)
    }

    rv$featureType=rep(featureType,length(rv))
    rv
}
