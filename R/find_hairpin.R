## a convience function to map hairpin
find_hairpin = function(thisRead,featureType="hairpin",
    feature=list(
        "fwd"="GTCGTATCCAGTGCAGGGTCCGAGGTATTCGCACTGGATACGAC",
        "rev"="GTCGTATCCAGTGCGAATACCTCGGACCCTGCACTGGATACGAC"
        ),tol=0.2){
    rv = matchPatternStranded(feature$fwd, feature$rev, thisRead$read, thisRead$qname, tol=tol)
    rv$featureType=featureType
    rv
}
