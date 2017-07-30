## extract polyA or polyT
find_polyN = function(x,polyN=paste(rep("A",10),collapse="")){
    GenomicRanges::reduce(ranges(matchPattern(polyN,x,max.mismatch=1, with.indels =TRUE)))
}
