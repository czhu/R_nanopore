## extract polyA or polyT
get_polyN = function(x,polyN=paste(rep("A",10),collapse="")){
    reduce(ranges(matchPattern(polyN,x,max.mismatch=1, with.indels =TRUE)))
}

