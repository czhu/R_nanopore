## extract polyA or polyT
find_polyN  = function(subject, pattern="AAAAAAAAAA",max.mismatch=1L) {
    if(!is(pattern,"DNAString")) { pattern = DNAString(pattern) }

    IRanges::reduce(vmatchPattern_direct(pattern, subject,
            max.mismatch=max.mismatch, min.mismatch=0L, with.indels=TRUE,
            fixed=c(TRUE,TRUE),algo="auto"))
}
