matchPatternStranded = function(query, revQuery, subject, qname, tol=0.2,with.indels=TRUE,...){
    require(Biostrings)
    # query is string, subject is DNAString or string,
    # qname is the read name
    fwd = matchPatternWithTol(query,subject,tol, with.indels =with.indels)
    rev = matchPatternWithTol(revQuery,subject,tol, with.indels =with.indels)

    return(c(
        GRanges(rep(qname,length(fwd)),fwd,rep("+",length(fwd))),
        GRanges(rep(qname,length(rev)),rev,rep("-",length(rev)))
        )
        )

    # if(length(fwd) ==0 & length(rev) ==0){
    #     return(NULL)
    # } else if (length(fwd) ==0) {
    #     return(GRanges(qname,rev,"-"))
    # } else if (length(rev) ==0) {
    #     return(GRanges(qname,fwd,"+"))
    # } else {
    #     return(c(GRanges(qname,fwd,"+"), GRanges(qname,rev,"-")))
    # }
}
