matchPatternStranded = function(query, revQuery, subject, qname, tol=0.2,with.indels=TRUE,...){
    require(Biostrings)
    # query is string, subject is DNAString or string,
    # qname is the read name

    c(GRanges(qname,
        matchPatternWithTol(query,subject,tol, with.indels =with.indels),"+"),
    GRanges(qname,
        matchPatternWithTol(revQuery,subject,tol, with.indels =with.indels),"-")
            )
}
