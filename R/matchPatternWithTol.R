#
matchPatternWithTol = function(query,subject,tol=0.2,with.indels=TRUE,...) {
    # if 0 < tol(erance) <1, frac of query length that's allowed for max.mismatch
    if(tol>0 & tol<1){tol = round(nchar(query) * tol)}
    ranges(matchPattern(query,subject,max.mismatch=tol, with.indels =with.indels))
}
