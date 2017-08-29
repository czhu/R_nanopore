## mutate sequencing with given error probablity
## limitation, only 1 nt insertion allowed
## vectorised version could be done using Biostrings::replaceAt
## https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/Ten_things_slides.pdf
mutate_seq = function(x, substitutionRate= 0.07,insertionRate=0.03,deletionRate=0.05) {
    require(Biostrings)

    chars <- strsplit(as.character(x), "")[[1]]
    nChars = length(chars)

    rv = ifelse(rbinom(nChars,1,substitutionRate+deletionRate),
       sample(c(DNA_BASES,""), nChars, replace = TRUE,
            prob=c(rep(substitutionRate/4,4),deletionRate)), chars)
    ## n + 1 position can be inserted
    toInsert = sample(c("",DNA_BASES), nChars+1, replace=TRUE,
        prob=c(1-insertionRate,rep(insertionRate/4,4)))
    paste(apply(cbind(c("",rv),toInsert),1,paste0), collapse="")
}
