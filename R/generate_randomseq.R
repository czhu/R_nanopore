generate_randomseq = function(n,l){
    require(Biostrings)
    apply(matrix(sample(DNA_BASES, l * n, replace = TRUE),
                   nrow = n), 1, paste, collapse = "")
}
