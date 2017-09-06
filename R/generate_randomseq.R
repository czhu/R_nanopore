generate_randomseq = function(n,l){
    ## n number of reads, l length of the read
    mat = matrix(sample(Biostrings::DNA_BASES, l * n, replace = TRUE), nrow = n)
    do.call(paste0,lapply(1:ncol(mat),function(i) mat[,i]))
}
