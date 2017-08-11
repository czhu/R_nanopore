vmatchPattern_direct <- function(pattern, subject,
                              max.mismatch=1L, min.mismatch=0L,
                              with.indels=FALSE, fixed=c(TRUE,TRUE),
                              algo="auto")
   {   ##https://support.bioconductor.org/p/58350/
       ## pattern must be DNAString
       ## subject must be DNAStringSet
     require(Biostrings)
     algo <- Biostrings:::selectAlgo(algo, pattern,
                                     max.mismatch, min.mismatch,
                                     with.indels, fixed)
     C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
                     max.mismatch, min.mismatch,
                     with.indels, fixed, algo,
                     "MATCHES_AS_RANGES",
                     PACKAGE="Biostrings")
     unlisted_ans <- IRanges(start=unlist(C_ans[[1L]],
         use.names=FALSE), width=unlist(C_ans[[2L]],use.names=FALSE))
     relist(unlisted_ans, C_ans[[1L]])
   }
