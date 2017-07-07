## x is a single ShortRead instance
# library(ShortRead)
# library(nanopore)
# infile = "mutant.fastq.gz"
# fq = readFastq(infile)
# x = fq[2]
# verbose=TRUE
# trim=FALSE
## this is a unfinished function
trim_adapter= function(x,verbose=FALSE,trim=TRUE){

    thisRead = list(
        read = as.character(sread(x[1])) ,
        qname = strsplit(as.character(id(x[1]))," ")[[1]][1]
        )
    if(verbose) cat("Processing",thisRead$qname,"\n")

    myfeatures = c(find_adapter(thisRead), find_polyA(thisRead))
    myfeatures
}
