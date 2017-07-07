## the c stop working properly, i get a list
## for fastq file with less then 1M reads no problem
stream_read_fastq = function(infile, FUN=identity) {
    # readFastq doesn't work
    # message: line too long /Users/czhu/data/DCM/pore/analysis_201704/reads.fastq.gz:55069
    fs= FastqStreamer(infile)
    rv = FUN(yield(fs))

    while (length(fq <- yield(fs))) {
         ## do work here
        rv  = c(rv,FUN(fq))
     }
    close(fs)
    rv
}
