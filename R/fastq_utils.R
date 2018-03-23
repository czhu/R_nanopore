## a collection of functions to manupilate fastq inputs
## supplement ShortRead package

change_header = function(x, header){
    ## x is a ShortReadQ class
    ## this is essentially the missing id<- method
    ShortReadQ(sread=sread(x), quality = quality(x), id = BStringSet(header))
}

get_seqid = function(x, header){
    ## x is a ShortReadQ class
    sub("^(\\S+).*", "\\1", ShortRead::id(x))
}
