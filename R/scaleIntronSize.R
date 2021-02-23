## scale intron size
# tx = asBED(GRangesList(
#     c(GRanges("x:1-10"),
#     GRanges("x:100-110"),
#     GRanges("x:200-210")),
#     c(GRanges("x:1-10"),
#     GRanges("x:100-110"),
#     GRanges("x:200-210"))
# ))
#plotRanges(ranges(unlist(blocks(rv))))
# untested
scaleIntronSize = function(tx, scaleFun=log10){
    ## plotRanges(ranges(unlist(blocks(tx))))
    nexons = function(x) lengths(blocks(x))
    intrn = extractIntrons(tx)
    intronSizeAfterScalling = scaleFun(width(blocks(intrn)))
    sizeShift = rep( 0, sum(nexons(tx)) )
    ## index of first exon of each transcript
    indexFirstExonEachTranscript = c(1, cumsum(nexons(tx))[-length(tx)] +1)
    sizeShift[-indexFirstExonEachTranscript] = unlist(width(blocks(intrn))) - unlist(intronSizeAfterScalling)
    sizeShift = unlist(cumsum(relist(sizeShift, blocks(tx))))

    txTmp = unlist(blocks(tx))
    grShifted = GRanges(seqnames(txTmp), IRanges(start=start(txTmp)-sizeShift, end=end(txTmp)-sizeShift),strand=strand(txTmp))
    rv = asBED(relist(grShifted, blocks(tx)))
    return(rv)
}
