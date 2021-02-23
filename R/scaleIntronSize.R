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
scaleIntronSize = function(tx, scaleFun=log10){
    ## plotRanges(ranges(unlist(blocks(tx))))
    intrn = extractIntrons(tx)
    intronSizeAfterScalling = scaleFun(width(blocks(intrn)))
    sizeShift = rep(0, sum(lengths(blocks(tx))))
    sizeShift[-c(1,lengths(blocks(tx))[-length(tx)] +1)] = unlist(width(blocks(intrn))) - unlist(intronSizeAfterScalling)

    sizeShift = unlist(cumsum(relist(sizeShift, blocks(tx))))

    txTmp = unlist(blocks(tx))
    grShifted = GRanges(seqnames(txTmp), IRanges(start=start(txTmp)-sizeShift, end=end(txTmp)-sizeShift),strand=strand(txTmp))
    rv = asBED(relist(grShifted, blocks(tx)))
    return(rv)
}
