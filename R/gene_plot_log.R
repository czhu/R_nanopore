## gene view with log scale representation of the supporting reads
## it's a convient wrapper to call chrom_plot

## NOTE the interface is from gene_plot in the hope that the function will
## be refact

gene_plot_log = function(GeneModel, txConsensus, txReads, txHighlight,
    title, plotTopToBottom=TRUE, axisHeight= 20, annotationTrackHeight=2,
    consensusTrackHeight = 5,consensusFeatureHeight = 3, trackHeights, readHeight) {

    myfunc = identity
    isPlusStrand = as.logical(unique(strand(txConsensus)) == "+")

    if( isPlusStrand )  myfunc = invertStrand

    plotDat = list(
        consensus = myfunc(txConsensus),
        geneModel = myfunc(GeneModel),
        highlight = txHighlight,
        reads = myfunc(txReads)
    )
    if(!missing(title)){
        plotDat$name = title
    }

    extraMargin = 1000 ## in bp
    allTx = GenomicRanges::reduce(c(granges(GeneModel), granges(txConsensus)))
    coord = allTx + extraMargin + width(allTx) * 0.1
    ## TODO change gene mode track height
    chrom_plot(plotDat,coord=c(start(coord), end(coord)),
        plotCountNum=TRUE,debug=FALSE, featureHeightPerRead = consensusFeatureHeight,
        config=generate_plot_config(list(readConsensus=list(height=consensusTrackHeight), geneModel=list(height=5))), annotationTrackHeight = annotationTrackHeight,
        spaceBetweenCluster = 4, shiftLabel=TRUE,geneNameTrackHeight=0, geneTrackHeight=0, singleStrand=TRUE, consensusNameFontsize = 5, highlightFontsize=5, genomeAxisHeight=axisHeight)
}
