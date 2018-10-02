## this is a complete rewrite of the existing gene_plot function as of 2018-09-18
## gene_plot_x is the a fexible function for plotting grouped reads in a gene view
## normally all reads will be plotted (natural scale)
## with additional wrapper gene_plot_x_log (to be developed) log level view is possible
## with only representative reads, in which selected number of reads is precomputed
## or with the wrapper gene_plot_x_multiplanel (to be developed) multiple samples
## can be compared for the same gene, in which case trackHeight for multiple samples
## is precomputed

## the input is a list of two elements
## plotData (list):
##      -- GeneModel (txRanges)
##      -- GroupedReads (list)
##          -- Gread1
##              -- consensus (txRanges)
##              -- reads (txRanges)
##          -- Gread2
##              -- consensus (txRanges)
##              -- reads (txRanges)
##           .
##           .
##           .
## for each txRanges we allow the following paremeter in the config in addition
## to trackHeight
## featureColor, featureHeight, spaceBetweenFeatures
## optionnal: highlight (txRanges with label slot)
##
## config (list):
##      -- Axis (list, if NULL not plotted)
##          -- trackHeight
##      -- GeneModel
##          -- trackHeight
##          -- featureColor
##          ..
##      -- GroupedReads
##          -- Gread 1
##          -- Gread 2
##
## additional global controls plotFromToBottom

## forget about all the above
## we just need two slots for the underlying plots
## data (txRanges)
## config

## test purpose
## single panel
# data(immt)
# gene_plot_single_panel(GeneModel=geneModelFull, txConsensus, txReads, txHighlight,title ="Mytest")
# gene_plot_single_panel(GeneModel=geneModelProtLinc, txConsensus, txReads, txHighlight)

## multi panel
# txReadsList = list(txReads, txReads, txReads)
# names(txReadsList) = c("WC","WC=NC","R634Q")
# sizeFactor = c(2,1,0.5)
# gene_plot_multi_panel(GeneModel=asBED(GRangesList(disjoin(geneModelProtLinc))),
#     txConsensus, txReadsList, txHighlight,
#     title="IMMT", plotTopToBottom=TRUE, axisHeight= 15, consensusTrackHeight = 5,
#     consensusFeatureHeight = 4, sizeFactor=sizeFactor)

plot_config_default = function(...){
    list(
        featureHeight = 2,
        featureAlpha = 0.75,
        doLine = TRUE,
        featureColor = "steelblue",
        lineAlpha = 0.5,
        lineType = "dotted",
        spaceBetweenFeatures=0,
        plotBottomToTop = FALSE,
        center = FALSE,
        lineWidth = 0.2,
        textLabelFront = NULL,
        highlightFontsize = 6,
        highlightColor = "black"
    )
}

add_title = function (title, titleHeight=10,titleFontSize = 7) {
    totalHeightInPoints = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE)
    message(totalHeightInPoints)
    pushViewport(
        viewport(
            layout = grid.layout(2, 1,
            height = unit(c(titleHeight,totalHeightInPoints-titleHeight),"points")),
        width=1,height=1))
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    grid.text(title, gp = gpar(fontsize=titleFontSize))
    popViewport()
}

GENEMODELFEATUREHEIGHT = 3
GENEMODELSPACEBETWEENFEATURES = 1
GENEMODELEXTRAMARGIN = 5

gene_plot_multi_panel = function(GeneModel, txConsensus, txReadsList, txHighlight,
    title, plotTopToBottom=TRUE, axisHeight= 15, consensusTrackHeight = 5,
    consensusFeatureHeight = 4, sizeFactor=rep(1, length(txReadsList)) ){

    if(!missing(title)){
        add_title(title, titleHeight=10)
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    } else {
        pushViewport(viewport(width=1, height=1))
    }

    ## all we need to do is figure out the track height
    ## and call gene_plot_single_panel with defined trackHeights

    geneModelTrackHeight = length(GeneModel) *
        (GENEMODELFEATUREHEIGHT + GENEMODELSPACEBETWEENFEATURES) + GENEMODELEXTRAMARGIN

    totalHeightInPoints = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE)
    trackHeightsConsensus = rep(consensusTrackHeight, length(txConsensus))

    spaceLeftForReads = totalHeightInPoints - axisHeight -
        sum(trackHeightsConsensus) - geneModelTrackHeight

    # if(spaceLeftForReads<100){
    #     stop("not enough space for plotting reads, considering inrease page height")
    # }
    ## tx by sample
    countMat = do.call(cbind, lapply(txReadsList, function(x) sapply(x, length)))
    countMatNorm = countMat / sizeFactor
    nMaxReadsPerTrack = apply(countMatNorm, 1, max)

    readHeightNorm = spaceLeftForReads/sum(nMaxReadsPerTrack)
    readHeightPerSample = readHeightNorm / sizeFactor

    trackHeightData = nMaxReadsPerTrack * readHeightNorm
    trackHeights =  c(geneModelTrackHeight, sapply(seq_len(length(txConsensus) * 2), function(i){
            if(i%%2){
                trackHeightsConsensus[i%/%2+1]
            } else {
                trackHeightData[i%/%2]
            }
        }))

    ### calling
    pushViewport(
        viewport( layout = grid.layout(nrow=1, ncol=length(txReadsList)) )
    )

    for( i in seq_len(length(txReadsList)) ) {
        pushViewport(viewport(layout.pos.col = i, layout.pos.row = 1))
        gene_plot_single_panel(
            GeneModel=GeneModel, txConsensus=txConsensus,
            txReadsList[[i]], txHighlight, title=names(txReadsList)[i],
            plotTopToBottom=TRUE, axisHeight= 15, consensusTrackHeight = 5,
            consensusFeatureHeight = 4,trackHeights=trackHeights,
            readHeight=readHeightPerSample[i])
        popViewport()
    }
    popViewport()

    if(!missing(title))
        popViewport()

    popViewport()
}

gene_plot_single_panel = function(GeneModel, txConsensus, txReads, txHighlight,
    title, plotTopToBottom=TRUE, axisHeight= 15, consensusTrackHeight = 5,
    consensusFeatureHeight = 4,trackHeights, readHeight, extraSpacingConsensus = 0 ) {

    if(!missing(title)){
        add_title(title, titleHeight=10)
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    } else {
        pushViewport(viewport(width=1, height=1))
    }

    ## Data tracks including Gene Model
    plotData = lapply(seq_len(length(txConsensus) * 2), function(i){
            if(i%%2){
                txConsensus[i%/%2+1]
            } else {
                txReads[[i%/%2]]
            }
        })
    # tmpdata = lapply( seq_len(length(txConsensus) *3), function(xx) as.numeric(NA) )
    # tmpdata[-seq(1, length(txConsensus) *3, by = 3)] = plotData

    plotData = c(list(GeneModel), plotData)

    ## config for Gene Model track
    geneModelConfig = plot_config_default()
    geneModelConfig$featureHeight = GENEMODELFEATUREHEIGHT
    geneModelConfig$spaceBetweenFeatures = GENEMODELSPACEBETWEENFEATURES
    geneModelConfig$featureColor = "gray40"
    geneModelConfig$center = TRUE
    geneModelConfig$plotBottomToTop = !plotTopToBottom

    ## default unit is in points
    if(missing(trackHeights)){
        geneModelTrackHeight = length(GeneModel) *
            (geneModelConfig$featureHeight + geneModelConfig$spaceBetweenFeatures) + 5

        totalHeightInPoints = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE)
        trackHeightsConsensus = rep(consensusTrackHeight, length(txConsensus))

        spaceLeftForReads = totalHeightInPoints - axisHeight -
            sum(trackHeightsConsensus) - geneModelTrackHeight - sum(rep(extraSpacingConsensus, length(txConsensus)))

        nreadsPerTrack = sapply(txReads, length)
        if(missing(readHeight)){
            readHeight = spaceLeftForReads/sum(nreadsPerTrack)
        }

        trackHeightData = nreadsPerTrack * readHeight
        trackHeights =  sapply(seq_len(length(txConsensus) * 2), function(i){
                if(i%%2){
                    trackHeightsConsensus[i%/%2+1]
                } else {
                    trackHeightData[i%/%2]
                }
            })

        trackHeights = c(geneModelTrackHeight, trackHeights)
    }

    plotConfig = lapply(seq_len(length(txConsensus) * 2), function(i){
            ans = plot_config_default()
            if(i%%2){
                ## for consensus tracks
                ans$featureHeight = 4
                ans$featureColor = txConsensus[i%/%2+1]$itemRgb
                ans$center = TRUE
                ans$textLabelFront = txHighlight[[i%/%2 + 1]]$shape
                thisTxName = granges(txConsensus[i%/%2+1])
                thisTxName$label = txConsensus[i%/%2+1]$name
                ans$highlight = thisTxName
                ans$plotBottomToTop = FALSE

                if(!is.null(txHighlight[[i%/%2 + 1]]$highlight)){
                    thisHighlight = granges(txHighlight[[i%/%2 + 1]]$highlight)
                    thisHighlight$label = txHighlight[[i%/%2 + 1]]$highlight$shape
                    ans$highlight = c(ans$highlight,thisHighlight)
                }
            } else {
                ## data tracks
                ans$featureHeight = readHeight
                ans$textLabelFront = length(txReads[[i%/%2]]) ## n reads
            }
            ans$plotBottomToTop = !plotTopToBottom
            ans
        })
    plotConfig = c(list(geneModelConfig), plotConfig)

    thisSpacing = rep(0, length(plotData))
    thisSpacing[seq(1,length(thisSpacing), by=2)] = extraSpacingConsensus

    txClass_plotter(plotData, plotConfig, axisHeight=axisHeight,
        trackHeights = trackHeights, extraSpacing = thisSpacing)
    if(!missing(title))
        popViewport()
    popViewport()
}

txClass_plotter = function(plotData, plotConfig, plotTopToBottom=TRUE,
    trackHeights, axisHeight, extraSpacing) {

    nDataTracks = length(plotData)

    if( nDataTracks != length( plotConfig) ){
        message(nDataTracks, " vs ", length(plotConfig))
        stop("Data and config and must have same number of elements")
    }


    defaultUnit = "points"
    extraMargin = 1000 ## in bp
    dataTrackPrefix = "txClass_"

    ## trackHeights
    #VP = c(axisHeight, sapply(plotConfig, "[[", "trackHeight") )

    if(missing(extraSpacing)){
        extraSpacing = rep(0, length(nDataTracks))
    }

    if( length(extraSpacing) != nDataTracks ){
        stop("extraSpacing needs to have the same length as plotData")
    }
    newTrackHeights = rep(0, nDataTracks*2)
    newTrackHeights[seq(1, length(newTrackHeights), by=2)] = trackHeights
    newTrackHeights[seq(2, length(newTrackHeights), by=2)] = extraSpacing

    names(newTrackHeights) = rep("spacing", length(newTrackHeights)/2)
    names(newTrackHeights)[seq(1, length(newTrackHeights), by=2)] = paste0(dataTrackPrefix, seq_len( nDataTracks ))

    VP = c(axisHeight, newTrackHeights)
    names(VP) = c("Axis",  names(newTrackHeights))

    if( !plotTopToBottom ) VP = rev(VP)

    grl = grid.layout(length(VP), 1, heights=unit(VP, defaultUnit))

    pushViewport(viewport(layout = grl, width=1, height=1))

    minStart = min(unlist(lapply(plotData, start)))
    maxEnd = max(unlist(lapply(plotData, end)))

    coord = c(minStart - extraMargin, maxEnd + extraMargin )
    coord = c(coord[1] - diff(coord) * 0.1, coord[2] + diff(coord) * 0.1)

    plot_coord( coord, vpr = which( names(VP) == "Axis" ) )

    for(i in seq_len( nDataTracks ) ){
        vprIndex = which(names(VP)==paste0(dataTrackPrefix, i))
        thisPlotConfig = plotConfig[[i]]
        thisPlotData = plotData[[i]]
        if( length(thisPlotData) > 0 ) {
            get_value = function(x) {
                if(is.null(thisPlotConfig[[x]]))
                    plot_config_default()[[x]]
                else
                    thisPlotConfig[[x]]
               }

            plot_feature_vpr(
                thisPlotData,
                vpr = vprIndex,
                coord = coord,
                featureHeight = get_value("featureHeight"),
                featureAlpha = get_value("featureAlpha"),
                doLine = get_value("doLine"),
                featureCols = get_value("featureColor"),
                lineAlpha = get_value("lineAlpha"),
                lineType = get_value("dotted"),
                spaceBetweenFeatures = get_value("spaceBetweenFeatures"),
                plotBottomToTop = get_value("plotBottomToTop"),
                center = get_value("center"),
                lineWidth = get_value("lineWidth"),
                textLabelFront = get_value("textLabelFront")
            )
            ## highlight
            plot_label = function(x) {
                plot_feature_text_vpr(
                    thisPlotConfig[[x]],
                    thisPlotConfig[[x]]$label,
                    vpr = vprIndex,
                    coord = coord,
                    fontsize = get_value("highlightFontsize"),
                    side = 2,
                    col = get_value("highlightColor"),
                    plotBottomToTop = !plotTopToBottom )
            }

            if( !is.null(thisPlotConfig[["highlight"]]) ){
                plot_label("highlight")
            }
        }
    }
    popViewport()
}
