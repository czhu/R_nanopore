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
# data(immt_tx_reads)
# gene_plot_single_panel(GeneModel=geneModelFull, txConsensus, txReads, txHighlight,title ="Mytest")
# gene_plot_single_panel(GeneModel=geneModelProtLinc, txConsensus, txReads, txHighlight)

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
        highlightFontsize = 4,
        highlightColor = "black"
    )
}

add_title = function (title, titleHeight=10,titleFontSize = 7) {
    totalHeightInPoints = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE)
    pushViewport(
        viewport(
            layout = grid.layout(2, 1,
            height = unit(c(titleHeight,totalHeightInPoints-titleHeight),"points")),
        width=0.95,height=0.95))
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    grid.text(title, gp = gpar(fontsize=titleFontSize))
    popViewport()
}

gene_plot_single_panel = function(GeneModel, txConsensus, txReads, txHighlight,
    title, plotTopToBottom=TRUE) {

    if(!missing(title)){
        add_title(title, titleHeight=10)
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    } else {
        pushViewport(viewport(width=0.95, height=0.95))
    }

    ## Data tracks including Gene Model
    plotData = lapply(seq_len(length(txConsensus) * 2), function(i){
            if(i%%2){
                txConsensus[i%/%2+1]
            } else {
                txReads[[i%/%2]]
            }
        })
    plotData = c(list(GeneModel), plotData)

    ## config for Gene Model track
    geneModelConfig = plot_config_default()
    geneModelConfig$featureHeight = 3
    geneModelConfig$spaceBetweenFeatures = 1
    geneModelConfig$featureColor = "gray40"
    geneModelConfig$center = TRUE
    geneModelConfig$plotBottomToTop = !plotTopToBottom

    ## default unit is in points
    axisHeight= 15
    consensusTrackHeight = 5
    consensusFeatureHeight = 4
    geneModelTrackHeight = length(GeneModel) *
        (geneModelConfig$featureHeight + geneModelConfig$spaceBetweenFeatures) + 5

    totalHeightInPoints = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE)
    trackHeightsConsensus = rep(consensusTrackHeight, length(txConsensus))

    spaceLeftForReads = totalHeightInPoints - axisHeight -
        sum(trackHeightsConsensus) - geneModelTrackHeight

    # if(spaceLeftForReads<100){
    #     stop("not enough space for plotting reads, considering inrease page height")
    # }
    nreadsPerTrack = sapply(txReads, length)
    readHeight = spaceLeftForReads/sum(nreadsPerTrack)

    trackHeightData = nreadsPerTrack * readHeight
    trackHeights =  c(geneModelTrackHeight, sapply(seq_len(length(txConsensus) * 2), function(i){
            if(i%%2){
                trackHeightsConsensus[i%/%2+1]
            } else {
                trackHeightData[i%/%2]
            }
        }))

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

    txClass_plotter(plotData, plotConfig, axisHeight=axisHeight,
        trackHeights = trackHeights)
    popViewport()
}


txClass_plotter = function(plotData, plotConfig, plotTopToBottom=TRUE,
    trackHeights, axisHeight) {
    if( length( plotData ) != length( plotConfig) )
        stop("Data and config and must have same number of elements")

    nDataTracks = length(plotData)

    defaultUnit = "points"
    extraMargin = 1000 ## in bp
    dataTrackPrefix = "txClass_"

    ## trackHeights
    #VP = c(axisHeight, sapply(plotConfig, "[[", "trackHeight") )
    VP = c(axisHeight, trackHeights)

    names(VP) = c("Axis",  paste0(dataTrackPrefix, seq_len( length(plotData) )) )

    if( !plotTopToBottom ) VP = rev(VP)

    grl = grid.layout(length(VP), 1, heights=unit(VP, defaultUnit))

    pushViewport(viewport(layout = grl, width=1, height=1))

    minStart = min(unlist(lapply(plotData, start)))
    maxEnd = max(unlist(lapply(plotData, end)))

    coord = c(minStart - extraMargin, maxEnd + extraMargin )
    coord = c(coord[1] - diff(coord) * 0.1, coord[2] + diff(coord) * 0.1)

    plot_coord( coord, vpr = which( names(VP) == "Axis" ) )

    for(i in seq_len( nDataTracks )){
        vprIndex = which(names(VP)==paste0(dataTrackPrefix, i))
        thisPlotConfig = plotConfig[[i]]
        thisPlotData = plotData[[i]]

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
            plotBottomToTop = !plotTopToBottom,
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
                side = 0,
                col = get_value("highlightColor"),
                plotBottomToTop = !plotTopToBottom )
        }

        if( !is.null(thisPlotConfig[["highlight"]]) ){
            plot_label("highlight")
        }
    }
    popViewport()
}