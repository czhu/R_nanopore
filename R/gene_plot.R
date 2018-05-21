### for plot clustered read with consensus of annotation
### a complete rewrite 
gene_plot = function(plotDat, title, plotTxLabel=TRUE,
    debug=FALSE, drawSpaceBetweenReads = TRUE,
    drawPanelRect = TRUE, drawReadCol = TRUE, drawReadBorder = FALSE, lineAlpha=0.2,lineWidth,
    doLine=TRUE, lineType= "dotted",readHeight){
    ## default settings
    ## shall we draw CDS
    ## line between exons
    readColor = "steelblue"
    readAltColor = "gray40"
     #"solid" ## line connecting reads
    ### all isoform variation
    plotBottomToTop = FALSE
    if(missing(readHeight)) readHeight = 2

    nDataTrack = length(plotDat$consensus)

    dataTracksHeight = rep(0.5, nDataTrack * 2)
    dataPerTrackHeights = sapply(plotDat$reads,length)
    spaceForData = 8

    names(dataTracksHeight)[seq(1,length(dataTracksHeight),by=2)] = paste0("annot_", seq_len(nDataTrack))
    names(dataTracksHeight)[seq(2,length(dataTracksHeight),by=2)] = paste0("data_", seq_len(nDataTrack))
    ## heigh based on normalised heights
    dataTracksHeight[seq(2,length(dataTracksHeight),by=2)] = spaceForData * (dataPerTrackHeights / sum(dataPerTrackHeights))

    if(is.null(plotDat$geneModel)){
        VP = c(coord = 0.5, consensus = 0.7, dataTracksHeight)
    } else {
        VP = c(coord = 0.5, GeneModel = 1, consensus = 0.7, dataTracksHeight)
    }
    if(plotBottomToTop) VP = rev(VP)


    ######### caculate plot range
    coord = c(min(start(plotDat$consensus)), max(end(plotDat$consensus)))
    extraSpace = min(1000, diff(coord)*0.1)
    coord = c(coord[1]- extraSpace, coord[2] + extraSpace)
    ##########

    ######### title extra space
    grl = grid.layout(length(VP), 1, heights=VP)
    if(!missing(title)){
        pushViewport(viewport(layout = grid.layout(2, 1, height = c(0.25,10)),width=0.95,height=0.95))
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
        grid.text(title)
        popViewport()
        pushViewport(viewport(layout = grl, layout.pos.col = 1, layout.pos.row = 2))
    } else {
        pushViewport(viewport(layout = grl, width=0.95, height=0.95))
    }

    #if(drawPanelRect) grid.rect(gp = gpar(lwd = 1,col="grey", alpha=0.4))
    ############ draw coord
    plot_coord( coord, vpr = which(names(VP)=="coord"))

    ############ draw GeneModel
    if(!is.null(plotDat$geneModel)){
        plot_feature(plotDat$geneModel, vpr=which(names(VP)=="GeneModel"), coord = coord,
            featureHeight = 4, featureAlpha = 0.8, doLine=TRUE, featureCols = "gray40",
            lineAlpha=0.5, lineType= "dotted", drawSpaceBetweenReads=TRUE,center=TRUE)
    }

    #############
    thisCols = if(is.null(plotDat$consensus$itemRgb)){"darkgreen"} else {plotDat$consensus$itemRgb}
    plot_feature(plotDat$consensus, vpr=which(names(VP)=="consensus"), coord = coord,
        featureHeight = 4, featureAlpha = 0.8, doLine=TRUE, featureCols = thisCols,
        lineAlpha=0.5, lineType= "dotted", drawSpaceBetweenReads=TRUE, center=TRUE)

    #### data
    for(ithCluster in 1:(nDataTrack)){
        ########## draw tx annotation

        plot_feature(plotDat$consensus[ithCluster],
            vpr=which(names(VP)==paste0("annot_",ithCluster)), coord = coord,
            featureHeight = 4, featureAlpha = 0.5, doLine=TRUE, featureCols = "firebrick",
            lineAlpha=0.5, lineType= "dotted", drawSpaceBetweenReads=FALSE,plotBottomToTop=FALSE, center=TRUE)
        plot_feature(plotDat$reads[[ithCluster]],
            vpr=which(names(VP)==paste0("data_",ithCluster)), coord = coord,
                featureHeight = readHeight, featureAlpha = 0.8, doLine=TRUE, featureCols = "steelblue",
                lineAlpha=0.5, lineType= "dotted", drawSpaceBetweenReads=TRUE, plotBottomToTop=FALSE,center=FALSE,lineWidth=0.2)
    }
    popViewport()
}
