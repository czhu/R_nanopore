##------------------------------------------------------------
## from tilingArray by Wolfgang Huber
##------------------------------------------------------------
alongChromTicks = function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}

## extract a few plotting parameters from the data if not present, eventually we let
## change this
process_plotDat = function(plotDat){
    ## for each panel of dataset we get the height and max height, same order as Tx
    if(is.null(plotDat$param$sizeFactor)) plotDat$param$sizeFactor = rep(1, length(plotDat$param$trackHeight))
    ## normalise by sizeFactor
    plotDat$param$countMat = do.call(cbind,lapply(plotDat$data, function(x) sapply(x,length)))
    plotDat$param$normCountMat = t(t(plotDat$param$countMat)/plotDat$param$sizeFactor)
    plotDat$param$trackHeight = rowMax(plotDat$param$normCountMat)
    plotDat$param$trackNames = names(plotDat$tx)
    plotDat$param$panelNames = colnames(plotDat$param$countMat)
    plotDat$param$nTrack = length(plotDat$tx)
    plotDat$param$nPanel = ncol(plotDat$param$countMat)

    plotDat
}
validate_plotDat = function(plotDat){
    if(is.null(plotDat$param$strand) | (!plotDat$param$strand %in% c("+","-")))
        stop("strand information needed for gene plot!")
}
convert_num_to_range = function(x,from=range(x),to=c(0.2,0.8)){
    ## mapping number to another range defined by to
    ## https://stackoverflow.com/questions/345187/math-mapping-numbers
    ## Y = (X-A)/(B-A) * (D-C) + C
    x[x<from[1]]=from[1]
    x[x>from[2]]=from[2]
    (x-from[1])/diff(from) * diff(to) + to[1]
}

num_to_color = function(x,from=range(x),to=c(0.2,0.8), cols=RColorBrewer::brewer.pal(9, "Blues")) {
    ## mapping numbers to range range, typically for color gradient etc
    mycolfun = colorRamp(cols)
    rgb(mycolfun(convert_num_to_range(x,from,to)), maxColorValue=255)
}

simplify_num_output = function(x){sprintf("%.01f",x)}

### an very ancient implementation should be removed at some point
compare_plot = function(plotDat, plotTxLabel=TRUE,doCDS = TRUE,debug=FALSE, drawSpaceBetweenReads = TRUE,
    drawPanelRect = TRUE, drawReadCol = TRUE, drawReadBorder = FALSE, lineAlpha=0.2,lineWidth,
    doLine=TRUE, lineType= "dotted",config){
    ## default settings
    ## shall we draw CDS
    ## line between exons
    readAltColor = "gray40"
    if(missing(config)){
        config = default_config()
    }
     #"solid" ## line connecting reads
    ### all isoform variation
    hasMultiPanel = plotDat$param$nPanel > 1
    if(!hasMultiPanel) {
        if(!is.null(plotDat$data[[1]])){
            plotDat$name = paste(plotDat$name, names(plotDat$data[1]))
        }
    }

    ######### title extra space
    pushViewport(viewport(layout = grid.layout(2, 1, height = c(0.25,10)),width=0.95,height=0.95))
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    grid.text(plotDat$name)
    popViewport()
    ##########

    ######### calculate the layout
    ##for per isoform, one annotation, one data, data totally take 10, annotation 1 each
    dataTracksHeight = rep(0.5, plotDat$param$nTrack*2)
    spaceForData = 10
    names(dataTracksHeight)[seq(1,length(dataTracksHeight),by=2)] = paste0("annot_", seq_len(plotDat$param$nTrack))
    names(dataTracksHeight)[seq(2,length(dataTracksHeight),by=2)] = paste0("data_", seq_len(plotDat$param$nTrack))
    ## heigh based on normalised heights
    dataTracksHeight[seq(2,length(dataTracksHeight),by=2)] = spaceForData * (plotDat$param$trackHeight / sum(plotDat$param$trackHeight))
    VP = c(coord = 0.5, "GeneModel" = 0.5, dataTracksHeight)
    if(plotDat$param$strand == "+") VP = rev(VP)

    if(hasMultiPanel) {VP = c(title = 0.25, VP)}

    ######### caculate plot range
    coord = c(min(start(plotDat$geneModel)), max(end(plotDat$geneModel)))
    extraSpace = min(1000, diff(coord)*0.1)
    coord = c(coord[1]- extraSpace, coord[2] + extraSpace)
    ##########

    grl = grid.layout(1, plotDat$param$nPanel)
    pushViewport(viewport(layout = grl, layout.pos.col = 1, layout.pos.row = 2))

    for(vpCol in 1:plotDat$param$nPanel) {
        pushViewport(viewport(layout.pos.col = vpCol, layout.pos.row = 1))
        if(drawPanelRect) grid.rect(gp = gpar(lwd = 1,col="grey", alpha=0.4))

        pushViewport(viewport(layout = grid.layout(length(VP), 1, height=VP), width=0.95,height=0.95))

        ## work in points: totalSpace * ratioOfSpaceForData / numEffCounts
        ## read height normalised
        readHeight = convertHeight(unit(1,"npc"),"points",valueOnly=TRUE) *
            spaceForData/sum(VP) / sum(plotDat$param$trackHeight) / plotDat$param$sizeFactor[vpCol]

        spaceBetweenReadsInPoint= ifelse(drawSpaceBetweenReads,readHeight/8,0)

        ############ extra title for the panels
        if(hasMultiPanel) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = which(names(VP)=="title")))
            grid.text(plotDat$param$panelNames[vpCol])
            popViewport()
        }
        ############

        ############ draw coord
        plot_coord(coord, vpr=which(names(VP)=="coord"))
        ############

        ############ draw GeneModel
        plot_feature_vpr(plotDat$geneModel, vpr=which(names(VP)=="GeneModel"),
            coord = coord,
            featureHeight = 5, featureAlpha = 0.8, doLine=FALSE,
            featureCols = config$geneModel$color,
            lineAlpha=0.5, lineType= "dotted",center=TRUE)

        ############ draw data per panel
        for(txIdx in 1:plotDat$param$nTrack){
            ########## draw tx annotation
            thisFeature = plotDat$tx[txIdx]
            thisCols = if(is.null(thisFeature$itemRgb)){config$readConsensus$color}
                else { thisFeature$itemRgb}

            if(plotTxLabel){
                plot_feature_vpr(thisFeature, vpr=which(names(VP)==paste0("annot_",txIdx)),
                    coord = coord,featureHeight = 4, featureAlpha = 0.8, doLine=TRUE,
                    featureCols = thisCols,lineAlpha=0.5, lineType= "dotted",
                    center=TRUE, keepOrder=FALSE)
            } else {
                plot_feature_vpr(thisFeature, vpr=which(names(VP)==paste0("annot_",txIdx)),
                    coord = coord,featureHeight = 4, featureAlpha = 0.8, doLine=TRUE,
                    featureCols = thisCols,lineAlpha=0.5, lineType= "dotted",
                    center=TRUE, keepOrder=FALSE)
            }

            ###### draw data for isoform x
            if(plotDat$param$countMat[txIdx,vpCol]>0){

                thisData = plotDat$data[[vpCol]][[txIdx]]

                if(drawReadCol){
                    mycols = rep(num_to_color(mcols(plotDat$data[[vpCol]][[txIdx]])$qscore,from=c(5,12)),
                        length(thisData))
                } else {mycols=config$read$color}

                plot_feature_vpr(thisData,
                    vpr = which(names(VP)==paste0("data_",txIdx)), coord = coord,
                        featureHeight = readHeight,
                        featureAlpha = 0.8, doLine=TRUE, featureCols = mycols,
                        lineAlpha=0.5, lineType= "dotted", spaceBetweenFeatures=0,
                        plotBottomToTop=FALSE,center=FALSE,lineWidth=0.2,
                        textLabelFront = length(thisData) )

                # if(debug)
                #     grid.text(
                #         paste("Raw", plotDat$param$countMat[txIdx,vpCol],
                #             "Norm", simplify_num_output(plotDat$param$normCountMat[txIdx,vpCol]),
                #             "SF", simplify_num_output(plotDat$param$sizeFactor[vpCol])),
                #         y = unit(ifelse(plotDat$param$strand == "+", 0.9, 0.1), "npc"), gp=gpar(cex=0.4))
                #mypos = pretty(c(1,length(myfeature)))
                #grid.yaxis(at=convertY(unit(mypos * eachReadSpace,"points"),"npc",valueOnly=TRUE),label=mypos,gp=gpar(col="red"))
            }

        }
        popViewport()
        popViewport()
    }

    popViewport(2)
}
