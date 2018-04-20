## from tiling array package
# this function sets up a new viewport.  It is used by plotAlongChromLegend,
# plotSegmentationHeatmap and plotSegmentationDots when they are called as
# stand-alone functions (ie when vpr is not specified)
new_vp = function(main, cexMain=1, dataPanelHeight=1, vpHeight=0.7, titleOffSet=0) {
     if(!missing(main)) {
        vpr = c("title"=0.1, "data"=dataPanelHeight)
        pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
        pushViewport(viewport(layout=grid.layout(length(vpr), 1, heights=vpr)))
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(vpr)=="title")))
        grid.text(label=main, x=0.5, y=1.1+titleOffSet, just="centre", gp=gpar(cex=cexMain))
        popViewport()
     } else {
        vpr = c("data"=dataPanelHeight)
        pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
        pushViewport(viewport(layout=grid.layout(length(vpr), 1, heights=vpr)))
     }
  return(which(names(vpr)=="data"))
}

plot_coord = function(coord, vpr) {
    if(missing(vpr)) {
        vpr = new_vp()
    }

    ############ draw coord
    pushViewport(dataViewport(xData=coord, yscale=c(0,0.3), extension=0, clip="off",
                              layout.pos.col=1, layout.pos.row=vpr))
    grid.lines(coord, c(0,0), default.units = "native")
    tck = alongChromTicks(coord)
    grid.text(label=formatC(tck, format="d"), x = tck, y = 0.1,
              just = c("centre", "bottom"), gp = gpar(cex=.5), default.units = "native")
    grid.segments(x0 = tck, x1 = tck, y0 = 0, y1 = 0.08,  default.units = "native")
    popViewport()

}

### this is core plotting function for plotting granges with exons (in blocks slot)
### FIXME: x could be GenomicRangesList, where each element in a list is a exon, this would be more general
plot_feature_vpr  = function(x, vpr, coord, lineWidth, featureCols="steelblue", featureAlpha=1, featureHeight=10,
    doLine=TRUE, lineAlpha=0.5, lineType= "dotted", plotBottomToTop  = FALSE, plotNames,
    drawSpaceBetweenReads=TRUE, center=FALSE) {
    ## x is a GRanges object with blocks
    ## conivence functon to call plot_feature with vpr
    if(missing(vpr)) {
        vpr = new_vp()
    }
    if(missing(coord)) {
        coord = c(min(start(x)), max(end(x)))
    }
    pushViewport(
        dataViewport(xData=coord, yscale=c(0,1), extension=0, clip="off",
        layout.pos.col=1,layout.pos.row=vpr))
    plot_feature(x=x, coord=coord, lineWidth=lineWidth,
            featureCols=featureCols, featureAlpha=featureAlpha, featureHeight=featureHeight,
            doLine=doLine, lineAlpha=lineAlpha, lineType= lineType, plotBottomToTop  = plotBottomToTop, plotNames,
            drawSpaceBetweenReads=drawSpaceBetweenReads, center=center)
    popViewport()
}


plot_feature  = function(x, coord, lineWidth, featureCols="steelblue", featureAlpha=1, featureHeight=10,
    doLine=TRUE, lineAlpha=0.5, lineType= "dotted", plotBottomToTop  = FALSE, plotNames,
    drawSpaceBetweenReads=TRUE, center=FALSE) {
    ## key function used to plot read and tx annotation
    ## x is a GRanges object with blocks
    ## plotBottomToTop TRUE for "+" strand FALSE for minus strand
    if(missing(coord)) {
        coord = c(min(start(x)), max(end(x)))
    }

    thisMaxHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)

    mybins = disjointBins(x, ignore.strand=TRUE)
    nfeature = max(mybins)
    marginSpace = thisMaxHeight - featureHeight * nfeature

    ### avoid read overflow i.e. drawing space cannot acommodate this many feature given the size
    featureHeight = min(featureHeight, thisMaxHeight/ nfeature)

    if(is.numeric(drawSpaceBetweenReads)){
        ## space cannot exceed featureHeight
        spaceBetweenReadsInPoint=  drawSpaceBetweenReads
        featureHeightInPoint = featureHeight
        featureHeight = featureHeightInPoint + spaceBetweenReadsInPoint
        message("Plotting with space ", drawSpaceBetweenReads)
    } else{
        spaceBetweenReadsInPoint = ifelse(drawSpaceBetweenReads,featureHeight/8,0)
        featureHeightInPoint = featureHeight - spaceBetweenReadsInPoint
    }

    myfeature = blocks(x)
    myx = unlist(start(myfeature))

    ## for - strand stack top to bottom, for + strand bottom to top
    if(plotBottomToTop){### usually for "+" strand
        yPerRead = (mybins-1) * featureHeight
        if(center) yPerRead = yPerRead + marginSpace/2
    } else{ ## usually for "-" strand
        yPerRead = thisMaxHeight - mybins * featureHeight
        if(center) yPerRead = yPerRead - marginSpace/2
    }
    nFeatureEach = lengths(myfeature)
    myy = rep(yPerRead, nFeatureEach)
    if(length(featureCols)>1){
        if(!(length(x) == length(featureCols))) stop("featureCols should have the same length as x or 1\n")
        featureCols = rep(featureCols, nFeatureEach)
    }

    grid.rect(myx, unit(myy,"points"), width=unlist(width(myfeature)),
        height=unit(featureHeightInPoint,"points"), gp=gpar(col = NA , fill = featureCols, alpha=featureAlpha),
        default.units="native", just=c("left","bottom"))

    cumLength = cumsum(elementNROWS(myfeature))
    myxStart = unlist(end(myfeature))[-cumLength]
    #lapply(end(myfeature),function(x) x[-length(x)])
    myxEnd = unlist(start(myfeature))[-c(1,cumLength[-length(cumLength)]+1)]
    ## lapply(start(myfeature),function(x) x[-1])
    myyLine = c(rep(yPerRead, nFeatureEach-1), rep(yPerRead, nFeatureEach-1))

    if(doLine & length(c(myxStart,myxEnd))>0){
        #penaltyFactorReadNumber = (1/log10(plotDat$param$normCountMat[txIdx,vpCol]))^2
        grid.polyline(
            x=unlist(c(myxStart,myxEnd)), y=unit(myyLine+ featureHeightInPoint/2,"points"),
            id = rep(1:length(unlist(myxStart)),2),
            gp=gpar(col=featureCols,
                lwd=if(missing(lineWidth)) unit(min(1,featureHeight/10),"points") else {
                unit(lineWidth,"points")},
            ## FIXME scale alpha depending on the number of reads
                alpha=lineAlpha,lty=lineType), ##lex=1/penaltyFactorReadNumber),
            default.units = "native")
    }
    if(!missing(plotNames)){
        nTrack = max(mybins)
        if(plotBottomToTop) {
            thisy = mybins * (1/nTrack) + 1/nTrack/2
        } else {
            thisy = 1 - mybins * (1/nTrack)  + 1/nTrack/2
        }

        grid.text(plotNames,
            x = unit((start(x) + end(x))/2,"native"),
            y = unit(thisy, "npc"),gp=gpar(cex=0.4))
    }

}
