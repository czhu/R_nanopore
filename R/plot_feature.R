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
    spaceBetweenFeatures, center=FALSE,textLabelFront,keepOrder=FALSE) {
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
            doLine=doLine, lineAlpha=lineAlpha, lineType= lineType, plotBottomToTop  = plotBottomToTop,
            plotNames=plotNames,spaceBetweenFeatures=spaceBetweenFeatures, center=center,keepOrder=keepOrder)
    if(!missing(textLabelFront)){
        s = as.character(textLabelFront)
        extendLeft=500
        grid.text(s,
            x= unit(convertX(unit(min(start(x)) - extendLeft,"native"),"npc",
                valueOnly=TRUE)-convertX(unit(1,"strwidth","s"),"npc",valueOnly=TRUE),"npc"),
            0.5,
        just=c("left","center"),gp=gpar(fontsize=6
            #convertY(unit(0.5,"npc"),"points",valueOnly=TRUE)
        ))
    }
    popViewport()
}


plot_feature  = function(x, coord, lineWidth, featureCols="steelblue", featureAlpha=1, featureHeight=10,
    doLine=TRUE, lineAlpha=0.5, lineType= "dotted", plotBottomToTop  = FALSE, plotNames,
    spaceBetweenFeatures, center=FALSE,keepOrder=FALSE) {
    ## key function used to plot read and tx annotation
    ## x is a GRanges object with blocks
    ## plotBottomToTop TRUE for "+" strand FALSE for minus strand
    if(missing(coord)) {
        coord = c(min(start(x)), max(end(x)))
    }

    thisMaxHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)

    if(keepOrder){
        mybins = seq_len(length(x))
    } else{
        mybins = disjointBins(x, ignore.strand=TRUE)
    }

    nfeature = max(mybins)
    marginSpace = thisMaxHeight - featureHeight * nfeature

    ### avoid read overflow i.e. drawing space cannot acommodate this many feature given the size
    featureHeight = min(featureHeight, thisMaxHeight/ nfeature)

    if(!missing(spaceBetweenFeatures)){
        ## space cannot exceed featureHeight
        spaceBetweenReadsInPoint=  spaceBetweenFeatures
        featureHeightInPoint = featureHeight
        featureHeight = featureHeightInPoint + spaceBetweenReadsInPoint
    } else{
        spaceBetweenReadsInPoint = featureHeight/8
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
        stop("plotNames is deprecated. The functionality is replaced by plot_feature_text")
        # nTrack = max(mybins)
        # if(plotBottomToTop) {
        #     thisy = mybins * (1/nTrack) + 1/nTrack/2
        # } else {
        #     thisy = 1 - mybins * (1/nTrack)  + 1/nTrack/2
        # }
        #
        # grid.text(plotNames,
        #     x = unit((start(x) + end(x))/2,"native"),
        #     y = unit(thisy, "npc"),gp=gpar(cex=0.4))
    }

}

# x = IRanges(start=c(10,40,20),end=c(20,50,80))
# mytext = c("a","b","c")
# pushViewport(viewport(width=0.5,height=0.5))
# pushViewport(dataViewport(xData=c(10,100), yscale=c(0,1), extension=0, clip="off"))
# grid.rect()
# plot_feature_text(x,mytext,fontsize=20,debug=TRUE,plotBottomToTop=FALSE)
plot_feature_text = function(x, text, fontsize=12, side=0, col="black",just="center", xjust=NULL, yjust=NULL, plotBottomToTop=TRUE, debug=FALSE){
    ## side: 0 center,1, left, 2,top 3,right 4 bottom
    ## xjust, yjust, shift in x or y if defined, currently ignored
    ## textSize: text height in points
    mybins = disjointBins(x)
    featureHeight = fontsize
    if(side!=0) stop("side other than center has not been implemented")
    myx = start(x)

    ## for - strand stack top to bottom, for + strand bottom to top
    myy = if(plotBottomToTop){### usually for "+" strand
         (mybins-1) * featureHeight}
         else { ## usually for "-" strand
        convertY(unit(1,"npc"),"points",valueOnly=TRUE) - mybins * featureHeight}
    if(debug){
        grid.rect(myx, unit(myy,"points"), width=width(x),
            height=unit(featureHeight,"points"), gp=gpar(col = "black" , fill = NA),
            default.units="native", just=c("left","bottom"))
    }
    ## for side other than 0 play with 1 strwidth and strheight
    ## use signif to make sure there are not too many digits after converting from npc to points
    if(side==0){
        xtext = (start(x) + end(x))/2
        ytext = myy + featureHeight/2
    } else if (side == 1) {
        ## substrct 1 strwidth
        xtext = (start(x) + end(x))/2
        ytext = myy + featureHeight/2
    } else if (side == 2) {
        xtext = (start(x) + end(x))/2
        #ytext = myy + featureHeight*1.5
    }

    grid.text(text, x =unit(xtext,"native"), y = unit(ytext,"points"),just=just,
        hjust = xjust, vjust = yjust, gp = gpar(col=col,fontsize=fontsize))
}

plot_feature_text_vpr  = function(x, text, vpr,coord, fontsize=12,side=0, col="black", just = "center", xjust=NULL, yjust=NULL, plotBottomToTop=TRUE, debug=FALSE) {
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
    plot_feature_text(x=x,text=text,fontsize=fontsize,side=side, col=col, just=just, xjust=xjust, yjust=yjust,
            plotBottomToTop=plotBottomToTop,debug=debug)
    popViewport()
}
