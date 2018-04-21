#mylog = function(x) log2(x+1)
mylog = function(x) log(x+1,1.5)

trsf_ct_nread = function(x) {
    ### transform count to
    mylog(x)
}

draw_rect = function(vpr) {
    pushViewport(viewport(layout.pos.row=vpr))
    grid.rect(gp=gpar(col="gray",alpha=0.3,lty="dashed"))
    popViewport()
}


### input data
validate_plotdat = function(x){
    stopifnot(!is.null(x$consensus$count))
    return( all(names(x$consensus) == names(x$reads)) )
}

chrom_plot = function(plotDat,coord, plotCountNum=TRUE,featureHeightPerRead = 3, spaceBetweenCluster = 5, debug = TRUE,
    doConsensus=TRUE){
    # x is plot data
    # plotDat = list(
    #     consensus = thisCluster, this should contian count
    #     geneModel = knownGene,
    #     reads = thisReads
    # )
    stopifnot(validate_plotdat(plotDat))
    ## a few plot paramter
    myannot = plotDat$geneModel
    doHighlight = !is.null(plotDat$highlight)
    genomeAxisHeight = 10

    extendLeft = 250 ## in bp
    extendRight = 250

    ct_to_drawPoint  = function(x){
        trsf_ct_nread(x)*featureHeightPerRead*2
    }

    spacePerAnnotTrack = rep(6,2)
    names(spacePerAnnotTrack) = c("Annot_plus","Annot_minus")

    ## height of the page
    wholeSpaceHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
    spacePerDataTrack = c(
        (wholeSpaceHeight - genomeAxisHeight)/2 - spacePerAnnotTrack["Annot_plus"],
        (wholeSpaceHeight - genomeAxisHeight)/2 - spacePerAnnotTrack["Annot_minus"])
    names(spacePerDataTrack) = c("Data_plus","Data_minus")

    ## ie, ei, coverage should really be optionnal
    # VP = c(
    #     spacePerDataTrack["Data_plus"], spacePerAnnotTrack["Annot_plus"],"coverage"=10,
    #     "axis"=genomeAxisHeight,"ie"=5,"ei"=5,
    #     spacePerAnnotTrack["Annot_minus"], spacePerDataTrack["Data_minus"])
    doHighlightGene = !is.null(plotDat$gene)
    if(doHighlightGene) {
        mygenes = plotDat$gene
        VP = c(
            spacePerDataTrack["Data_plus"], spacePerAnnotTrack["Annot_plus"], Gene_plus = 5,
            "axis"=genomeAxisHeight, Gene_minus = 5, spacePerAnnotTrack["Annot_minus"], spacePerDataTrack["Data_minus"])
    } else {
        VP = c(
            spacePerDataTrack["Data_plus"], spacePerAnnotTrack["Annot_plus"],
            "axis"=genomeAxisHeight, spacePerAnnotTrack["Annot_minus"], spacePerDataTrack["Data_minus"])
    }

    ######### plot title
    mytitle=plotDat$name
    doTitle = !is.null(mytitle)
    if(doTitle){
        titleVerticalHeight = unit(2,"line")
        pushViewport(viewport(x=0,y=1,width=1,height=titleVerticalHeight,just=c("left","top")))
        grid.text(mytitle)
        popViewport()
        pushViewport(viewport(x=0,y=unit(1,"npc")-titleVerticalHeight,
            width=1,height=unit(1,"npc")-titleVerticalHeight, just=c("left","top")))
    }
    ##########

    grl=grid.layout(length(VP), 1, heights=unit(VP,"points"))
    pushViewport(viewport(layout=grl))
    #grid.show.layout(glo,newpage=FALSE)
    if(debug) sapply(1:length(VP),draw_rect)

    plot_coord(coord,which(names(VP)=="axis"))

    # regionName = paste0(as.character(seqnames(myregion)),":", coord[1],"-", coord[2])
    #
    # make_bar_plot( myregion = regionName, coord, tgzFile = shortreadFile,
    #     vpr = which(names(VP)=="coverage"),col="darkgreen",trsf=identity,lwd=0.1)
    #
    # make_bar_plot(
    #     myregion = regionName, coord, tgzFile = splicesiteFileEI,
    #     vpr = which(names(VP)=="ei"),col="purple4",trsf=log2,lwd=0.25)
    #
    # make_bar_plot(
    #     myregion = regionName, coord, tgzFile = splicesiteFileIE,
    #     vpr = which(names(VP)=="ie"),col="saddlebrown",trsf=log2,lwd=0.25)

    ## annotation track
    if(doHighlightGene){
        if(any(strand(mygenes)=="+")) {
            plot_feature_vpr(mygenes[strand(mygenes)=="+"],vpr=which(names(VP)=="Gene_plus"),coord=coord,
            featureHeight=0.5, plotBottomToTop=TRUE, featureCols="black",doLine=FALSE,center=TRUE,spaceBetweenFeatures=1)
        }
        if(any(strand(mygenes)=="-")) {
            plot_feature_vpr(mygenes[strand(mygenes)=="-"],vpr=which(names(VP)=="Gene_minus"),coord=coord,
            featureHeight=0.5, plotBottomToTop=FALSE,featureCols="black",doLine=FALSE, center=TRUE,spaceBetweenFeatures=1)
        }
    }

    if(any(strand(myannot)=="+")) {
        plot_feature_vpr(plotDat$geneModel[strand(plotDat$geneModel)=="+"],vpr=which(names(VP)=="Annot_plus"),coord=coord,
        featureHeight=3, plotBottomToTop=TRUE, featureCols="firebrick",doLine=FALSE,center=TRUE)
    }
    if(any(strand(myannot)=="-")) {
        plot_feature_vpr(plotDat$geneModel[strand(plotDat$geneModel)=="-"],vpr=which(names(VP)=="Annot_minus"),coord=coord,
        featureHeight=3, plotBottomToTop=FALSE,featureCols="firebrick",doLine=FALSE, center=TRUE)
    }

    ## draw data track
    for(thisStrd in c("+","-")){
        isThisStrand = strand(plotDat$consensus) == thisStrd
        if(any(isThisStrand)){
            pushViewport(viewport(xscale=coord,clip="off",
                layout.pos.col=1,layout.pos.row=which(names(VP)==ifelse(thisStrd=="+","Data_plus","Data_minus")))
            )
            #if(debug) grid.rect()
            # if(as.vector(strand(higlightRegion)==thisStrd)){
            #     grid.lines(x=c(start(higlightRegion),start(higlightRegion)),y=unit(c(0,1),"npc"),
            #         default.units="native",gp=gpar(col="firebrick",alpha=0.4,lwd=0.8))
            #     grid.lines(x=c(end(higlightRegion),end(higlightRegion)),y=unit(c(0,1),"npc"),
            #         default.units="native",gp=gpar(col="firebrick",alpha=0.4,lwd=0.8))
            # }

            thisSpaceHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
            #nMaxClusterPlotted = floor(thisSpaceHeight/minHeightPerClusterFeature)

            subconsensusFeature = plotDat$consensus[isThisStrand]

            ## find transcript unit
            expressedUnit = GenomicRanges::reduce(subconsensusFeature)
            ovlps = as.list(findOverlaps(expressedUnit, subconsensusFeature))

            ## plotting by transcript unit
            for(expressedUnitIndex in 1:length(expressedUnit)){
                subconsensusFeatureWithinEU = subconsensusFeature[ovlps[[expressedUnitIndex]]]

                mybins = disjointBins(subconsensusFeatureWithinEU)

                countPerTrack = tapply(subconsensusFeatureWithinEU$count,mybins,max)
                countPerTrack = countPerTrack[order(countPerTrack,decreasing=TRUE)]
                trackHeightShift = cumsum(ct_to_drawPoint(countPerTrack)+spaceBetweenCluster)

                for(trackIndex in 1:length(countPerTrack)) {
                    thisBin = as.integer(names(countPerTrack)[trackIndex])
                    for(clusterIndex in which(mybins == thisBin)) {
                        thisCluster = subconsensusFeatureWithinEU[clusterIndex]
                        thisReads = plotDat$reads[[names(thisCluster)]]
                        thisx = unit(start(thisCluster),"native")
                        thisy = unit(
                            if(thisStrd=="-"){
                                thisSpaceHeight- trackHeightShift[trackIndex]
                            } else {
                                if( trackIndex==1 ) {spaceBetweenCluster} else {
                                    trackHeightShift[trackIndex-1]+spaceBetweenCluster
                                }
                            },"points")
                        thisHeight = unit(ct_to_drawPoint(thisCluster$count),"points")
                        pushViewport(
                            viewport(
                            x = thisx, y = thisy,
                            width = unit(width(thisCluster),"native"),
                            clip="off",just=c("left","bottom"),
                            xscale=c(start(thisCluster),end(thisCluster)),
                            height=thisHeight
                            )
                        )
                        if(doConsensus){
                            featureHeightConsensus = unit(5,"points")
                            if(thisStrd=="+"){
                                newy1 = 0
                                newy2 = featureHeightConsensus
                            } else{
                                newy1 = unit(round(convertY(unit(1,"npc"),"points",valueOnly=TRUE) - as.integer(featureHeightConsensus)),"points")
                                newy2 = 0
                            }
                            # message(thisStrd,"\t",convertY(thisHeight-featureHeightConsensus,"points",valueOnly=TRUE),
                            #     "\t",convertY(unit(1,"npc"),"points",valueOnly=TRUE),"\t", as.integer(newy1))
                            pushViewport(
                                viewport(
                                x = thisx,
                                y = newy1,
                                width = unit(width(thisCluster),"native"),
                                clip="off",just=c("left","bottom"),
                                xscale=c(start(thisCluster),end(thisCluster)),
                                height= featureHeightConsensus
                                ))
                            plot_feature(thisCluster,
                                featureCols=ifelse(is.null(thisCluster$itemRgb),"black", thisCluster$itemRgb),
                                    featureHeight=as.integer(featureHeightConsensus),
                                    doLine=TRUE,lineAlpha=1,lineType= "dotted",
                                    plotBottomToTop = ifelse(thisStrd=="+",TRUE,FALSE),
                                    center=TRUE)
                            if(doHighlight){
                                ## loop because thisCluster could have multiple hits
                                for(thisName in names(thisCluster)){
                                    wh = which(names(plotDat$consensus) == thisName)
                                    if(!is.null(plotDat$highlight[[wh]])){
                                        grid.text(plotDat$highlight[[wh]]$shape,
                                            x= unit(convertX(unit(start(plotDat$consensus[wh]) - extendLeft,"native"),"npc",
                                                valueOnly=TRUE)-convertX(unit(1,"strwidth","s"),"npc",valueOnly=TRUE),"npc"),
                                            0.5,just=c("left","center"),gp=gpar(fontsize=4))
                                        if(!is.null(plotDat$highlight[[wh]]$highlight)){
                                            thisStart = start(plotDat$highlight[[wh]]$highlight)
                                            thisEnd = end(plotDat$highlight[[wh]]$highlight)
                                            plot_feature(plotDat$highlight[[wh]]$highlight,
                                                featureCols="gray4",
                                                    featureHeight=as.integer(featureHeightConsensus)-2,
                                                    doLine=FALSE,featureAlpha=0.5,
                                                    plotBottomToTop = ifelse(thisStrd=="+",TRUE,FALSE),
                                                    center=TRUE)
                                            grid.text(plotDat$highlight[[wh]]$highlight$shape,
                                                x = convertX(unit((thisStart+thisEnd)/2,"native"),"npc"),
                                                y=0.5,just="center",gp=gpar(fontsize=4))
                                        }
                                    }
                                }
                            }
                            popViewport()
                            pushViewport(
                                viewport(
                                x= thisx,
                                y = newy2,
                                width = unit(width(thisCluster),"native"),
                                clip="off",just=c("left","bottom"),
                                xscale=c(start(thisCluster),end(thisCluster)),
                                height = thisHeight - featureHeightConsensus
                                )
                            )
                        }
                        thisMaxHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
                        nReadsToPlot = floor(thisMaxHeight/featureHeightPerRead)
                        nReadsToPlot = ifelse(nReadsToPlot> length(thisReads), length(thisReads), nReadsToPlot)

                        x = thisReads[order(extract_qscore(thisReads$name),decreasing=TRUE)][1:nReadsToPlot]
                        mycols = rep(num_to_color(mylog(thisCluster$count),from=mylog(c(1,256)),to=c(0.4,1)))
                        plot_feature(x,
                            featureCols=mycols,
                            featureHeight=featureHeightPerRead,
                            doLine=TRUE,lineAlpha=1,lineType= "dotted",
                            plotBottomToTop = ifelse(thisStrd=="+",TRUE,FALSE),
                            center=TRUE)
                        if(doConsensus) popViewport()

                        if(plotCountNum){
                            s = as.character(thisCluster$count)
                            grid.text(s,
                                x= unit(convertX(unit(min(start(x)) - extendLeft,"native"),"npc",
                                    valueOnly=TRUE)-convertX(unit(1,"strwidth","s"),"npc",valueOnly=TRUE),"npc"),
                                0.5,
                            just=c("left","center"),gp=gpar(fontsize=4
                                #convertY(unit(0.5,"npc"),"points",valueOnly=TRUE)
                            ))
                        }

                        if(debug) {
                            #grid.rect(gp=gpar(col="gray",alpha=0.4,lineType= "dashed"))
                            x0 = min(start(x)) - extendLeft
                            x1 = max(end(x)) + extendRight
                            y0 = convertHeight(unit(0,"npc"),"native",valueOnly=TRUE)
                            y1 = convertHeight(unit(1,"npc"),"native",valueOnly=TRUE)

                            grid.lines(x=c(x0,x0),y=c(y0,y1),default.units="native",gp=gpar(col="black",lineType= "dashed",alpha=0.6,lwd=0.4))
                            grid.lines(x=c(x0,x1),y=c(y0,y0),default.units="native",gp=gpar(col="gray",lineType= "dashed",alpha=0.2,lwd=0.2))
                            grid.lines(x=c(x1,x1),y=c(y0,y1),default.units="native",gp=gpar(col="gray",lineType= "dashed",alpha=0.2,lwd=0.2))
                            grid.lines(x=c(x0,x1),y=c(y1,y1),default.units="native",gp=gpar(col="gray",lineType= "dashed",alpha=0.2,lwd=0.2))
                        }

                        popViewport()
                    }
                }
            }
            popViewport()
        }
    }
    popViewport()
    if(doTitle) popViewport()
}
