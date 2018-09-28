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

### red
# default_config = function(...){
#     ### height in points
#     list(
#         gene = list(color="black"),
#         geneModel=list(color="#B22222", height=3),
#         geneModel_reduced = list( color = "#EF2D2D", height=3 ),
#         readConsensus=list(color="#4daf4a", height=5),
#         read=list(color="steelblue")
#     )
# }

## purple
default_config = function(...){
    ### height in points
    list(
        gene = list(color="black"),
        geneModel=list(color="#984ea3", height=3),
        geneModel_reduced = list( color = "#B771C2", height=3 ),
        readConsensus=list(color="#4daf4a", height=5),
        read=list(color="steelblue")
    )
}

## tan
# default_config = function(...){
#     ### height in points
#     list(
#         gene = list(color="black"),
#         geneModel=list(color="#CD853F", height=3),
#         geneModel_reduced = list( color = "#F1A668", height=3 ),
#         readConsensus=list(color="#4daf4a", height=5),
#         read=list(color="steelblue")
#     )
# }

### input data
validate_plotdat = function(x){
    stopifnot(!is.null(x$consensus$count))
    return( all(names(x$consensus) == names(x$reads)) )
}

HIGHLIGHT_FONTSIZE = 4
CONSENSUS_NAME_FONTSIZE = 3
DO_CONSENSUS_NAME = TRUE

chrom_plot = function(plotDat,coord, plotCountNum=TRUE,featureHeightPerRead = 3,
    spaceBetweenCluster = 5, debug = TRUE,doConsensus=TRUE, config, singleStrand=FALSE){
    # x is plot data
    # plotDat = list(
    #     consensus = thisCluster, this should contian count
    #     geneModel = knownGene,
    #     reads = thisReads
    # )
    ## CHANGED 2018-09-28
    ## add singleStrand for single strand mode for gene_plot_log
    ## all data and annotation should have minus strand FIXME this is currently not checked

    stopifnot(validate_plotdat(plotDat))
    ## a few plot paramter
    myannot = plotDat$geneModel
    doHighlight = !is.null(plotDat$highlight)
    doHighlightGene = !is.null(plotDat$gene)
    doRedcuedGene = !is.null(plotDat$geneModel_reduced)
    if(missing(config)){
        config = default_config()
    }
    genomeAxisHeight = 10

    extendLeft = 200 ## in bp
    extendRight = 200

    ct_to_drawPoint  = function(x){
        trsf_ct_nread(x)*featureHeightPerRead*2
    }

    spacePerAnnotTrack = rep(2,2)
    names(spacePerAnnotTrack) = c("Annot_plus","Annot_minus")

    ## height of the page
    wholeSpaceHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
    # spacePerDataTrack = c(
    #     (wholeSpaceHeight - genomeAxisHeight)/2 - spacePerAnnotTrack["Annot_plus"],
    #     (wholeSpaceHeight - genomeAxisHeight)/2 - spacePerAnnotTrack["Annot_minus"])
    # names(spacePerDataTrack) = c("Data_plus","Data_minus")

    ## ie, ei, coverage should really be optionnal
    # VP = c(
    #     spacePerDataTrack["Data_plus"], spacePerAnnotTrack["Annot_plus"],"coverage"=10,
    #     "axis"=genomeAxisHeight,"ie"=5,"ei"=5,
    #     spacePerAnnotTrack["Annot_minus"], spacePerDataTrack["Data_minus"])
    ## extra 2 px spacing between data and annotation
    spaceBetweenAnnotationAndData = 2

    if(singleStrand){
        VP = c(
            "axis"=genomeAxisHeight,
            Gene_minus = ifelse(doHighlightGene, 5, 0),
            Gene_name_minus = ifelse(doHighlightGene, 5, 0),
            spacePerAnnotTrack["Annot_minus"], spaceBetweenAnnotationAndData)
        VP = c(VP, "Data_minus" = wholeSpaceHeight - sum(VP) )
    } else {
        VP = c(
            spaceBetweenAnnotationAndData,  spacePerAnnotTrack["Annot_plus"],
            Gene_name_plus=ifelse(doHighlightGene, 5, 0), Gene_plus = ifelse(doHighlightGene, 5, 0),
            "axis"=genomeAxisHeight,
            Gene_minus = ifelse(doHighlightGene, 5, 0), Gene_name_minus =ifelse(doHighlightGene, 5, 0),
            spacePerAnnotTrack["Annot_minus"], spaceBetweenAnnotationAndData)
        VP = c("Data_plus" = (wholeSpaceHeight - sum(VP))/2, VP, "Data_minus" = (wholeSpaceHeight - sum(VP))/2 )
    }

    if(doRedcuedGene){
        if( !singleStrand ){
            indexToInsert=which(names(VP) =="Annot_plus")-1
            VP = append(VP, spacePerAnnotTrack["Annot_plus"],indexToInsert)
            names(VP)[indexToInsert+1] = "Annot_reduced_plus"
        }

        indexToInsert=which(names(VP) =="Annot_minus")
        VP = append(VP, spacePerAnnotTrack["Annot_minus"],indexToInsert)
        names(VP)[indexToInsert+1] = "Annot_reduced_minus"
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

    grl = grid.layout(length(VP), 1, heights=unit(VP,"points"))
    pushViewport(viewport(layout=grl, clip="off"))
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

    ## draw data track
    if(singleStrand){
        strds = "minus"
        names(strds) =  "-"
    } else {
        strds = c("plus", "minus")
        names(strds) = c("+", "-")
    }

    for(thisStrd in names(strds) ){

        if(doHighlightGene & length(plotDat$gene)>0){
            if(any(strand(plotDat$gene)== thisStrd)) {
                isThiStrd = as.logical(strand(plotDat$gene)== thisStrd)
                plot_feature_vpr(plotDat$gene[isThiStrd],vpr=which(names(VP)== paste0("Gene_", strds[thisStrd]) ),
                coord=coord, featureHeight=0.5, plotBottomToTop=TRUE, featureCols= config$gene$color,
                    doLine=FALSE,center=TRUE,spaceBetweenFeatures=1)
                plot_feature_text_vpr(plotDat$gene[isThiStrd],
                    plotDat$gene$name[isThiStrd], vpr=which(names(VP)== paste0("Gene_name_", strds[thisStrd])),
                    coord, fontsize=HIGHLIGHT_FONTSIZE, side=0, col="black",xjust=unit(0,"npc"), yjust=y(0,"npc"),
                    plotBottomToTop=TRUE,debug=FALSE)
            }
        }

        if(any(strand(myannot)== thisStrd)) {
            thisFeature = subset(myannot, strand== thisStrd)
            plot_feature_vpr( thisFeature, vpr=which(names(VP)== paste0("Annot_", strds[thisStrd])),coord=coord,
                featureHeight=config$geneModel$height, plotBottomToTop=TRUE,
                featureCols= if(is.null(thisFeature$itemRgb)) {config$geneModel$color} else {thisFeature$itemRgb},
                doLine=FALSE,center=TRUE)
        }

        if(doRedcuedGene){
            myannot2 = plotDat$geneModel_reduced
            if(any(strand(myannot2)== thisStrd)) {
                plot_feature_vpr(subset(myannot2, strand== thisStrd),vpr=which(names(VP)== paste0("Annot_reduced_", strds[thisStrd])),
                coord=coord, featureHeight=config$geneModel_reduced$height, plotBottomToTop=TRUE,
                featureCols=config$geneModel_reduced$color,
                doLine=FALSE,center=TRUE)
            }
        }

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
                            featureHeightConsensus = unit(config$readConsensus$height, "points")
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
                                featureCols = if(is.null(thisCluster$itemRgb)) {"black"} else {thisCluster$itemRgb},
                                    featureHeight=as.integer(featureHeightConsensus),
                                    doLine=TRUE,lineAlpha=1,lineType= "dotted",
                                    plotBottomToTop = ifelse(thisStrd=="+",TRUE,FALSE),
                                    center=TRUE)
                            if(DO_CONSENSUS_NAME){
                                plot_feature_text(
                                    thisCluster,
                                    thisCluster$name, fontsize=CONSENSUS_NAME_FONTSIZE, side=0, col="black",
                                    xjust=unit(0,"npc"), yjust=y(0,"npc"),
                                    plotBottomToTop = (thisStrd =="+"), debug=FALSE)
                            }
                            if(doHighlight){
                                ## loop because thisCluster could have multiple hits
                                for(thisName in names(thisCluster)){
                                    wh = which(names(plotDat$consensus) == thisName)
                                    if(!is.null(plotDat$highlight[[wh]])){
                                        grid.text(plotDat$highlight[[wh]]$shape,
                                            x= unit(convertX(unit(start(plotDat$consensus[wh]) - extendLeft,"native"),"npc",
                                                valueOnly=TRUE)-convertX(unit(1,"strwidth","s"),"npc",valueOnly=TRUE),"npc"),
                                            0.5,just=c("left","center"),gp=gpar(fontsize=HIGHLIGHT_FONTSIZE))
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
                                                y=0.9,just="center",gp=gpar(fontsize=HIGHLIGHT_FONTSIZE))
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
                                    valueOnly=TRUE)-convertX(unit(1,"strwidth",s),"npc",valueOnly=TRUE),"npc"),
                                0.5,
                            just=c("left","center"),gp=gpar(fontsize=HIGHLIGHT_FONTSIZE
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
