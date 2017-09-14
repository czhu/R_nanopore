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

gene_plot = function(plotDat, plotTxLabel=TRUE,doCDS = TRUE,debug=FALSE, drawSpaceBetweenReads = TRUE,
    drawPanelRect = TRUE, drawReadCol = TRUE, drawReadBorder = FALSE,
    doLine=TRUE, lineType= "dotted" ){
    ## default settings
    ## shall we draw CDS
    ## line between exons
    readColor = "steelblue"
    readAltColor = "gray40"
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
    coord = c(min(start(plotDat$exons)), max(end(plotDat$exons)))
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
        pushViewport(dataViewport(xData=coord, yscale=c(0,0.3), extension=0, clip="off",
                                  layout.pos.col=1, layout.pos.row=which(names(VP)=="coord")))
        grid.lines(coord, c(0,0), default.units = "native")
        tck = alongChromTicks(coord)
        grid.text(label=formatC(tck, format="d"), x = tck, y = 0.1,
                  just = c("centre", "bottom"), gp = gpar(cex=.4), default.units = "native")
        grid.segments(x0 = tck, x1 = tck, y0 = 0, y1 = 0.08,  default.units = "native")
        popViewport()
        ############

        ############ draw GeneModel  ## FIXME write a function?
        myfeature = plotDat$exons
        pushViewport(dataViewport(xData=coord, yscale=c(0,1), extension=0,
            clip="off", layout.pos.col=1, layout.pos.row=which(names(VP)=="GeneModel")))
        thisMaxHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
        thisHeight = min(10, thisMaxHeight/2)
        y0= (thisMaxHeight - thisHeight)/2 ## start in this middle

        grid.rect(x=end(myfeature),y=unit(y0, "points"),width=width(myfeature),
            height=unit(thisHeight, "points"), gp=gpar(col = NA ,fill = "gray40", alpha=0.8),
            default.units="native", just=c("right","bottom"))
        ### no line for disjoint exons
        #myx = c(end(myfeature)[-length(myfeature)],start(myfeature)[-1])
        #grid.polyline(x=myx, y=unit(rep((y0+thisHeight)/2,length(myx)), "points"),id = rep(1:(length(myfeature)-1),2),
        #    gp=gpar(col="grey", lwd=2,alpha=0.5), default.units = "native")
        popViewport()
        #############

        #### data
        for(txIdx in 1:plotDat$param$nTrack){
            ########## draw tx annotation
            pushViewport(dataViewport(xData=coord, yscale=c(0,1), extension=0, clip="off",
                layout.pos.col=1,layout.pos.row=which(names(VP)==paste0("annot_",txIdx))))
            myfeature= plotDat$tx[[txIdx]]  ## !!!sorted GRanges!!!!

            thisMaxHeight = convertY(unit(1,"npc"),"points",valueOnly=TRUE)
            thisHeight = min(10, thisMaxHeight/2)
            y0= (thisMaxHeight - thisHeight)/2

            grid.rect(end(myfeature),unit(y0,"points"),width=width(myfeature),height=unit(thisHeight,"points"),
                gp=gpar(col = NA ,fill = "firebrick",alpha=0.5), default.units="native", just=c("right","bottom"))
            ### grid line for annotattion
            myx = c(end(myfeature)[-length(myfeature)],start(myfeature)[-1])
            if(length(myx)>0){
                grid.polyline(x=myx, y=unit(rep(y0+thisHeight/2, length(myx)),"points"),
                    id = rep(1:(length(myfeature)-1),2),
                    gp=gpar(col="firebrick", lwd=1,alpha=0.5,lty=lineType),
                    default.units = "native")
            }
            id(plotTxLabel) grid.text(plotDat$param$trackNames[txIdx],y = unit(0.7, "npc"), gp=gpar(cex=0.4))

            if(doCDS) {
                myfeature=plotDat$cds[[txIdx]]
                if(length(myfeature)>0){
                    extraHeight = min(5,thisMaxHeight/8)
                    grid.rect(end(myfeature),unit(y0 - extraHeight,"points"),
                    width=width(myfeature),
                    height=unit(thisHeight + 2*extraHeight,"points"),
                        gp=gpar(col = NA ,fill = "firebrick",alpha=0.5),
                            default.units="native", just=c("right","bottom"))
                }
            }
            popViewport()

            ###### draw data for isoform x
            if(plotDat$param$countMat[txIdx,vpCol]>0){
                pushViewport(dataViewport(xData=coord, yscale=c(0,1), extension=0, clip="on",
                    layout.pos.col=1,layout.pos.row=which(names(VP)==paste0("data_",txIdx))))
                myfeature= sort(grglist(plotDat$data[[vpCol]][[txIdx]], drop.D.ranges=FALSE))
                #eachReadSpace = min(2,convertY(unit(1,"npc"),"points",valueOnly=TRUE)/length(myfeature))
                readHeightInPoint = readHeight - spaceBetweenReadsInPoint
                myx = unlist(end(myfeature))
                ## for - strand stack top to bottom, for + strand bottom to top
                if(plotDat$param$strand == "+"){
                    yPerRead = seq(0, by=readHeight, length.out=length(myfeature))
                    } else{
                    yPerRead = seq(convertHeight(unit(1,"npc"),"points",valueOnly=TRUE),
                        by=-readHeight, length.out=length(myfeature))
                }
                nFeatureEach = sapply(end(myfeature),length)

                myy = rep(yPerRead, nFeatureEach)

                if(drawReadCol){
                    # mycols = rep(
                    #     ifelse(mcols(plotDat$data[[vpCol]][[txIdx]])$is_fulllength, readColor,readAltColor),
                    #     nFeatureEach) ### FIXME this should be a option
                    mycols = rep(num_to_color(mcols(plotDat$data[[vpCol]][[txIdx]])$qscore,from=c(5,12)), nFeatureEach)
                } else {mycols="steelblue"}

                grid.rect(myx,unit(myy,"points"), width=unlist(width(myfeature)),
                    height=unit(readHeightInPoint,"points"), gp=gpar(col = NA , fill = mycols, alpha=0.8),
                    default.units="native", just=c("right","bottom"))

                myxStart = lapply(end(myfeature),function(x) x[-length(x)])
                myxEnd = lapply(start(myfeature),function(x) x[-1])
                myy = c(rep(yPerRead, sapply(myxStart,length)), rep(yPerRead, sapply(myxEnd,length)))

                if(doLine & length(unlist(c(myxStart,myxEnd)))>0){
                    penaltyFactorReadNumber = (1/log10(plotDat$param$normCountMat[txIdx,vpCol]))^2
                    grid.polyline(
                        x=unlist(c(myxStart,myxEnd)), y=unit(myy+readHeightInPoint/2,"points"),
                        id = rep(1:length(unlist(myxStart)),2),
                        gp=gpar(col=mycols, lwd=unit(min(1,readHeight/3)*penaltyFactorReadNumber,"points"),
                        ## FIXME scale alpha depending on the number of reads
                            alpha=0.2*penaltyFactorReadNumber,lty=lineType), ##lex=1/penaltyFactorReadNumber),
                        default.units = "native")
                }
                if(debug)
                    grid.text(
                        paste("Raw", plotDat$param$countMat[txIdx,vpCol],
                            "Norm", simplify_num_output(plotDat$param$normCountMat[txIdx,vpCol]),
                            "SF", simplify_num_output(plotDat$param$sizeFactor[vpCol])),
                        y = unit(ifelse(plotDat$param$strand == "+", 0.9, 0.1), "npc"), gp=gpar(cex=0.4))
                #mypos = pretty(c(1,length(myfeature)))
                #grid.yaxis(at=convertY(unit(mypos * eachReadSpace,"points"),"npc",valueOnly=TRUE),label=mypos,gp=gpar(col="red"))
                popViewport()
            }

        }
        popViewport()
        popViewport()
    }

    popViewport(2)
}
