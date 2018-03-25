genome_barplot = function(x, coord, vpr, trsf=identity,col="darkgreen",alpha=1,lwd=1){
    ### x is data.frame or tibble in the format of bedGraph
    ### i.e chr, start, end, value,
    ### trsf function to transform value i.e log10 or other function

    if(missing(vpr)) {
        vpr = new_vp()
    }

    xmat = as.matrix(cbind(x[,2],x[,2:3],x[,3]))
    ymat = as.matrix(cbind(0, trsf(x[,4]),trsf(x[,4]),0))

    if(missing(coord)) {
        coord = c(min(xmat), max(xmat))
    }

    pushViewport(dataViewport(xData=coord,yData=as.vector(ymat),
            extension=0, clip="off",layout.pos.col=1, layout.pos.row=vpr))

    grid.polygon( as.vector(t(xmat)),as.vector(t(ymat)), default.units="native",
            gp=gpar(col=col,fill=col,alpha=alpha,lwd=lwd) )

    popViewport()
}

make_bar_plot = function(myregion,coord,
    tgzFile, vpr,trsf=identity,col="darkgreen",alpha=1,lwd=1){
    ## convience function to work with bedGraph file for coverag
    if(missing(tgzFile)) {
        stop("A bedGraph for coverage must be provided")
    }
    mytmpfile = tempfile()
    system(paste("tabix", tgzFile, myregion,">", mytmpfile))
    tsv = suppressMessages(read_tsv(mytmpfile,col_names = FALSE,progress=FALSE))
    genome_barplot(tsv, coord = coord, vpr=vpr, trsf = trsf,col=alpha, alpha=alpha,lwd=lwd)
}
