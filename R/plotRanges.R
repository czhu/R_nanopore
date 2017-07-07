plotRanges = function(x,ypos,add=FALSE, col= "black", xlim,ylim, sep = 0.2 ,height=1, ...){
    if(missing(ypos)) {ypos=0}

    bins <- disjointBins(IRanges(start(x), end(x)))
    ybottom <- ypos + bins * (sep + height) - height

    if(!add) {
        if(missing(xlim)) xlim = c(min(start(x))-1, max(end(x))+1)
        if(missing(ylim)) ylim = c(min(ybottom) - 1, max(ybottom+height)+1)
        plot.new()
        plot.window(xlim, ylim)
    }

    rect(start(x), ybottom, end(x), ybottom + height, col = col,border=NA, ...)

    if(!add){
        axis(1)
    }
    max(ybottom + height)+1
}
