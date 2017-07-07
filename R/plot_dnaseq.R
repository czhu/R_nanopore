plot_dnaseq = function(myseq,ypos,add=FALSE,height=1,do.legend=TRUE,
    dnaCols = c(A="#00F728",T="#FF2A19",G="#000000",C="#1D49FB")) {
    myrle = rle(strsplit(myseq,"")[[1]])

    if(missing(ypos)){ypos=0.5}

    if(!add){
        plot(c(0, nchar(myseq)+1), c(ypos-1, ypos+height+1), type= "n", xlab = "", ylab = "")
    }

    xleft = cumsum(c(1,myrle$length))[-(length(myrle$length)+1)]
    xright = c(cumsum(myrle$length)[-1],nchar(myseq)+1)

    rect(xleft, ypos, xright, ypos+height, col=dnaCols[myrle$values],border = NA)

    if(do.legend){
        legend("topleft",names(dnaCols),fill=dnaCols,ncol=2)
    }
    ypos+height
}
