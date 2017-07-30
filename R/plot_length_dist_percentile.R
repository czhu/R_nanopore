# Author: czhu
###############################################################################

plot_length_dist_percentile = function(ss,...){

    highPoint = 1e4
    stepSize = 200

    breakPoints = c(seq(0,highPoint,by=stepSize),max(ss$sequence_length_template))
    lbp = length(breakPoints)
    mylabels = c(paste0("<",breakPoints[2]),paste(breakPoints[-c(1,lbp-1,lbp)],
        breakPoints[-c(1:2,lbp)],sep="-"), paste0(">",breakPoints[lbp-1]))

    xat = barplot(
        table(
            cut(ss$sequence_length_template,
                breakPoints,
                labels=mylabels
            )
        )/nrow(ss),
        #xaxt="n",xpd=TRUE
        xlab="Bins", main=paste0("Percentile of reads per bin\n","total ",nrow(ss)),las=2,
        cex.names=0.65,...)

    #axis(1,xat,mylabels,srt=45,cex=0.5,xpd=TRUE)
}
