# Author: czhu
###############################################################################

plot_length_dist_percentile = function(ss,...){
    dd = ss %>% mutate(length = replace(sequence_length_template,sequence_length_template<lowEnd,lowEnd)) %>%
        mutate(length = replace(length,length>highEnd,highEnd))

    highPoint = 1e4
    stepSize = 200

    breakPoints = c(seq(0,highPoint,by=stepSize),max(ss$sequence_length_template))
    lbp = length(breakPoints)
    mylabels = c(paste0("<",breakPoints[1]),paste(breakPoints[-c(1,lbp-1,lbp)],
        breakPoints[-c(1:2,lbp)],sep="-"), paste0(">",breakPoints[lbp-1]))

    xat = barplot(
        table(
            cut(ss$sequence_length_template,
                breakPoints,
                labels=mylabels
            )
        )/nrow(ss),
        xaxt="n",
        xlab="Bins", main=paste0("Percentile of reads per bin\n","total ",nrow(ss)),...)

    axis(1,xat,mylabels,srt=45,cex=0.5,xpd=TRUE)
}
