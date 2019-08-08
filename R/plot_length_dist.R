plot_length_dist = function(ss,main="", xlim=c(200,10000), midPoint=4000,stepSizeBeforeMidPoint=200, stepSizeAfterMidPoint=1000, nbins=500, doMedianLine=TRUE){
    ## https://stackoverflow.com/questions/24646594/how-to-improve-the-aspect-of-ggplot-histograms-with-log-scales-and-discrete-valu
    if(is.vector(ss)){
        dd = tibble(length=ss)
    } else if(is.data.frame(ss)) {
        dd = ss %>% mutate(length=sequence_length_template)
    } else {
        stop("Input either a vector or a data frame with the column sequence_length_template like the sequencing_summary file")
    }

    breaks = c(
        seq(xlim[1],midPoint,by=stepSizeBeforeMidPoint),
        seq(midPoint,xlim[2],by=stepSizeAfterMidPoint)[-1])

    # x.p25.log <- quantile(dd$length, 0.25)
    # x.median.log <- median(dd$length)
    # x.mean.log <- mean(dd$length)
    # x.p75.log <- quantile(dd$length, 0.75)

    p = ggplot(dd, aes(x=length)) +
      #geom_density(color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE) +
      geom_histogram(color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE,bins=nbins) +
      #stat_density(aes(y=..count..), color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE, n=1000) +
      scale_x_continuous(breaks=breaks,trans="log10",limits=xlim, oob=scales::squish) +
      #scale_y_continuous(breaks=c(0,125,250,375,500,625,750), expand=c(0,0)) +
      theme_bw() + ggtitle(main) +  labs(x="Read length in bp") +
      theme(axis.text.x = element_text(angle = 45,hjust = 1), plot.title = element_text(hjust = 0.5))
    if(doMedianLine){
        p = p + geom_vline(xintercept = median(dd$length), color = "firebrick", alpha = 1)
    }
      #geom_vline(xintercept = c(x.p25.log), color = "chocolate3", alpha = 1) +

      #geom_vline(xintercept = c(x.mean.log), color = "royalblue", alpha = 1) +
      #geom_vline(xintercept = c(x.p75.log), color = "chocolate3", alpha = 1) +
      # geom_text(x = x.p25.log, y = 200, label = paste("P25", round(x.p25.log, digits = 2), sep = '\n'), color = "chocolate3") +
      #geom_text(x = x.median.log, y = 200, label = paste("Median", round(x.median.log, digits = 2), sep = '\n'), color = "red") +
    #   geom_text(x = x.mean.log, y = 200, label = paste("Mean", round(x.mean.log, digits = 2), sep = '\n'), color = "royalblue") +
    #   geom_text(x = x.p75.log, y = 200, label = paste("P75", round(x.p75.log, digits = 2), sep = '\n'), color = "chocolate3")
    return(p)
}
