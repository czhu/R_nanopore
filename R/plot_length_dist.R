plot_length_dist = function(ss,main=""){
    ## https://stackoverflow.com/questions/24646594/how-to-improve-the-aspect-of-ggplot-histograms-with-log-scales-and-discrete-valu
    lowEnd = 200
    highEnd = 10000
    midPoint = 4000 # below this go with stepSizeShort above this with stepSizeLong for breaks
    stepSizeShort = 200
    stepSizeLong = 1000
    breaks = c(seq(lowEnd,midPoint,by=stepSizeShort),seq(midPoint,highEnd,by=stepSizeLong)[-1])

    dd = ss %>% mutate(length = replace(sequence_length_template,sequence_length_template<lowEnd,lowEnd)) %>%
        mutate(length = replace(length,length>highEnd,highEnd))

    x.p25.log <- quantile(dd$length, 0.25)
    x.median.log <- median(dd$length)
    x.mean.log <- mean(dd$length)
    x.p75.log <- quantile(dd$length, 0.75)

    ggplot(dd, aes(x=length)) +
      #geom_density(color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE) +
      geom_histogram(color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE,bins=100) +
      #stat_density(aes(y=..count..), color="steelblue", fill="steelblue", alpha=0.8,na.rm=TRUE, n=1000) +
      scale_x_continuous(breaks=breaks,trans="log10") +
      #scale_y_continuous(breaks=c(0,125,250,375,500,625,750), expand=c(0,0)) +
      theme_bw() + ggtitle(main) +  labs(x = "Read length in BP") +
      theme(axis.text.x = element_text(angle = 45,hjust = 1), plot.title = element_text(hjust = 0.5))+
      #geom_vline(xintercept = c(x.p25.log), color = "chocolate3", alpha = 1) +
      geom_vline(xintercept = c(x.median.log), color = "firebrick", alpha = 1) + xlim(c(lowEnd,highEnd))
      #geom_vline(xintercept = c(x.mean.log), color = "royalblue", alpha = 1) +
      #geom_vline(xintercept = c(x.p75.log), color = "chocolate3", alpha = 1) +
      # geom_text(x = x.p25.log, y = 200, label = paste("P25", round(x.p25.log, digits = 2), sep = '\n'), color = "chocolate3") +
      #geom_text(x = x.median.log, y = 200, label = paste("Median", round(x.median.log, digits = 2), sep = '\n'), color = "red") +
    #   geom_text(x = x.mean.log, y = 200, label = paste("Mean", round(x.mean.log, digits = 2), sep = '\n'), color = "royalblue") +
    #   geom_text(x = x.p75.log, y = 200, label = paste("P75", round(x.p75.log, digits = 2), sep = '\n'), color = "chocolate3")

}
