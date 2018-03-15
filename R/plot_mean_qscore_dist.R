plot_mean_qscore_dist = function(ss,main=""){
    ## https://stackoverflow.com/questions/24646594/how-to-improve-the-aspect-of-ggplot-histograms-with-log-scales-and-discrete-valu
    #x.p25.log <- log10(quantile(dd$mean_qscore_template, 0.25))
    x.median <- median(ss$mean_qscore_template)
    #x.mean.log <- log10(mean(dd$mean_qscore_template))
    #x.p75.log <- log10(quantile(dd$mean_qscore_template, 0.75))

    ggplot(ss, aes(x=mean_qscore_template)) +
      geom_histogram(color="firebrick", fill="firebrick", alpha=0.6,na.rm=TRUE,bins=100)+
      #stat_density(aes(y=..count..), color="black", fill="steelblue", alpha=0.5,na.rm=TRUE) +
      scale_x_continuous() +
      #scale_y_continuous(breaks=c(0,125,250,375,500,625,750), expand=c(0,0)) +
      theme_bw() + ggtitle(main) + labs(x = "Mean quality score per read", y = "Count") +
      theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept = c(x.median), color = "darkgreen", alpha = 1)
}
