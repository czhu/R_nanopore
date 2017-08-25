bivar_plot_heatscatter = function(x,y,log) {
    # commonTheme = list(labs(color="Density",fill="Density",
    #                      x="RNA-seq Expression",
    #                      y="Microarray Expression"),
    #                 theme_bw(),
    #                 theme(legend.position=c(0,1),
    #                       legend.justification=c(0,1)))
    # mycolramp = function(colpal){
    #     palette = switch(colpal,heat = c("grey", "darkblue",
    #             "red", "orange", "gold"), crazyred = c("#940000",
    #             "#A50000", "#FF5C5C", "#FFB9B9"), crazygreen = c("dark green",
    #             "#009700", "green", "#C0F5D0"), crazyblue = c("dark blue",
    #             "blue", "#7390EE", "light blue")
    #             )
    #     colorRampPalette(palette)
    # }
    ## heatscatter
    ## https://wresch.github.io/2012/11/06/ggplot2-smoothscatter.html
    ##
    dens = densCols(x, y, colramp = colorRampPalette(LSD::colorpalette("heat",10)),nbin=512)
    df = tibble(x=x,y=y, dens=dens)
    p= ggplot(df, aes(x, y)) + geom_point(aes(col=dens),size = 1) + scale_color_identity() + theme_bw()

    if(!missing(log)){
        if(grepl("x",log))  p = p + scale_x_log10()
        if(grepl("y",log))  p = p + scale_y_log10()
    }
    return(p)
}

bivar_plot_contour = function(x,y,addPoint=FALSE,log) {
    # https://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay
    df = tibble(x=x,y=y)
    p = ggplot(data=df,aes(x,y))  +
      stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',color=FALSE,contour=TRUE,bins=128) +
      scale_fill_continuous(high="#006400",low="#C0F5D0") + guides(fill="none",alpha="none") +
      theme_bw()

    if(addPoint) p = p +  geom_point(size=1,alpha=0.5)
    if(!missing(log)){
        if(grepl("x",log))  p = p + scale_x_log10()
        if(grepl("y",log))  p = p + scale_y_log10()
    }
    return(p)
}

### test purpose
# bivar_plot_contour(df$x,df$y)
# df <- data.frame(x = rnorm(1e3, 50, 10), y = rnorm(1e3, 50, 10))
# library(gridExtra)
# grid.arrange(bivar_plot_contour(df$x,df$y), bivar_plot_heatscatter(df$x,df$y) ,ncol=2)
#    ggExtra::ggMarginal(p, type = "histogram",bins=100)
