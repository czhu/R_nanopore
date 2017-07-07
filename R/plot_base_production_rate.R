plot_base_production_rate = function (ss, main=""){
    ggplot(ss, aes(x = start_time%/%60, y = sequence_length_template/duration)) +
        geom_point(alpha = 0.3)  + xlab("Minute") + theme_bw() + ggtitle(main) +
        theme(plot.title = element_text(hjust = 0.5)) + ylab("Basecalled/duration")
}
