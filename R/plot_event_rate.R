plot_event_rate = function (ss, main=""){
    ggplot(ss, aes(x = start_time%/%60, y = num_events/duration)) +
        geom_point(alpha = 0.3)  + xlab("Minute") +
        theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + ggtitle(main)

}
