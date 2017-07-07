plot_active_channels = function (ss, main="") {
    ## function from IONISR
    startEndSummary <- mutate(ss, first = start_time%/%60,
        last = (start_time + duration)%/%60)
    tab <- data_frame(minute = unlist(apply(startEndSummary,
        1, function(x) {
            x["first"]:x["last"]
        }))) %>% count(minute)
    ggplot(tab, aes(x = minute/60, y = n)) + geom_point() + xlab("Hour") +
        ylab("Active channels") + theme_bw() + ggtitle(main) +
        theme(plot.title = element_text(hjust = 0.5))
}
