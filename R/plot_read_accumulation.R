plot_read_accumulation = function(ss,main=""){
    readAccumulation <- group_by(ss, minute = start_time%/%60) %>%
        summarise(new_reads = n()) %>% mutate(accumulation = order_by(minute,cumsum(new_reads)))
    ggplot(readAccumulation, aes(x = minute/60, y = accumulation)) +
        geom_point() + xlab("hour") + ylab("reads produced") + theme_bw() + ggtitle(main) +
        theme(plot.title = element_text(hjust = 0.5))
}
