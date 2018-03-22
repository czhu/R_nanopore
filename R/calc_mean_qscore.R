calc_mean_qscore = function(x) {
    ## https://community.nanoporetech.com/posts/difference-in-quality-scor
    ## qscore is in log, mean of qscores doesn't work
    ## mean of error rate works  first scale back and take mean and taking log again

    ## isssue1 :the calculated mean_qscore is lower than the one given in sequece_summary
    ## ~ 0.1 lower tried to find the difference, no luck
    ## I think it's numeric error in the python packages used by ONT for albacore
    ## this implementation returns the same value as the my R implementation
    ## a similar issue is observed here, not exactly the same
    ## https://gigabaseorgigabyte.wordpress.com/2017/07/14/calculated-average-quality-vs-albacore-summary/
    log10(BiocGenerics::mean(10^(x/(-10))))*(-10)
}
