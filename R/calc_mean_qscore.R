calc_mean_qscore = function(x) {
    ## https://community.nanoporetech.com/posts/difference-in-quality-scor
    ## qscore is in log, mean of qscores doesn't work
    ## mean of error rate works  first scale back and take mean and taking log again
    log10(BiocGenerics::mean(10^(x/(-10))))*(-10)
}
