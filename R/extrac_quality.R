extrac_quality = function(sr, index=1:100,rev=FALSE) {
    l=as(PhredQuality(quality(fq)),"IntegerList")
    hist(sapply(l,median),100)


    mat = as.matrix(l)
    boxplot(mat[,1:1000])

}
