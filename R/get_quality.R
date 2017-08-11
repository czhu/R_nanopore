get_quality = function(x){
    as(PhredQuality(quality(x)),"IntegerList")
}
