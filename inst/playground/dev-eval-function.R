infoVcf <- info(evcf)

colnames(infoVcf)
grep("_AF", colnames(infoVcf), value = T)

AFcut <- 0.05
popMinRatio <- 2/3

f <- function(rd){
    df <- as.data.frame(rd)
    print(df)
    print(colnames(df))
    popFreq <- grep("_AF", colnames(df))
    print(popFreq)
    print(df[,popFreq])
    print(AFcut)
    print(df[,popFreq] > AFcut)
    counts <- rowSums(df[,popFreq] > AFcut)
    print(counts)
    popCutOff <- popMinRatio * length(popFreq)
    testRes <- counts > popCutOff
    print(testRes)
    return(testRes)
}

filters <- FilterRules(f)

str(filters)
filters[[1]]

eval(filters, infoVcf)
