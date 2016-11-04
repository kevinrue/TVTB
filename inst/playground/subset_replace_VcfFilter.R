vcfRules
v123 <- vcfRules[1:3]
v5 <- vcfRules[[5]]
v56 <- vcfRules[5:6]

str(v123)
v123[[2]] <- v5

vep1 <- VcfVepRules(exprs = list(hello = expression(test == 2)))
v56
str(v56)
v56[[2]] <- vep1

slot(vep1, "listData")

v56@listData[2] <- slot(vep1, "listData")
names(v56)[2] <- names(vep1)
