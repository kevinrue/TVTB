vff <- VcfFixedFilter(name = "FILTER", condition = "==", value = "PASS")

vif <- VcfInfoFilter(name = "MAF", condition = "<=", value = 1e-6)

vvf <- VcfVepFilter(
    name = "IMPACT", condition = "%in%", value = c("MODERATE", "HIGH"))

vl <- list(vff, vif, vvf)
vlBad <- list(vff, vif, vvf, 2)

which(sapply(X = vl, FUN = "class") == "VcfFixedFilter")
sapply(X = vl, FUN = "class")
sapply(X = vl, FUN = "typeof")
sapply(X = vl, FUN = "name")
