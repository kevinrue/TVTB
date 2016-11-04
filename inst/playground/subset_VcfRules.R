library(TVTB)
library(S4Vectors)
library(VariantAnnotation)

# Sample data ----

# Make a FilterRules (for reference)
f <- function(envir){
    return(envir[,"FILTER"] == "PASS")
}

f2 <- function(env){
    return(env[,"FILTER"] == "FAIL")
}

fr <- FilterRules(exprs = list(
    ruleExpr = expression(QUAL > 20),
    ruleFUN = f
    ))
active(fr)[2] <- FALSE
fr
str(fr)

# Make a VcfFixedRules (specific VCF filter rules)
vfixed <- VcfFixedRules(exprs = list(
    fixedExpr = expression(QUAL > 20),
    fixedFUN = f
))
active(vfixed)[2] <- FALSE
vfixed
str(vfixed)

# Make a VcfInfoRules (specific VCF filter rules)
vinfo <- VcfInfoRules(exprs = list(
    infoExpr = expression(QUAL > 20),
    infoFUN = f
))
active(vinfo)[2] <- FALSE
vinfo
str(vinfo)

# Make a VcfInfoRules (specific VCF filter rules)
vvep <- VcfInfoRules(exprs = list(
    vepExpr = expression(QUAL > 20),
    vepFUN = f
))
active(vvep)[2] <- FALSE
vvep
str(vvep)

# Make a VcfFilterRules (with multiple types of VCF filter rules)
vfilter <- VcfFilterRules(vfixed, vinfo, vvep)
vfilter
str(vfilter)

# Test [[ ----

str(fr[[1]]) # expression
str(fr[[1:2]]) # error: double brackets can only extract a single element

str(vfixed[[1]]) # expression: OK
str(vfixed[[2]]) # FilterClosure: OK
str(vinfo[[1]]) # expression: OK
str(vinfo[[2]]) # FilterClosure: OK
str(vvep[[1]]) # expression: OK
str(vvep[[2]]) # FilterClosure: OK
str(vfilter[[1]]) # expression: OK
str(vfilter[[2]]) # FilterClosure: OK
vfixed[[1:2]] # error: OK
vinfo[[1:2]] # error: OK
vvep[[1:2]] # error: OK
vfilter[[1:2]] # error: OK

# Test [[<- ----

fr[["ruleExpr"]] <- f2 # does not change the name
fr[[2]] <- expression(B < 4)
str(fr)
fr[[1:2]] <- c(expression(C == 4), expression(D == 8)) # error: object of type 'closure' is not subsettable

vfixed[["fixedExpr"]] <- f2
vfixed[[2]] <- expression(NewField > 4)
str(vfixed) # replacement: OK
vfixed[[1:2]] <- c(f2, expression(C == 4)) # error: OK

vinfo[["infoExpr"]] <- f2
vinfo[[2]] <- expression(NewField > 4)
str(vinfo) # replacement: OK
vinfo[[1:2]] <- c(f2, expression(C == 4)) # error: OK

vvep[["vepExpr"]] <- f2
vvep[[2]] <- expression(NewField > 4)
str(vvep) # replacement: OK
vvep[[1:2]] <- c(f2, expression(C == 4)) # error: OK

vfilter[["fixedExpr"]] <- f
vfilter[[2]] <- expression(NewField > 4)
str(vfilter) # replacement: OK
vfilter[[1:2]] <- c(f2, expression(C == 4)) # error: OK

names(vfixed)[1] <- "hello"
names(vfixed)

names(vinfo)[1] <- "hello"
names(vinfo)

names(vvep)[1] <- "hello"
names(vvep)

names(vfilter)[1] <- "hello"
names(vfilter)

# Test [ ----

str(fr)
str(fr[1:2]) # FilterRules
str(fr[names(fr)]) # FilterRules

str(vfixed[1]) # VcfFixedRules: OK
str(vfixed["fixedFUN"]) # VcfFixedRules: OK

str(vinfo[1]) # VcfFixedRules: OK
str(vinfo["infoFUN"]) # VcfFixedRules: OK
