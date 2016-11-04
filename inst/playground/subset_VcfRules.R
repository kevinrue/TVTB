
library(S4Vectors)

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

# Test [[ ----

str(fr[[1]]) # expression
tryCatch(
    str(fr[[1:2]]),
    error = function(e) geterrmessage()
)
# error: double brackets can only extract a single element

# Test [[<- ----

fr[["ruleExpr"]] <- f2 # does not change the name
fr[[2]] <- expression(B < 4)
str(fr)
tryCatch(
    fr[[1:2]] <- c(expression(C == 4), expression(D == 8)),
    error = function(e) geterrmessage()
)
# (above) error: object of type 'closure' is not subsettable

# Test [ ----

str(fr)
str(fr[1:2]) # FilterRules
str(fr[names(fr)]) # FilterRules

# Test [<- ----

str(fr)
active(fr) <- c(TRUE, FALSE)
fr[2] <- FilterRules(exprs = list(a = expression(E == 5)), active = TRUE)
# The name of the new rule is not transferred
# The active state of the new rule is not transferred
fr[1:2] <- FilterRules(exprs = list(
     a = expression(F == 5),
     b = f))

tryCatch(
    fr[1:2] <- FilterRules(exprs = list(
        a = expression(G == 5))),
    error = function(err) geterrmessage()
)
# The length of value must match the number of elements to replace

tryCatch(
    fr[1] <- NULL,
    error = function(err) geterrmessage()
)
# NULL is not allowed (deactivate the rule instead!)
# For the sake of TVTB I would say that NULL should remove the rule



tryCatch(
    fr[1:2] <- expression(A==3),
    error = function(e) geterrmessage()
)
# error: "'value' must be a FilterRules object (or coercible to a FilterRules object)"


