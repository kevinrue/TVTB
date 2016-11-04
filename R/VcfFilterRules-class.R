
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFixedRules"),
    definition = function(.Object, exprs = list(), ..., active = TRUE){

        parentRules <- FilterRules(exprs = exprs, ..., active = TRUE)

        .Object@listData <- slot(parentRules, "listData")
        .Object@active <- slot(parentRules, "active")
        .Object@elementType <- slot(parentRules, "elementType")
        .Object@elementMetadata <- slot(parentRules, "elementMetadata")
        .Object@metadata <- slot(parentRules, "metadata")

        validObject(.Object)
        return(.Object)
    })

VcfFixedRules <- function(exprs = list(), ..., active = TRUE){
    new(
        Class = "VcfFixedRules",
        exprs = exprs, ..., active = active)
}

setMethod(
    f = "initialize",
    signature = c("VcfInfoRules"),
    definition = function(.Object, exprs = list(), ..., active = TRUE){

        parentRules <- FilterRules(exprs = exprs, ..., active = TRUE)

        .Object@listData <- slot(parentRules, "listData")
        .Object@active <- slot(parentRules, "active")
        .Object@elementType <- slot(parentRules, "elementType")
        .Object@elementMetadata <- slot(parentRules, "elementMetadata")
        .Object@metadata <- slot(parentRules, "metadata")

        validObject(.Object)
        return(.Object)
    })

VcfInfoRules <- function(exprs = list(), ..., active = TRUE){
    new(
        Class = "VcfInfoRules",
        exprs = exprs, ..., active = active)
}

setMethod(
    f = "initialize",
    signature = c("VcfVepRules"),
    definition = function(
        .Object, exprs = list(), ..., active = TRUE, vep = "CSQ"){

        parentRules <- FilterRules(exprs = exprs, ..., active = TRUE)

        .Object@listData <- slot(parentRules, "listData")
        .Object@active <- slot(parentRules, "active")
        .Object@elementType <- slot(parentRules, "elementType")
        .Object@elementMetadata <- slot(parentRules, "elementMetadata")
        .Object@metadata <- slot(parentRules, "metadata")

        .Object@vep <- vep

        validObject(.Object)
        return(.Object)
    })

VcfVepRules <- function(exprs = list(), ..., active = TRUE, vep = "CSQ"){
    new(
        Class = "VcfVepRules",
        exprs = exprs, ..., active = active, vep = vep)
}

setMethod(
    f = "initialize",
    signature = c("VcfFilterRules"),
    definition = function(.Object, ...){

        filterList <- list(...)

        new.listData <- do.call(
            c, lapply(filterList, function(x){slot(x, "listData")}))
        new.active <- do.call(
            c, lapply(filterList, function(x){slot(x, "active")}))
        new.elementMetadata <- do.call(
            c, lapply(filterList, function(x){slot(x, "elementMetadata")}))

        .Object@listData <- new.listData
        .Object@active <- new.active
        .Object@elementMetadata <- new.elementMetadata

        filterType <- unlist(lapply(
            X = filterList,
            FUN = function(x){
                    switch(
                        class(x),
                        VcfFixedRules = rep("fixed", length(x)),
                        VcfInfoRules = rep("info", length(x)),
                        VcfVepRules = rep("vep", length(x)),
                        VcfFilterRules = slot(x, "type"),
                        stop("Invalid filter")
                    )
                }
            ))

        .Object@type <- filterType

        vepKey <- unique(na.exclude(unlist(lapply(
            X = filterList,
            FUN = function(x){
                switch(
                    class(x),
                    VcfFixedRules = rep(NA_character_, length(x)),
                    VcfInfoRules = rep(NA_character_, length(x)),
                    VcfVepRules = vep(x),
                    VcfFilterRules = slot(x, "vep"),
                    stop("Invalid filter")
                )
            }
        ))))

        if (length(vepKey) > 1)
            stop("All VcfVepFilter objects must have the same vep slot")
        else if (length(vepKey) == 1)
            .Object@vep <- vepKey
        else
            .Object@vep <- NA_character_

        validObject(.Object)
        return(.Object)
    })

# Accessors ----

setMethod(
    f = "type",
    signature = c("VcfFilterRules"),
    definition = function(x){
        # Apply filters to the fixed slot
        slot(x, "type")
    }
)

setMethod(
    f = "vep",
    signature = c("VcfVepRules"),
    definition = function(x){
        # Apply filters to the fixed slot
        slot(x, "vep")
    }
)

setMethod(
    f = "vep",
    signature = c("VcfFilterRules"),
    definition = function(x){
        # Apply filters to the fixed slot
        slot(x, "vep")
    }
)

# eval ----

setMethod(
    f = "eval",
    signature = c("VcfFixedRules", "ExpandedVCF"),
    definition = function(expr, envir){
        # Apply filters to the fixed slot
        eval(expr = expr, envir = fixed(envir))
    }
)

setMethod(
    f = "eval",
    signature = c("VcfInfoRules", "ExpandedVCF"),
    definition = function(expr, envir){
        # Apply filters to the info slot
        eval(expr = expr, envir = info(envir))
    }
)

# TODO info.key taken from TVTBParam
.evalVepFilter <- function(expr, envir){
    csq <- parseCSQToGRanges(
        x = envir, VCFRowID = rownames(envir), info.key = vep(expr))
    # Apply filters to the VEP predictions
    csqBoolean <- eval(expr = expr, envir = csq)
    # Index of variants with 1+ prediction passing the filter
    vcfPass <- mcols(csq)[csqBoolean,"VCFRowID", drop = TRUE]
    # Mark each variant with 1+ prediction passing the filter
    vcfBoolean <- rep(FALSE, length(envir))
    vcfBoolean[vcfPass] <- TRUE
    return(vcfBoolean)
}

setMethod(
    f = "eval",
    signature = c("VcfVepRules", "ExpandedVCF"),
    definition = function(expr, envir){
        .evalVepFilter(expr = expr, envir = envir)
    }
)

setMethod(
    f = "eval",
    signature = c("VcfFilterRules", "ExpandedVCF"),
    definition = function(expr, envir){
        # Apply each type of filter to the relevant slot
        resultList <- lapply(
            X = c("fixed", "info", "vep"),
            FUN = function(validType, expr, envir){
                filters <- expr[which(type(expr) == validType)]

                switch(
                    validType,
                    fixed = eval(expr = filters, envir = fixed(envir)),
                    info = eval(expr = filters, envir = info(envir)),
                    vep = .evalVepFilter(expr = filters, envir = envir))
            },
            expr = expr,
            envir = envir)
        # Combine the boolean DataFrame objects
        resultMatrix <- do.call(cbind, resultList)
        # TRUE if all row is TRUE
        apply(X = resultMatrix, MARGIN = 1, FUN = all)
        }
)

# [ ----

setMethod(
    f = "[",
    signature = c(x="VcfFilterRules", i="ANY", j="ANY", drop="missing"),
    definition = function(x, i, j, ..., drop = FALSE){

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            x@type <- slot(x, "type")[i]
            x <- callNextMethod(x, i)
        }

        x <- .dropVcfFilterRules(x)

        return(x)
    }
)

.dropVcfFilterRules <- function(x){

    filterType <- unique(type(x))

    if (length(filterType) > 1)
        return(x)
    else if (length(filterType) == 0)
        return(VcfFixedRules()) # a default empty filter

    listData <- slot(x, "listData")

    res <- switch (filterType,
        fixed = VcfFixedRules(exprs = listData, active = slot(x, "listData")),
        info = VcfInfoRules(exprs = listData, active = slot(x, "listData")),
        vep = VcfVepRules(
            exprs = listData,
            active = slot(x, "listData"),
            vep = slot(x, "vep")),
        stop("Invalid type")
    )

    res@active <- slot(x, "active")
    res@elementType <- slot(x, "elementType")
    res@elementMetadata <- slot(x, "elementMetadata")
    res@metadata <- slot(x, "metadata")

    return(res)
}

# [[ ----

setMethod(
    f = "[[",
    signature = c(x="VcfFilterRules", i="ANY", j="ANY"),
    definition = function(x, i, j, ...){

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            x <- x[i]
        }

        return(x)
    }
)

setMethod(
    f = "[[<-",
    signature = c(x="VcfFilterRules", i="ANY", j="ANY"),
    definition = function(x, i, j, ..., value) {

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        # Mark new rules with corresponding type
        x@type[[i]] <- switch(
            class(value),
            VcfFixedRules = "fixed",
            VcfInfoRules = "info",
            VcfVepRules = "vep",
            VcfFilterRules = type(value),
            stop("Invalid filter type")
        )

        if (!missing(i)) {
            # New rules are active by default
            x@active[i] <- TRUE
            x@listData[i] <- slot(value, "listData")
            x@elementMetadata[i] <- slot(value, "elementMetadata")
            x@metadata <- slot(value, "metadata")
            names(x)[i] <- names(value)
        }

        x <- .dropVcfFilterRules(x)

        return(x)
    }
)

setMethod(
    f = "[[<-",
    signature = c(x="VcfFixedRules", i="ANY", j="ANY"),
    definition = function(x, i, j, ..., value) {

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            # New rules are active by default
            x <- .replaceVcfBasicFilter(x, i, value)
        }

        return(x)
    }
)

setMethod(
    f = "[[<-",
    signature = c(x="VcfInfoRules", i="ANY", j="ANY"),
    definition = function(x, i, j, ..., value) {

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            # New rules are active by default
            x <- .replaceVcfBasicFilter(x, i, value)
        }

        return(x)
    }
)

setMethod(
    f = "[[<-",
    signature = c(x="VcfVepRules", i="ANY", j="ANY"),
    definition = function(x, i, j, ..., value) {

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            # New rules are active by default
            x <- .replaceVcfBasicFilter(x, i, value)
        }

        return(x)
    }
)


.replaceVcfBasicFilter <- function(x, i, value){

    x@active[i] <- TRUE
    x@listData[[i]] <- value

    return(x)
}

# coerce ----

setMethod(
    f = "coerce",
    signature = c(from="VcfFixedRules", to="VcfFilterRules"),
    definition = function(from, to){
        VcfFilterRules(from)
    }
)

setMethod(
    f = "coerce",
    signature = c(from="VcfInfoRules", to="VcfFilterRules"),
    definition = function(from, to){
        VcfFilterRules(from)
    }
)

setMethod(
    f = "coerce",
    signature = c(from="VcfVepRules", to="VcfFilterRules"),
    definition = function(from, to){
        VcfFilterRules(from)
    }
)

# Combine ----

setMethod(
    f = "c",
    signature = c(x="VcfFixedRules"),
    definition = function(x, ..., recursive = FALSE){
        VcfFilterRules(x, ...)
    }
)

setMethod(
    f = "c",
    signature = c(x="VcfInfoRules"),
    definition = function(x, ..., recursive = FALSE){
        VcfFilterRules(x, ...)
    }
)

setMethod(
    f = "c",
    signature = c(x="VcfVepRules"),
    definition = function(x, ..., recursive = FALSE){
        VcfFilterRules(x, ...)
    }
)

setMethod(
    f = "c",
    signature = c(x="VcfFilterRules"),
    definition = function(x, ..., recursive = FALSE){
        VcfFilterRules(x, ...)
    }
)
