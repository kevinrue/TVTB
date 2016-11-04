
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFixedRules"),
    definition = function(.Object, exprs = list(), ..., active = TRUE){

        parentRules <- FilterRules(exprs = exprs, ..., active = active)

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

        parentRules <- FilterRules(exprs = exprs, ..., active = active)

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

        parentRules <- FilterRules(exprs = exprs, ..., active = active)

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

        new.listData <- as.list(do.call(
            c, lapply(filterList, function(x){slot(x, "listData")})))
        new.active <- as.logical(do.call(
            c, lapply(filterList, function(x){slot(x, "active")})))
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

        .Object@type <- as.character(filterType)

        vepKey <- as.character(unique(unlist(lapply(
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

        if (length(vepKey) > 0)
            vepKey <- as.character(na.omit(vepKey))

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

# Setters ----

setReplaceMethod(
    f = "vep", c("VcfVepRules", "character"),
    function(x, value){

        slot(x, "vep") <- value
        validObject(x)
        return(x)
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

setReplaceMethod(
    f = "vep", c("VcfFilterRules", "character"),
    function(x, value){

        slot(x, "vep") <- value
        validObject(x)
        return(x)
    }
)

# eval ----

setMethod(
    f = "eval",
    signature = c("VcfFixedRules", "VCF"),
    definition = function(expr, envir, enclos){
        # Apply filters to the fixed slot
        eval(expr = expr, envir = fixed(envir), enclos = enclos)
    }
)

setMethod(
    f = "eval",
    signature = c("VcfInfoRules", "VCF"),
    definition = function(expr, envir, enclos){
        # Apply filters to the info slot
        eval(expr = expr, envir = info(envir), enclos = enclos)
    }
)

# TODO info.key taken from TVTBParam
.evalVepFilter <- function(expr, envir, enclos){
    # If empty filter
    if (length(expr) == 0)
        return(rep(TRUE, length(envir)))
    # This filter requires unique rownames in the ExpandedVCF
    if (any(duplicated(rownames(envir))))
        rownames(envir) <- paste(
            rownames(envir),
            mcols(envir)[,"ALT"],
            sep = "_")
    if (any(duplicated(rownames(envir))))
        stop("<rownames(x)>_<ALT> is not unique")
    # In case users called expand(..., row.names = FALSE)
    if (is.null(rownames(envir)))
        rownames(envir) <- 1:length(envir)
    csq <- parseCSQToGRanges(
        x = envir, VCFRowID = rownames(envir), info.key = vep(expr))
    # Apply filters to the VEP predictions
    csqBoolean <- eval(expr = expr, envir = mcols(csq), enclos = enclos)
    # Index of variants with 1+ prediction passing the filter
    vcfPass <- unique(mcols(csq)[csqBoolean, "VCFRowID", drop = TRUE])
    # Mark each variant with 1+ prediction passing the filter
    vcfBoolean <- rep(FALSE, length(envir))
    vcfBoolean[vcfPass] <- TRUE
    return(vcfBoolean)
}

setMethod(
    f = "eval",
    signature = c("VcfVepRules", "VCF"),
    definition = function(expr, envir, enclos){
        .evalVepFilter(expr = expr, envir = envir, enclos = enclos)
    }
)

setMethod(
    f = "eval",
    signature = c("VcfFilterRules", "VCF"),
    definition = function(expr, envir, enclos){
        if (length(expr) == 0)
            return(rep(TRUE, length(envir)))

        # Apply each type of filter to the relevant slot
        resultList <- lapply(
            X = c("fixed", "info", "vep"),
            FUN = function(validType, expr, envir){
                filters <- expr[which(type(expr) == validType)]

                switch(
                    validType,
                    fixed = eval(
                        expr = filters, envir = fixed(envir), enclos = enclos),
                    info = eval(
                        expr = filters, envir = info(envir), enclos = enclos),
                    vep = .evalVepFilter(
                        expr = filters, envir = envir, enclos = enclos))
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

# After subsetting/replacing, down-type a VcfFilterRules to one of the
# more specialised classes
.dropVcfFilterRules <- function(x){

    filterType <- unique(type(x))

    if (length(filterType) > 1)
        return(x)
    else if (length(filterType) == 0)
        return(VcfFilterRules()) # a default empty filter

    listData <- slot(x, "listData")

    res <- switch (
        filterType,
        fixed = VcfFixedRules(exprs = listData, active = slot(x, "active")),
        info = VcfInfoRules(exprs = listData, active = slot(x, "active")),
        vep = VcfVepRules(
            exprs = listData,
            active = slot(x, "active"),
            vep = slot(x, "vep")),
        stop("Invalid type")
    )

    res@active <- slot(x, "active")
    res@elementType <- slot(x, "elementType")
    res@elementMetadata <- slot(x, "elementMetadata")
    res@metadata <- slot(x, "metadata")

    return(res)
}

.namesToIndex <- function(indexNames, object){
    i <- match(indexNames, names(object))
    iMissing <- which(is.na(i))
    if (length(iMissing) > 0)
        stop(
            "names not found: ",
            paste(indexNames[iMissing], collapse = ", "))
    return(i)
}

setMethod(
    f = "[",
    signature = c(x="VcfFilterRules", i="ANY", j="ANY", drop="ANY"),
    definition = function(x, i, j, ..., drop = TRUE){

        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        if (!missing(i)) {
            if (is.character(i)){
                i <- .namesToIndex(indexNames = i, object = x)
            }

            # callNextMethod will cause validity check, filter type before
            x@type <- slot(x, "type")[i]
            x <- callNextMethod(x, i)
        }

        # Down-type to a more specialised class if requested and applicable
        if (drop)
            x <- .dropVcfFilterRules(x)

        return(x)
    }
)

# [<- ----

setReplaceMethod(
    f = "[", c("VcfFixedRules", "ANY", "ANY", "ANY"),
    function(x, i, j, value){
        # Check that no j is provided
        if (!missing(j))
            stop("invalid subsetting")
        # Convert i from character to index
        if (!missing(i)) {
            if (is.character(i)){
                i <- .namesToIndex(indexNames = i, object = x)
            }

            # If NULL, remove rule(s)
            if (is.null(value)){
                x@active <- x@active[-i]
                x@listData <- x@listData[-i]
            } else {
                stopifnot(class(value) == class(x))
                if (length(i) != length(value))
                    stop(paste(
                        "number of items to replace is not a multiple of",
                        "replacement length"
                    ))
                # Transfer the active status(es)
                active(x)[i] <- active(value)
                # Transfer the name(s)
                names(x)[i] <- names(value)
                # Tranfer the expression(s)
                x@listData[i] <- value@listData
                # maybe elementMetadata too
            }
        }
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "[", c("VcfInfoRules", "ANY", "ANY", "ANY"),
    function(x, i, j, value){
        # Check that no j is provided
        if (!missing(j))
            stop("invalid subsetting")
        # Convert i from character to index
        if (!missing(i)) {
            if (is.character(i)){
                i <- .namesToIndex(indexNames = i, object = x)
            }

            # If NULL, remove rule(s)
            if (is.null(value)){
                x@active <- x@active[-i]
                x@listData <- x@listData[-i]
            } else {
                stopifnot(class(value) == class(x))
                if (length(i) != length(value))
                    stop(paste(
                        "number of items to replace is not a multiple of",
                        "replacement length"
                    ))
                # Transfer the active status(es)
                active(x)[i] <- active(value)
                # Transfer the name(s)
                names(x)[i] <- names(value)
                # Tranfer the expression(s)
                x@listData[i] <- value@listData
                # maybe elementMetadata too
            }
        }
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "[", c("VcfVepRules", "ANY", "ANY", "ANY"),
    function(x, i, j, value){
        # Check that no j is provided
        if (!missing(j))
            stop("invalid subsetting")
        # Convert i from character to index
        if (!missing(i)) {
            if (is.character(i)){
                i <- .namesToIndex(indexNames = i, object = x)
            }

            # If NULL, remove rule(s)
            if (is.null(value)){
                x@active <- x@active[-i]
                x@listData <- x@listData[-i]
            } else {
                stopifnot(class(value) == class(x))
                if (length(i) != length(value))
                    stop(paste(
                        "number of items to replace is not a multiple of",
                        "replacement length"
                    ))
                if (vep(value) != vep(x))
                    stop(
                        "Incompatible vep slots: ",
                        vep(x), ", ", vep(value))
                # Transfer the active status(es)
                active(x)[i] <- active(value)
                # Transfer the name(s)
                names(x)[i] <- names(value)
                # Tranfer the expression(s)
                x@listData[i] <- value@listData
                # maybe elementMetadata too
            }
        }
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "[", c("VcfFilterRules", "ANY", "ANY", "ANY"),
    function(x, i, j, value, drop = TRUE){

        # Check that no j is provided
        if (!missing(j))
            stop("invalid subsetting")
        # Check that value is a valid class
        validValueClass <- c(
            "NULL",
            "VcfFixedRules",
            "VcfInfoRules",
            "VcfVepRules",
            "VcfFilterRules")
        if (!class(value) %in% validValueClass)
            stop(
                "class(value) must be one of",
                paste(validValueClass, collapse = ", "))

        # Convert i from character to index
        if (!missing(i)) {
            if (is.character(i)){
                i <- .namesToIndex(indexNames = i, object = x)
            }

            # If NULL, remove rule(s)
            if (is.null(value)){
                x@active <- x@active[-i]
                x@type <- x@type[-i]
                x@listData <- x@listData[-i]
                # maybe elementMetadata too
            } else {
                if (length(i) != length(value))
                    stop(paste(
                        "number of items to replace is not a multiple of",
                        "replacement length"
                    ))

                # Like FilterRules:
                # Only transfer the expression
                # Do not transfer the active state
                # Do not transfer the name of new rules
                x <- switch(
                    class(value),
                    VcfFixedRules = {
                        x@type[i] <- "fixed"
                        x
                    },
                    VcfInfoRules = {
                        x@type[i] <- "info"
                        x
                    },
                    VcfVepRules = {
                        if (vep(value) != vep(x))
                            stop(
                                "Incompatible vep slots: ",
                                vep(x), ", ", vep(value))
                        x@type[i] <- "vep"
                        x
                    },
                    VcfFilterRules = {
                        x@type[i] <- value@type;
                        x
                    }
                )
                # Transfer the active status(es)
                active(x)[i] <- active(value)
                # Transfer the name(s)
                names(x)[i] <- names(value)
                # Tranfer the expression(s)
                x@listData[i] <- value@listData
            }

            # Down-type to a more specialised class if applicable
            if (drop)
                x <- .dropVcfFilterRules(x)
        }
        validObject(x)
        return(x)
    }
)

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
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(x = res)
        return(res)
    }
)

setMethod(
    f = "c",
    signature = c(x="VcfVepRules"),
    definition = function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(x = res)
        return(res)
    }
)

setMethod(
    f = "c",
    signature = c(x="VcfFilterRules"),
    definition = function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(x = res)
        return(res)
    }
)

# append ----

setMethod(
    f = "append",
    signature = c(x="VcfFixedRules", values="FilterRules"),
    definition = function(x, values, after = length(x))
        base::append(x, values, after)
)
setMethod(
    f = "append",
    signature = c(x="VcfInfoRules", values="FilterRules"),
    definition = function(x, values, after = length(x))
        base::append(x, values, after)
)
setMethod(
    f = "append",
    signature = c(x="VcfVepRules", values="FilterRules"),
    definition = function(x, values, after = length(x))
        base::append(x, values, after)
)
setMethod(
    f = "append",
    signature = c(x="VcfFilterRules", values="FilterRules"),
    definition = function(x, values, after = length(x))
        base::append(x, values, after)
)
