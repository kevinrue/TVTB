
# Constructors ----

.VcfBasicRules <- function(filterClass, exprs = list(), ..., active = TRUE){

    # Use the parent constructor to parse exprs and ...
    parentRules <- FilterRules(exprs, ..., active = active)

    # Initalise a prototype object of the class
    newRules <- new(filterClass)

    # Transfer the slot from the parent object
    newRules@active <- parentRules@active
    newRules@listData <- parentRules@listData
    newRules@elementType <- parentRules@elementType
    newRules@elementMetadata <- parentRules@elementMetadata
    newRules@metadata <- parentRules@metadata

    return(newRules)
}

# VcfFixedRules ----

VcfFixedRules <- function(exprs = list(), ..., active = TRUE){

    # Populate the basic slots (same as FilterRules)
    newRules <- .VcfBasicRules("VcfFixedRules", exprs, ..., active = active)

    # Check validity & return
    validObject(newRules)
    return(newRules)
}

# VcfInfoRules ----

VcfInfoRules <- function(exprs = list(), ..., active = TRUE){

    # Populate the basic slots (same as FilterRules)
    newRules <- .VcfBasicRules("VcfInfoRules", exprs, ..., active = active)

    # Check validity & return
    validObject(newRules)
    return(newRules)
}

# VcfVepRules ----

VcfVepRules <- function(exprs = list(), ..., active = TRUE, vep = "CSQ"){

    # First populate the basic slots
    newRules <- .VcfBasicRules("VcfVepRules", exprs, ..., active = active)

    # Then populated the vep slot
    newRules@vep <- vep

    # Check validity & return
    validObject(newRules)
    return(newRules)
}

# VcfFilterRules ----

VcfFilterRules <- function(...){

    filterList <- list(...)

    # Concatenate the basic slots of the objects supplied
    new.active <- as.logical(do.call(
        c, lapply(filterList, function(x){x@active})))
    new.listData <- as.list(do.call(
        c, lapply(filterList, function(x){x@listData})))
    new.listData <- as.list(do.call(
        c, lapply(filterList, function(x){x@listData})))
    # Do not transfer elementType: let the constructor set it
    new.elementMetadata <- do.call(
        c, lapply(filterList, function(x){x@elementMetadata}))
    # Do not transfer metadata: it is object-specific

    # Identify the VEP key(s) defined across the rules
    # NOTE: use accessor that handles classes without vep slot
    vepKey <- as.character(unique(unlist(
        lapply(filterList, function(x) na.omit(vep(x)))
    )))

    # Concatenate the type of each rule supplied
    # Same as above: use the accessor that handles classes without slot
    filterType <- as.character(unlist(lapply(filterList, type)))

    # Initialise an empty object
    newRules <- new("VcfFilterRules")
    # Populate the basic slots of the object
    newRules@listData <- new.listData
    newRules@active <- new.active
    newRules@elementMetadata <- new.elementMetadata

    # Populate the vep slot (default value is NA)
    if (length(vepKey) > 0){
        newRules@vep <- vepKey
    }

    # Populate the type slot
    newRules@type <- filterType

    # Check validity & return
    validObject(newRules)
    return(newRules)
}

# Accessors ----

# type ----
# Returned as a named vector

setMethod("type", "VcfFilterRules", function(x){
    t <- x@type
    names(t) <- names(x)
    return(t)
})

setMethod("type", "VcfVepRules", function(x){
    t <- rep("vep", length(x))
    names(t) <- names(x)
    return(t)
})

setMethod("type", "VcfInfoRules", function(x){
    t <- rep("info", length(x))
    names(t) <- names(x)
    return(t)
})

setMethod("type", "VcfFixedRules", function(x){
    t <- rep("fixed", length(x))
    names(t) <- names(x)
    return(t)
})

setMethod("type", "FilterRules", function(x){
    t <- rep("filter", length(x))
    names(t) <- names(x)
    return(t)
})

# vep ----

setMethod("vep", "VcfFilterRules", function(x) x@vep)

setMethod("vep", "VcfVepRules", function(x) x@vep)

setMethod("vep", "VcfInfoRules", function(x) NA_character_)

setMethod("vep", "VcfFixedRules", function(x) NA_character_)

setMethod("vep", "FilterRules", function(x) NA_character_)

# Setters ----

setReplaceMethod(
    "vep", c("VcfVepRules", "character"),
    function(x, value){
        x@vep <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    "vep", c("VcfFilterRules", "character"),
    function(x, value){
        x@vep <- value
        validObject(x)
        return(x)
    }
)

# eval ----

setMethod(
    "eval", c("VcfFixedRules", "VCF"),
    function(expr, envir, enclos){
        # Apply filters to the fixed slot
        return(eval(expr, fixed(envir), enclos))
    }
)

setMethod(
    "eval", c("VcfInfoRules", "VCF"),
    function(expr, envir, enclos){
        # Apply filters to the info slot
        return(eval(expr, info(envir), enclos))
    }
)

# Evaluate VepFilterRules
.evalVepFilter <- function(expr, envir, enclos){

    # If empty filter, return original object
    if (length(expr) == 0){return(rep(TRUE, length(envir)))}

    # This filter requires unique rownames in the VCF
    rownames(envir) <- paste0("R", seq_along(envir))

    # Extract VEP prediction to evaluate filter
    csq <- parseCSQToGRanges(envir, rownames(envir), info.key = vep(expr))

    # Apply filters to the VEP predictions
    csqBoolean <- eval(expr, mcols(csq), enclos)

    # Index of variants with 1+ prediction passing the filter
    vcfPass <- unique(mcols(csq)[csqBoolean, "VCFRowID", drop = TRUE])

    # Mark each variant with 1+ prediction passing the filter
    vcfBoolean <- rep(FALSE, length(envir))
    vcfBoolean[vcfPass] <- TRUE
    return(vcfBoolean)
}

setMethod(
    "eval", c("VcfVepRules", "VCF"),
    function(expr, envir, enclos){
        return(.evalVepFilter(expr, envir, enclos))
    }
)

setMethod(
    "eval", c("VcfFilterRules", "VCF"),
    function(expr, envir, enclos){
        if (length(expr) == 0){return(rep(TRUE, length(envir)))}

        # TODO: consider parallel processing
        # Apply each type of filter to the relevant slot
        # NOTE: probably faster to process all rules of the same type together
        # than dispatching <N> VcfFixedRules independently, for instance
        resultList <- lapply(
            X = c("filter", "fixed", "info", "vep"),
            FUN = function(validType, expr, envir){
                filters <- expr[which(type(expr) == validType)]
                # Use TVTB VCF filter rules dispatch
                eval(filters, envir, enclos)
            },
            expr = expr,
            envir = envir)
        # Combine the boolean DataFrame objects
        resultMatrix <- do.call(cbind, resultList)

        # TRUE if all row is TRUE
        return(apply(resultMatrix, 1, all))
    }
)

# [ ----

# Down-type VcfFilterRules of a single type
.dropVcfFilterRules <- function(x){

    # Empty object cannot be down-typed
    if (length(x) == 0){return(x)}

    # Identify the type(s) of VCF filter rules present in the object
    filterType <- unique(x@type)

    # Multiple types cannot be down-typed, by definition
    if (length(filterType) > 1){return(x)}

    # Down-type to appropriate class
    newRules <- switch(
        filterType,
        filter = FilterRules(exprs = x@listData, active = x@active),
        fixed = VcfFixedRules(exprs = x@listData, active = x@active),
        info = VcfInfoRules(exprs = x@listData, active = x@active),
        vep = VcfVepRules(exprs = x@listData, active = x@active, vep = x@vep),
        stop("Invalid filter type")
    )

    newRules@active <- slot(x, "active")
    newRules@elementType <- slot(x, "elementType")
    newRules@elementMetadata <- slot(x, "elementMetadata")
    newRules@metadata <- slot(x, "metadata")

    return(newRules)
}

setMethod(
    "vertical_slot_names", "VcfFilterRules",
    function(x) c("type", callNextMethod()))

# Override method to subset the additional "type" slot
setMethod(
    "[", c("VcfFilterRules", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE){

        x <- callNextMethod() # x, i, j, ..., drop = drop

        # Down-type to a more specialised class if requested and applicable
        if (drop){
            x <- .dropVcfFilterRules(x)
        }

        return(x)
    }
)

# [<- ----

# Substitution:

setReplaceMethod(
    "[", c("VcfVepRules", "numeric", "missing", "VcfVepRules"),
    function(x, i, j, ..., value){

        # Refuse conflict of vep slots
        stopifnot(value@vep == x@vep)

        x <- callNextMethod()
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    "[", c("VcfFilterRules", "numeric", "missing", "VcfFilterRules"),
    function(x, i, j, ..., value){

        # Refuse conflict of vep slots
        # NOTE: NA is expected when fixed/info rules are coerced to VCF rule
        stopifnot(value@vep %in% c(x@vep, NA_character_))

        # Replace inherited and parallel slots
        x <- callNextMethod()

        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    "[", c("VcfFilterRules", "numeric", "missing", "VcfFixedRules"),
    function(x, i, j, ..., value){
        x[i] <- VcfFilterRules(value)
        return(x)
    }
)

setReplaceMethod(
    "[", c("VcfFilterRules", "numeric", "missing", "VcfInfoRules"),
    function(x, i, j, ..., value){
        x[i] <- VcfFilterRules(value)
        return(x)

    }
)

setReplaceMethod(
    "[", c("VcfFilterRules", "numeric", "missing", "VcfVepRules"),
    function(x, i, j, ..., value){
        x[i] <- VcfFilterRules(value)
        return(x)
    }
)


# setAs ----

setAs("VcfFixedRules", "VcfFilterRules", function(from) VcfFilterRules(from))

setAs("VcfInfoRules", "VcfFilterRules", function(from) VcfFilterRules(from))

setAs("VcfVepRules", "VcfFilterRules", function(from) VcfFilterRules(from))

setAs("FilterRules", "VcfFilterRules", function(from) VcfFilterRules(from))

### Coercions from SimpleList to VcfFixedRules/VcfInfoRules/VcfVepRules work
### out-of-the-box but silently return a broken object! The problem is that
### these coercions are performed by dummy coercion methods that are
### automatically defined by the methods package and that often do the wrong
### thing (like here). Furthermore, they don't bother to validate the object
### they return. So we overwrite them. The reason it's important that these
### coercions work properly is that, starting with S4Vectors >= 0.17.4, the
### "setListElement" method for SimpleList objects was replaced with a method
### for List objects that uses these coercions internally.
setAs("SimpleList", "VcfFixedRules", function(from) VcfFixedRules(from))
setAs("SimpleList", "VcfInfoRules", function(from) VcfInfoRules(from))
setAs("SimpleList", "VcfVepRules", function(from) VcfVepRules(from))

# setListElement ---

### Because the "type" of the filters stored in a SimpleList or FilterRules
### instance is not known, it doesn't seem to be possible to define a
### coercion method from SimpleList to VcfFilterRules. Therefore, in order
### to make `[[<-` work again on VcfFilterRules objects, we implement a
### "setListElement" method for them. Note that this is the method that was
### defined for SimpleList objects prior to S4Vectors 0.17.4 so it actually
### restores the behavior of `[[<-` on VcfFilterRules objects, which,
### unfortunately, was (and is still) broken in case of appending, that is,
### when doing 'x[["new_name"]] <- value'. Furthermore, it also seems to do
### the wrong thing in situations like 'x[[i]] <- y[[i]]' when 'x' and 'y'
### are both VcfFilterRules objects and the type of the i-th element in 'x'
### is not the same as the type of the i-th element in 'y'. These problems
### are quite deep and perhaps reveal a possible flaw in the current design
### of the VcfFilterRules class.
setMethod(
    "setListElement", "VcfFilterRules",
    function(x, i, value){
        x@listData[[i]] <- value
        return(x)
    }
)

# Combine ----

setMethod(
    "c", "VcfFixedRules",
    function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(res)
        return(res)
    }
)

setMethod(
    "c", "VcfInfoRules",
    function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(res)
        return(res)
    }
)

setMethod(
    "c", "VcfVepRules",
    function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(res)
        return(res)
    }
)

setMethod(
    "c", "VcfFilterRules",
    function(x, ..., recursive = FALSE){
        res <- VcfFilterRules(x, ...)
        res <- .dropVcfFilterRules(res)
        return(res)
    }
)

# append ----

setMethod(
    "append", c("VcfFixedRules", "FilterRules"),
    function(x, values, after = length(x))
        base::append(x, values, after)
)
setMethod(
    "append", c("VcfInfoRules", "FilterRules"),
    function(x, values, after = length(x))
        base::append(x, values, after)
)
setMethod(
    "append", c("VcfVepRules", "FilterRules"),
    function(x, values, after = length(x))
        base::append(x, values, after)
)
