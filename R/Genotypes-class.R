
# Constructors ----

Genotypes <- function(
    ref = NA_character_, het = NA_character_, alt = NA_character_,
    suffix = c(ref="REF", het="HET", alt="ALT")){
    return(new(
        "Genotypes",
        ref = ref, het = het, alt = alt,
        suffix = suffix[c("ref", "het", "alt")]
    ))
}

# Accessor ----

# genos ----

setMethod("genos", "Genotypes", function(x) c(x@ref, x@het, x@alt))

# ref ----

setMethod("ref", "Genotypes", function(x) x@ref)

setReplaceMethod(
    "ref", c("Genotypes", "character"),
    function(x, value){
        x@ref <- value
        validObject(x)
        return(x)
    }
)

# het ----

setMethod("het", "Genotypes", function(x) x@het)

setReplaceMethod(
    "het", c("Genotypes", "character"),
    function(x, value){
        x@het <- value
        validObject(x)
        return(x)
    }
)

# alt ----

setMethod("alt", "Genotypes", function(x) x@alt)

setReplaceMethod(
    "alt", c("Genotypes", "character"),
    function(x, value){
        x@alt <- value
        validObject(x)
        return(x)
    }
)

# carrier ----

setMethod("carrier", "Genotypes", function(x) c(x@het, x@alt))

# suffix ----

setMethod("suffix", "Genotypes", function(x) c(x@suffix))
