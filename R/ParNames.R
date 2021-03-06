mvnormParNames <- function(iN, ScalarParameters = FALSE) {
    if (ScalarParameters) {

        foo = c("mu", "sigma", "rho")

    } else {

        foo = c(paste("mu", 1:iN, sep = ""), paste("sigma", 1:iN, sep = ""))
        baz = RhoNames(iN)
        foo = c(foo, baz)

    }
    return(foo)
}

mvtParNames <- function(iN, ScalarParameters = FALSE) {

    if (ScalarParameters) {

        foo = c("mu", "phi", "rho", "nu")

    } else {

        foo = c(paste("mu", 1:iN, sep = ""), paste("phi", 1:iN, sep = ""))
        baz = RhoNames(iN)
        foo = c(foo, baz, "nu")

    }

    return(foo)
}

RhoNames <- function(iN) {
    baz = numeric(iN * (iN - 1L)/2L)
    iC = 1L
    for (i in 1:iN) {
        for (j in i:iN) {
            if (i != j) {
                baz[iC] = paste("rho", i, j, sep = "")
                iC = iC + 1L
            }
        }
    }
    return(baz)
}

FullNamesCoefMulti <- function(iN, Dist, CoefName, ScalarParameters) {
    if (ScalarParameters) {

        if (Dist == "mvnorm")
            vNames = c("mu", "sigma", "rho")
        if (Dist == "mvt")
            vNames = c("mu", "phi", "rho", "nu")

    } else {

        vNames = FullNamesMulti(iN, Dist)

    }

    vOut = paste(CoefName, vNames, sep = ".")
    return(vOut)

}

FullNamesMulti <- function(iN, Dist) {
    vRhoNames = RhoNames(iN)
    if (Dist == "mvnorm")
        vNames = c(paste("mu", 1:iN, sep = ""), paste("sigma", 1:iN, sep = ""), vRhoNames)
    if (Dist == "mvt")
        vNames = c(paste("mu", 1:iN, sep = ""), paste("phi", 1:iN, sep = ""), vRhoNames, "nu")
    return(vNames)
}

FullNamesUni <- function(Dist) {
    vNames = c("location", "scale", "skewness", "shape", "shape2")

    if (Dist == "norm")
        vNames = vNames[c(1L, 2L)]
	if (Dist == "lnorm") 
        vNames = vNames[c(1L, 2L)]
    if (Dist == "snorm")
        vNames = vNames[c(1L, 2L, 3L)]
    if (Dist == "std")
        vNames = vNames[c(1L, 2L, 4L)]
    if (Dist == "sstd")
        vNames = vNames[c(1L, 2L, 3L, 4L)]
    if (Dist == "ast")
        vNames = vNames
    if (Dist == "ast1")
        vNames = vNames[c(1L, 2L, 3L, 4L)]
    if (Dist == "ald")
        vNames = vNames[c(1L, 2L, 3L)]
    if (Dist == "ghskt")
      vNames = vNames[c(1L, 2L, 3L, 4L)]
    if (Dist == "poi")
        vNames = vNames[1L]
    if (Dist == "ber")
        vNames = vNames[1L]
    if (Dist == "gamma")
        vNames = vNames[c(2L, 4L)]
	if (Dist == "gumbel")
        vNames = vNames[c(2L, 4L)]
	if (Dist == "weibull")
        vNames = vNames[c(2L, 4L)]
    if (Dist == "exp")
        vNames = vNames[1L]
    if (Dist == "beta")
        vNames = vNames[c(2L, 4L)]
    if (Dist == "negbin")
      vNames = vNames[c(1L, 2L)]
    if (Dist == "skellam")
      vNames = vNames[c(1L, 2L)]
    return(vNames)
}

getParNamesUni <- function(object) {
    Dist = getDist(object)
    parNames = FullNamesUni(Dist)
    return(parNames)
}
getParNamesMulti <- function(object) {
    Dist = getDist(object)
    iN = object@ModelInfo$iN
    parNames = FullNamesMulti(iN, Dist)
    return(parNames)
}

getParNames <- function(object) {
    if (is(object, "uGASFit") | is(object, "uGASSim"))
        parNames = getParNamesUni(object)
    if (is(object, "mGASFit") | is(object, "mGASSim"))
        parNames = getParNamesMulti(object)
    return(parNames)
}

TypeOfParameters <- function(ScalarParameters) {
    if (ScalarParameters)
        return("Scalars")
    if (!ScalarParameters)
        return("Diagonals")
}
