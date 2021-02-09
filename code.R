####################################################################
## R code to accompany manuscript:
## Reduced multidimensional scaling
## by Emmanuel Paradis

## Algorithm 1 (four implementations)

## Version from Paradis (2018) J. Comput. Graph. Statist.
## X: matrix of data frame
## m: number of observations to select
## FUN: function (with 2 arguments) to compute distances between two observations
select_m_2 <- function(X, m, FUN = function(a, b) sqrt(sum((a - b)^2)))
{
    if (m < 2) stop("'m' must be at least 2")
    res <- numeric(m)
    i <- 1 # sample(nrow(D), size = 1)
    o <- 1
    res[o] <- i
    repeat {
        ##cat(o, " ")
        D <- apply(X, 1, function(x) FUN(x, X[i, ]))
        D[res[1:o]] <- NA_real_
        k <- which.max(D)
        o <- o + 1
        res[o] <- k
        if (o == m) break

        D <- apply(X, 1, function(x) FUN(x, X[k, ]))
        D[res[1:o]] <- NA_real_

        ## two alternatives:
        #i <- which.min(abs(D - median(D, na.rm = TRUE)))
        i <- sample(which(!is.na(D)), size = 1)

        o <- o + 1
        res[o] <- i
        if (o == m) break
    }
    res
}

## this version uses the Euclidean distance with a (way faster) vectorized code
select_m_2BIS <- function(X, m = 100)
{
    mysample <- function(n) ceiling(runif(1, 0, n)) # a bit faster than: sample.int(n, size = 1), even if replace = TRUE
    if (m < 2) stop("'m' must be at least 2")
    res <- numeric(m)
    XT <- t(X)
    ## randomize the 1st observation picked
    i <- mysample(nrow(X))
    o <- 1
    res[o] <- i
    repeat {
        ##cat("\r", o)
        D <- colSums((XT - X[i, ])^2)
        D[res[1:o]] <- NA_real_
        k <- which.max(D)
        o <- o + 1
        res[o] <- k
        if (o == m) break

        D <- colSums((XT - X[k, ])^2)
        D[res[1:o]] <- NA_real_

        ## two alternatives:
        #i <- which.min(abs(D - median(D, na.rm = TRUE)))
        i <- sample(which(!is.na(D)), size = 1)

        o <- o + 1
        res[o] <- i
        if (o == m) break
    }
    res
}

## if there is a running installation of R:
system("R CMD SHLIB rmds.c")
dyn.load("rmds.so") # under Windows: substitute .so by .dll
select_m_Call <- function(X, m = 100) .Call("select_Euclidean", X, m)
select_m_DNA <- function(X, m = 100) .Call("select_DNAbin", X, m)

## Function to perform 2-D projection
##
## MODIFIED FROM Paradis 2018
## X: the data (matrix or data frame)
## m: number of observations to select
## k: number of dimensions
## FUN: the function used to compute the distance matrix
foo.projectMDS2d.PORT <- function(X, m = 100, k = 2, FUN = dist)
{
    projectMDS2d.PORT <- function(x) {
        D <- 0
        ## vectorizing this loop and the one below is not faster
        for (i in 1:ncol(X)) D <- D + (x[i] - xs[, i])^2
        D <- sqrt(D)
        ##D <- sqrt(colSums((x - xsT)^2))
        f <- function(p) {
            d <- 0
            for (i in 1:k) d <- d + (p[i] - z[, i])^2
            d <- sqrt(d)
            ##d <- sqrt(colSums((p - zT)^2))
            sum((D - d)^2)
        }
        nlminb(rep(0, k), f)$par
    }
    n <- nrow(X)
    ## is <- sample(n, size = m)
    is <- select_m_Call(X, m)
    xs <- X[is, ]
    z <- cmdscale(FUN(xs), k)
    ##zT <- t(z)
    xp <- X[-is, ] # the 'unknown' points
    ##xsT <- t(xs)
    PROJ <- t(apply(xp, 1, projectMDS2d.PORT))
    OUT <- matrix(NA, n, k)
    OUT[is, ] <- z
    OUT[(1:n)[-is], ] <- PROJ
    OUT
}

