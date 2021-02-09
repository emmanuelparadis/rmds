## if needs compilation:
## system("R CMD SHLIB rmds.c")
dyn.load("rmds.so")

source("common_files.R")

## function to select 'm' observations
select_m <- function(X, m = 100) .Call("select_Euclidean", X, m)

## assume 1000 or less observations will be enough:
m0 <- 1000L
## Besides MDS with n > 1000 points is difficult.

## choose one of the two below:
## S <- sample.int(n, size = m0)
S <- select_m(x, m = m0)

d0 <- as.matrix(dist(x[S, ])) # m0 x m0 distance matrix

crit <- opt.m(d0, S)

m <- length(crit)
s <- S[1:m]
d0 <- d0[s, s]

## reduced MDS
## assume not more than 20 dimensions will be needed
MDS <- cmdscale(d0, k = 20, eig = TRUE)

## assess the number of dimensions:
barplot(MDS$eig[1:20])
q <- 2 # select the appropriate value from the previous barplot
## rerun cmdscale() if q > 20
mds <- MDS$points[, 1:q, drop = FALSE] # use 'drop = FALSE' in case q = 1

## prepare to use Gower's (1968) formula
tmds <- t(mds)
B <- mds %*% tmds
dB <- diag(B)
A <- 0.5 * solve(tmds %*% mds) %*% tmds

n <- # number of observations
PROJ <- numeric(n * q)
dim(PROJ) <- c(n, q)

## the "reduced" (or reference) sample:
xs <- x[s, ]

for (i in (1:n)[-s]) {
    tmp <- matrix(x[i, ], m, ncol(x), byrow = TRUE)
    d <- rowSums((xs - tmp)^2) # assumes Euclidean distance
    ## d <- sqrt(d)
    PROJ[i, ] <- A %*% (dB - d) # avoid squaring d
}
PROJ[s, ] <- mds

plot(PROJ)
