source("RMDS_DNA.R")

fl <- "HIV1_ALL_2018_gag_DNA.fasta"

library(ape)
x <- read.dna("HIV1_ALL_2018_gag_DNA.fasta", "f")
dim(x)
x <- del.colgapsonly(x)
dim(x)

## if you want to analyse only the complete columns of the alignment:
##x <- x[, del.colgapsonly(x, freq.only = TRUE) == 0]
##dim(x)

m0 <- 1000
system.time(S <- select_m_DNA(x, m0)) # 18 sec
d0 <- as.matrix(dist.dna(x[S, ], "N")) # m0 x m0 distance matrix

crit <- opt.m(d0, S)

## asses convergence:
plot(crit, type = "b")

m <- length(crit)
if (m < m0) {
    s <- S[1:m]
    d0 <- d0[s, s]
} else {
    s <- S
}

## reduced MDS
## assume no more than 20 dimensions will be needed
MDS <- cmdscale(d0, k = 20, eig = TRUE)

## assess the number of dimensions:
barplot(MDS$eig[1:20])
## q = 5 seems enough

vars <- MDS$eig
## drop the values <= 0:
vars <- vars[vars > 0]
varpct <- 100 * vars / sum(vars)
varpct[1:5]
## [1] 47.690447 25.009822 12.197244  7.746209  1.765797
cumsum(varpct)[1:5]
## [1] 47.69045 72.70027 84.89751 92.64372 94.40952

q <- 5
mds <- MDS$points[, 1:q, drop = FALSE] # use 'drop = FALSE' in case q = 1

## compute the matrices dB and A to use Gower's (1968) formula
tmds <- t(mds)
B <- mds %*% tmds
dB <- diag(B)
A <- 0.5 * solve(tmds %*% mds) %*% tmds

## prepare the final projections PROJ:
n <- nrow(x) # total number of sequences
PROJ <- numeric(n * q)
dim(PROJ) <- c(n, q)

## the reference sequences used for the reduced MDS above:
mat <- x[s, ]

for (i in (1:n)[-s]) {
    vec <- x[i, ]
    d <- fdist(vec, mat)
    PROJ[i, ] <- A %*% (dB - d^2)
}
PROJ[s, ] <- mds

plot(PROJ) # the first two dimensions
pairs(PROJ)
