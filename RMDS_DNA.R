## if needs compilation:
## system("R CMD SHLIB rmds.c")
dyn.load("rmds.so")

source("common_files.R")

## function to select 'm' sequences
select_m_DNA <- function(X, m = 100) .Call("select_DNAbin", X, m)

## compute distances from a single sequence (vec) to a matrix of sequences (mat)
fdist <- function(vec, mat) .Call("dist_dna_vector_to_matrix", vec, mat)
