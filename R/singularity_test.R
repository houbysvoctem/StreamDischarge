#testing singularity of matrix
#13. 3. 2018 E
#z webu https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
#returns FALSE, if the matrix is not invertable, which means the matrix is singular

testinv <- function(m) class(try(solve(m),silent=T))=="matrix"
