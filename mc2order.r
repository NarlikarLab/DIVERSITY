library(MCMCpack)


getSequence <- function(n ,mat){
    lst <- sample(c(0, 1, 2, 3), size = 2, replace = TRUE)
    i = 2
    while(n > 2){
        i = i + 1
        x <- sample(c(0, 1, 2, 3), size = 1, prob = mat[lst[i - 2] + 1, lst[i - 1] + 1, ])
        lst <- c(lst, x)
        n = n - 1
    }
    print(lst)
}

## nuc <- 4
## mat <- array(0, dim = c(nuc, nuc, nuc))

## for(i in 1:nuc){
##     for(j in 1:nuc){
##         t <- rdirichlet(1, rep(1, nuc))
##         print(t)
##         for(k in 1:nuc)
##             mat[i, j, k] <- t[k]
##     }
## }

## getSequence(5, mat)
