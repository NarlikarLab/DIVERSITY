n <- 1500

kcount <- 3

#set.seed(10)
set.seed(5)
#set.seed(15)

k <- rdirichlet(1, rep(2, kcount))
ksize <- unlist(Map(function(x) round(x*n), head(as.vector(k), kcount - 1)))
ksize <- c(ksize, n - sum(ksize))
print(ksize)
l <- 100
d <- 8

## zoops = 1
## pos <- Map(function(x) sample(c(1:(l - d + 1), l + 1), size = x, replace = TRUE), ksize)

##zoops = 0
## pos <- Map(function(x) sample(1:(l - d + 1), size = x, replace = TRUE), ksize)

dval <- 0.3

##print(pos)
