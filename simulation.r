library(MCMCpack)

vals <- c("A", "C", "G", "T")

sampling <- function(dist, n){
  t = Map(function(x) sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = dist[[x]]), 1:length(dist))
  l <- length(dist)
  return(Map(function(x) unlist(Map(function(y) t[[y]][x], 1:l)), 1:n))
}

generateMCmatrix <- function(arch){
    nuc <- 4
    lst <- list()
    for(i1 in 1:arch){
        mat <- array(0, dim = c(nuc, nuc, nuc))

        for(i in 1:nuc){
            for(j in 1:nuc){
                t <- c(0.25, 0.25, 0.25, 0.25)
##                t <- rdirichlet(1, rep(1, nuc))
                for(k in 1:nuc)
                    mat[i, j, k] <- t[k]
            }
        }
        lst <- append(lst, list(mat))
    }
    return(lst)
}

getdst <- function(n, dval){
  ## dval <- runif(1, 0, 0.7)
  ## y <- Map(function(x) rdirichlet(1, rep(dval, 4)), 1:n)

  ## y <- Map(function(x) rdirichlet(1, rep(runif(1, 0, 1), 4)), 1:n)
  y <- Map(function(x) rdirichlet(1, rep(runif(1, 0, 0.1), 4)), 1:n)
##  print(y)
#  print(y)
  for(i in 1:n){
    y[[i]] <- unlist(Map(function(x) x, y[[i]]))
  }
  return(y)
}

generateData <- function(n, k, d, pos, dval){

  sizeK <- length(k)
  sizeG <- unlist(Map(function(x) round(n*k[x]), 1:(sizeK-1)))
  sizeG <- c(sizeG, (n-sum(sizeG)))
  ## dval1 <- runif(1, 0.7, 1)
  ## undst <- Map(function(x) rdirichlet(1, rep(dval1, 4)), 1:d)

  ## undst <- Map(function(x) rdirichlet(1, rep(runif(1, 2, 4), 4)), 1:d)
  ## undst <- Map(function(x) rdirichlet(1, rep(runif(1, 2, 4), 4)), 1:d)
  undst <- Map(function(x) rep(0.25, 4), 1:d)
  for(i in 1:d)
    undst[[i]] <- unlist(Map(function(x) x, undst[[i]]))
  undst <- Map(function(x) undst, 1:sizeK)
  counts <- unlist(Map(function(x) length(x), pos))
  for(i in 1:sizeK){
    dst <- getdst(counts[i], dval)
    p <- 1
    for(j in pos[[i]]){
      undst[[i]][j] <- dst[p]
      p <- p + 1
    }
  }
  dst <- undst
  data <- Map(function(x) sampling(dst[[x]], sizeG[x]), 1:sizeK)
  return(list(data, sizeK, sizeG))
}

shuffle <- function(data, startSites, sizeG, sizeK, n){

    ## positions <- sample(1:n, size = n, replace = FALSE)
    ## print(positions)
    groupList <- unlist(Map(function(x) rep(x, times = sizeG[x]), 1:sizeK))
    groupList <- sample(groupList, size = n, replace = FALSE, prob = NULL)
    groups <- data
    gStart <- startSites
##    groups <- Map(function(x) sample(data[[x]], size = sizeG[x], replace = FALSE, prob = NULL), 1:sizeK)
##    print(groups)
    lst <- list()
    strt <- c()
    i = 1
##    print(data)
    while(i<=n){
        lst <- append(lst, head(groups[[groupList[i]]], 1))
        strt <- c(strt, head(gStart[[groupList[i]]], 1))
        groups[[groupList[i]]] <- tail(groups[[groupList[i]]], -1)
        gStart[[groupList[i]]] <- tail(gStart[[groupList[i]]], -1)
        i <- i + 1
    }
    return(list(groupList, lst, gStart))
}

## shuffle <- function(data, sizeG, sizeK, n){
##   groupList <- unlist(Map(function(x) rep(x, times = sizeG[x]), 1:sizeK))
##   groupList <- sample(groupList, size = n, replace = FALSE, prob = NULL)
##   groups <- Map(function(x) sample(data[[x]], size = sizeG[x], replace = FALSE, prob = NULL), 1:sizeK)

##   lst <- list()
##   i = 1
##   while(i<=n){
##     lst <- append(lst, head(groups[[groupList[i]]], 1))
##     groups[[groupList[i]]] <- tail(groups[[groupList[i]]], -1)
##     i <- i + 1
##   }
##   print(groupList)
##   return(list(groupList, lst))
## }

getMarkovSequence <- function(n ,mat){
    if(n == 1){
        lst <- sample(c(0, 1, 2, 3), size = 1)
        return(lst)
    }
    lst <- sample(c(0, 1, 2, 3), size = 2, replace = TRUE)
    i = 2
    while(n > 2){
        i = i + 1
        x <- sample(c(0, 1, 2, 3), size = 1, prob = mat[lst[i - 2] + 1, lst[i - 1] + 1, ])
        lst <- c(lst, x)
        n = n - 1
    }
    return(lst)
}

getStartSites <- function(n, totalLength, windowSize){
    lst <- sample(1:(totalLength - windowSize + 1), size = n, replace = TRUE)
    return(lst)
}

getNoiseLeft <- function(mat, startSites){
    lst <- list()
    for(i in startSites){
        if(i == 1){
            lst <- append(lst, list(c()))
        }
        else{
##            lst <- append(lst, list(getMarkovSequence(i - 1, mat)))
            tmp <- sample(c(0, 1, 2, 3), size = (i - 1), replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
            lst <- append(lst, list(tmp))
        }
    }
    print(c(length(lst[[1]]), startSites[1]))
    return(lst)
}

getNoiseRight <- function(n, mat, startSites, windowSize, totalLength){
    lst <- list()
    for(i in startSites){
        if((i + windowSize - 1) >= totalLength)
            lst <- append(lst, list(c()))
        else{
##            lst <- append(lst, list(getMarkovSequence(totalLength - i - windowSize + 1, mat)))
            tmp <- sample(c(0, 1, 2, 3), size = (totalLength - i - windowSize + 1), replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
            lst <- append(lst, list(tmp))
        }
    }
    return(lst)
}


getAll <- function(n, k, kcount, ksize, d, pos, dval, mats, totalLength){

    dst <- Map(function(x) getdst(d, dval), ksize)
    data <- Map(function(x) sampling(dst[[x]], ksize[x]), 1:kcount)
    startSites <- pos
    noiseLeft <- Map(function(x) getNoiseLeft(mats[[1]], startSites[[x]]), 1:kcount)
##    print(noiseLeft)
    noiseRight <- Map(function(x) getNoiseRight(ksize[x], mats[[1]], startSites[[x]], d, totalLength), 1:kcount)
    minStart <- min(unlist(startSites))
    maxStart <- max(unlist(startSites))
    m1 <- matrix(0, n, d)
    i <- 1
    j <- 1
    for(i1 in 1:kcount){
        for(j1 in 1:ksize[i1]){
            for(k1 in 1:d){
                m1[i, j] = data[[i1]][[j1]][k1]
                if(j == d){
                    j = 1
                    i = i + 1
                }
                else j = j + 1
            }
        }
    }

    {
    if(startSites[[1]][1] == maxStart){
        lftN <- c()
    }
    else {
        lftN <- rep(4, maxStart - startSites[[1]][1])
    }
    }
    {
    if(startSites[[1]][1] == minStart){
        rghtN <- c()
    }
    else{
        rghtN <- rep(4, startSites[[1]][1] - minStart)
    }
    }
    tm2 <- c(lftN, noiseLeft[[1]][[1]])
    {
    if(startSites[[1]][1] > (totalLength - d + 1)){
        tm2 <- c(tm2, head(m1[1, ], totalLength - startSites[[1]][1] + 1))
    }
    else{
        tm2 <- c(tm2, m1[1, ])
    }
    }
    tm2 <- c(tm2, noiseRight[[1]][[1]], rghtN)
    size <- length(tm2)
##    m2 <- matrix(0, n, 2*(maxStart - minStart))
    m2 <- matrix(0, n, size)
    
    i = 1 
    for(i1 in 1:kcount){
        for(j1 in 1:ksize[i1]){
            if(startSites[[i1]][j1] == maxStart)
                lftN <- c()
            else
                lftN <- rep(4, maxStart - startSites[[i1]][j1])
            if(startSites[[i1]][j1] == minStart)
                rghtN <- c()
            else
                rghtN <- rep(4, startSites[[i1]][j1] - minStart)
            tm2 <- c(lftN, noiseLeft[[i1]][[j1]])
            tm1 <- c(noiseLeft[[i1]][[j1]])
            if(startSites[[i1]][j1] > (totalLength - d + 1)){
                tm2 <- c(tm2, head(m1[i, ], totalLength - startSites[[i1]][j1] + 1))
                tm1 <- c(tm1, head(m1[i, ], totalLength - startSites[[i1]][j1] + 1))
            }
            else{
                tm2 <- c(tm2, m1[i, ])
                tm1 <- c(tm1, m1[i, ])
            }
            tm2 <- c(tm2, noiseRight[[i1]][[j1]], rghtN)
            tm1 <- c(tm1, noiseRight[[i1]][[j1]])
            m2[i, ] <- tm2
            data[[i1]][[j1]] <- tm1
            i = i + 1
        }
    }
    
######################################################################


##    print(m2[1, ])
    png(paste("original", "png", sep = "."), width = 400, height = 1000)
    image(t(m2), col = c('#008000', '#0000ff', '#ffa500', '#ff0000', '#ffffff'))
#    pdf(paste("originalCropped", "pdf", sep = "."), width = 2, height = 10)
#    image(t(m1), col = c('#008000', '#0000ff', '#ffa500', '#ff0000'))
    a1 <- c(ksize[1])
    for(i in 2:kcount)
        a1 <- c(a1, a1[i - 1] + ksize[i])
    a2 <- unlist(Map(function(x) x/n, a1))
    abline(h=a2, col="black", lwd = 4)
    ## abline(v=maxStart/(size + d), col="black", lwd = 2)
    ## abline(v=(maxStart + d - 1)/(size + d), col="black", lwd = 2)
    dev.off()
    
    l <- shuffle(data, startSites, ksize, kcount, n)
#    print(l)
 
    write(unlist(lapply(l[[1]], paste, collapse = " ")), paste("labels", "txt", sep = "."))

################################## Save labels and start sites #################################

    write(paste(l[[1]][1], "\t", startSites[[l[[1]][1]]][1]), "info.txt", append = FALSE)
    startSites[[l[[1]][1]]] <- tail(startSites[[l[[1]][1]]], -1)
    for(i in tail(l[[1]], -1)){
        write(paste(i, "\t", startSites[[i]][1]), "info.txt", append = TRUE)
        startSites[[i]] <- tail(startSites[[i]], -1)
    }
################################################################################################
    
    datafile <- paste("data", "txt", sep = ".")
    if(file.exists(datafile)) file.remove(datafile)
    m = matrix(0, n, d)
    
    
    for(i in 1:n)
        for(j in 1:d)
            m[i, j] = l[[2]][[i]][j]
    
    ## png(paste("jumbled", "png", sep = "."), col = c("red", "blue", "green", "yellow"))
    ## image(t(m))
    ## dev.off()

    for(i in 1:n){
        lft <- noiseLeft[[l[[1]][i]]][1]
        noiseLeft[[l[[1]][i]]] <- tail(noiseLeft[[l[[1]][i]]], -1)
        rght <- noiseRight[[l[[1]][i]]][1]
        noiseRight[[l[[1]][i]]] <- tail(noiseRight[[l[[1]][i]]], -1)


	############################# Apply reverse #####################
	strandInfo <- sample(c(0, 1), size = 1)
	if(strandInfo == 0){
        	      s <- paste(unlist(Map(function(x) vals[x + 1], l[[2]][[i]])), collapse = "")
	}
	else{
        	      s <- paste(unlist(Map(function(x) vals[3 - x + 1], l[[2]][[i]])), collapse = "")
	}

	#################################################################
        write(paste(">seq", i, sep = ""), "data.txt", append = TRUE)
        write(s, "data.txt", append = TRUE)
    }
#####################################################################################################
}

source("config.R")
mats <- generateMCmatrix(1)
## print(mats)
getAll(n, k, kcount, ksize, d, pos, dval, mats, l)
