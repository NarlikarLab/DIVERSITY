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
                t <- rdirichlet(1, rep(1, nuc))
                for(k in 1:nuc)
                    mat[i, j, k] <- t[k]
            }
        }
        lst <- append(lst, list(mat))
    }
    return(lst)
}

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

getNoiseLeft <- function(n, mat, startSites){
    lst <- list()
    for(i in startSites){
        if(i == 1){
            lst <- append(lst, list(c()))
        }
        else
            lst <- append(lst, list(getMarkovSequence(i - 1, mat)))
    }
    return(lst)
}

getNoiseRight <- function(n, mat, startSites, windowSize, totalLength){
    lst <- list()
    for(i in startSites){
        if((i + windowSize - 1) == totalLength)
            lst <- append(lst, list(c()))
        else
            lst <- append(lst, list(getMarkovSequence(totalLength - i - windowSize + 1, mat)))
    }
    return(lst)
}

makeMatrix <- function(s){
    m <- unlist(Filter(function(x) x != "", unlist(strsplit(s, "\n"))))
    l <- length(unlist(strsplit(m[1], " ")))
    mat <- matrix(0, length(m), l)
    for(i in 1:length(m))
        mat[i, ] <- unlist(Map(function(x) as.numeric(x), unlist(strsplit(m[i], " "))))
    return(mat)
}

makeDataFrame <- function(l){
    n <- NULL
    pattern <- "([a-z]|[A-Z])+"
    lst <- list()
    for(i in l){
        while(i[1] == "\n")
            i = tail(i, -1)
        reg <- regexpr(pattern, i)
        n <- c(n, regmatches(i, reg))
        lst <- append(lst, list(makeMatrix(regmatches(i, reg, invert = TRUE)[[1]][2])))
    }
    names(lst) <- n
    return(lst)
}

motifPos <- function(motifMat, motifs, lSeq){
    startPos <- NULL
    start <- 1
    end <- lSeq
    lMotifs <- unlist(Map(function(x) length(motifMat[[x]][ , 1]), motifs))
    lm <- length(lMotifs)
    while(lm > 0){
        samp <- sample(start:(end - sum(lMotifs) + 1), size = 1)
        start <- samp + as.vector(lMotifs[1])
        startPos <- c(startPos, samp)
        lMotifs <- tail(lMotifs, -1)
        lm <- lm - 1
    }
    return(startPos)
}


generateData <- function(n, k, d, startPos, motifMat, motifList){
    sizeK <- length(k)
    sizeG <- unlist(Map(function(x) round(n*k[x]), 1:(sizeK-1)))
    sizeG <- c(sizeG, (n-sum(sizeG)))
    undst <- Map(function(x) rep(0.25, 4), 1:d)
    for(i in 1:d)
        undst[[i]] <- unlist(Map(function(x) x, undst[[i]]))
    undst <- Map(function(x) undst, 1:sizeK)
    for(i in 1:sizeK){
        for(j in 1:length(startPos[[i]])){
            for(k1 in 1:(length(motifMat[[motifList[[i]][j]]][, 1]))){
                undst[[i]][startPos[[i]][j] + k1 - 1] <- list(motifMat[[motifList[[i]][j]]][k1, ])
            }
        }
    }
    data <- Map(function(x) sampling(undst[[x]], sizeG[x]), 1:sizeK)
    return(list(data, sizeK, sizeG))
}

readMotifs <- function(filename){
    txt <- paste(readLines(filename), collapse = "\n")
    pattern <- "(\n)*([a-z]|[A-Z])+\n((([0-9]+.[0-9]+)| )+|0|1|\n)+"
    leftOut <- FALSE
    lst <- list()
    while(txt != ""){
        reg <- regexpr(pattern, txt)
        lst <- append(lst, list(regmatches(txt, reg)))
        txt <- regmatches(txt, reg, invert = TRUE)[[1]][2]
        leftOut <- (txt == "")
    }
    motifs <- makeDataFrame(lst)
}

shuffle <- function(data, sizeG, sizeK, n){
  groupList <- unlist(Map(function(x) rep(x, times = sizeG[x]), 1:sizeK))
  groupList <- sample(groupList, size = n, replace = FALSE, prob = NULL)
  groups <- Map(function(x) sample(data[[x]], size = sizeG[x], replace = FALSE, prob = NULL), 1:sizeK)

  lst <- list()
  i = 1
  while(i<=n){
    lst <- append(lst, head(groups[[groupList[i]]], 1))
    groups[[groupList[i]]] <- tail(groups[[groupList[i]]], -1)
    i <- i + 1
  }
  return(list(groupList, lst))
}

getDataMatrix <- function(n, d, data){
    i <- 1
    j <- 1
    m1 <- matrix(0, n, d)
    for(i1 in 1:data[[2]]){
        for(j1 in 1:data[[3]][i1]){
            for(k1 in 1:d){
                m1[i, j] = data[[1]][[i1]][[j1]][k1]
                if(j == d){
                    j = 1
                    i = i + 1
                }
                else j = j + 1
            }
        }
    }
    return(m1)
}

shuffleAndSave <- function(data, n, startSites, noiseLeft, noiseRight){
    l <- shuffle(data[[1]], data[[3]], data[[2]], n)
    
#    write(unlist(lapply(l[[1]], paste, collapse = " ")), paste("labels", "txt", sep = "."))
    write(paste(l[[1]][1], "\t", startSites[[l[[1]][1]]][1]), "info.txt", append = FALSE)
    startSites[[l[[1]][1]]] <- tail(startSites[[l[[1]][1]]], -1)
    for(i in tail(l[[1]], -1)){
        write(paste(i, "\t", startSites[[i]][1]), "info.txt", append = TRUE)
        startSites[[i]] <- tail(startSites[[i]], -1)
    }
    
    datafile <- paste("data", "txt", sep = ".")
    if(file.exists(datafile)) file.remove(datafile)
    for(i in 1:n){
        lft <- noiseLeft[[l[[1]][i]]][1]
        noiseLeft[[l[[1]][i]]] <- tail(noiseLeft[[l[[1]][i]]], -1)
        rght <- noiseRight[[l[[1]][i]]][1]
        noiseRight[[l[[1]][i]]] <- tail(noiseRight[[l[[1]][i]]], -1)
        s <- paste(unlist(Map(function(x) vals[x + 1], c(lft, l[[2]][[i]], rght))), collapse = "")
        write(paste(">seq", i, sep = ""), "data.txt", append = TRUE)
        write(s, "data.txt", append = TRUE)
    }
}

getArgs <- function(){
    args<-commandArgs(trailingOnly = T)
    if(length(args) < 5){
        print("Usage: Rscript generate.R <motif file> <number of sequences> <sequence length> <seq1, seq2> <seq1, seq1> ...")
        print("Example: Rscript generate.R motifs.txt 1000 100 50 A,B C,C,D D,B A")
        stop()
    }
    lst <- list(args[1], as.numeric(args[2]), as.numeric(args[3]), as.numeric(args[4]))
    m <- list()
    for(i in tail(args, -4))
        m <- append(m, strsplit(i, ","))
    return(append(lst, list(m)))
}

getFullMatrix <- function(minStart, maxStart, data, startSites, noiseLeft, noiseRight, m1){
    m2 <- matrix(0, n, 2*(maxStart - minStart) + d)
    i = 1 
    for(i1 in 1:data[[2]]){
        for(j1 in 1:data[[3]][i1]){
            if(startSites[[i1]][j1] == maxStart)
                lftN <- c()
            else
                lftN <- rep(4, maxStart - startSites[[i1]][j1])
            if(startSites[[i1]][j1] == minStart)
                rghtN <- c()
            else
                rghtN <- rep(4, startSites[[i1]][j1] - minStart)
            m2[i, ] = c(lftN, noiseLeft[[i1]][[j1]], m1[i, ], noiseRight[[i1]][[j1]], rghtN)
            i = i + 1
        }
    }
    return(m2)
}


set.seed(5)
lst <- getArgs()
filename <- lst[[1]]
n <- lst[[2]]
l <- lst[[3]]
d <- lst[[4]]
kcount <- length(lst[[5]])
k <- rdirichlet(1, rep(2, kcount))
motifList <- lst[[5]]

motifs <- readMotifs(filename)
mats <- generateMCmatrix(kcount)
sp <- Map(function(x) as.vector(motifPos(motifs, x, d)), motifList)
data <- generateData(n, k, d, sp, motifs, motifList)
startSites <- Map(function(x) getStartSites(data[[3]][x], l, d), 1:data[[2]])
noiseLeft <- Map(function(x) getNoiseLeft(data[[3]][x], mats[[x]], startSites[[x]]), 1:data[[2]])
noiseRight <- Map(function(x) getNoiseRight(data[[3]][x], mats[[x]], startSites[[x]], d, l), 1:data[[2]])
minStart <- min(unlist(startSites))
maxStart <- max(unlist(startSites))

m1 <- getDataMatrix(n, d, data)
m2 <- getFullMatrix(minStart, maxStart, data, startSites, noiseLeft, noiseRight, m1)

png(paste("imageWindow", "png", sep = "."), width = 600, height = 1000)
image(t(m1), col = c('#008000', '#0000ff', '#ffa500', '#ff0000'))
a1 <- c(data[[3]][1])
for(i in 2:data[[2]])
    a1 <- c(a1, a1[i-1] + data[[3]][i])
a2 <- unlist(Map(function(x) x/n, a1))
abline(h=a2, col="black", lwd = 2)
dev.off()

png(paste("imageFull", "png", sep = "."), width = 600, height = 1000)
image(t(m2), col = c('#008000', '#0000ff', '#ffa500', '#ff0000', '#ffffff'))
abline(h=a2, col="black", lwd = 2)
abline(v=maxStart/(2*(maxStart - minStart) + d), col="black", lwd = 2)
abline(v=(maxStart + d - 1)/(2*(maxStart - minStart) + d), col="black", lwd = 2)
dev.off()


shuffleAndSave(data, n, startSites, noiseLeft, noiseRight)
#print(tmp)
