args <- commandArgs(trailingOnly = T)
dataFile <- args[1]
labelsFile <- args[2]
startFile <- args[3]
motifWidthFile <- args[4]

getData <- function(data, index, start, width){
    l <- length(data[[index]])
    if(start < l)
        return(data[[index]][start:(start + width - 1)])
    i = 1
    j = l
    lst <- unlist(Map(function(x) if((3-x) < 0) 4 else (3-x), rev(data[[index]])))
    return(lst[(start - l - 1):(start + width - l - 1 - 1)])
        
}


#motifL <- c(8, 8, 8)
mapV <- c(0, 1, 2, 3, 0, 1, 2, 3, 4)
names(mapV) <- c("a", "c", "g", "t", "A", "C", "G", "T", "N")
l1 <- readLines(dataFile)
data <- list()
for(i in 1:(length(l1)/2)){
    t <- as.vector(unlist(Map(function(x) mapV[[x]], unlist(strsplit(l1[[2*i]], "")))))
    data <- append(data, list(t))
}

labels <- as.vector(unlist(Map(function(x) as.numeric(x), readLines(labelsFile))))
labels <- unlist(Map(function(x) x + 1, labels))
startSites <- as.vector(unlist(Map(function(x) as.numeric(x) + 1, readLines(startFile))))
motifL <- as.vector(unlist(Map(function(x) as.numeric(x), readLines(motifWidthFile))))

#print(startSites)
n <- length(Filter(function(x) x > 0, startSites))
positions <- unlist(Filter(function(x) startSites[x] > 0, 1:length(labels)))
data <- Map(function(x) data[[x]], positions)
labels <- unlist(Map(function(x) labels[x], positions))
startSites <- unlist(Map(function(x) startSites[x], positions))
                 

## for(i in 1:n){
##     if(labels[i] == 3) print(startSites[i])
## }
motif <- Map(function(x) getData(data, x, startSites[x], motifL[labels[x]]), 1:n)
## motif <- Map(function(x) data[[x]][(startSites[x]):(startSites[x] + motifL - 1)], 1:n)
pos <- Map(function(x) Filter(function(y) labels[y] == x, 1:length(labels)), 1:max(labels))

motif <- Map(function(x) Map(function(y) motif[[y]], x), pos)

maxMotifL <- max(motifL)

m <- matrix(0, n, maxMotifL)
k <- 1

for(i in motif){
    for(j in i){
        lst <- c(j)
#        print(length(lst))
        if(length(lst) < maxMotifL) for(i1 in length(lst):(maxMotifL-1)) lst <- c(lst, 4)
#        print(length(lst))
#        print(maxMotifL)
        m[k, ] <- lst
        k <- k + 1
    }
}

#print(m)

sizes <- unlist(Map(function(x) length(Filter(function(y) labels[y] == x, 1:length(labels))), 1:max(labels)))

cols <- c('#008000', '#0000ff', '#ffa500', '#ff0000')
if(max(motifL) != min(motifL)) cols <- c(cols, '#ffffff')
png("learnedImage.png")
image(t(m), col = cols)
a1 <- c(sizes[1])
for(i in 2:max(labels))
    a1 <- c(a1, a1[i - 1] + sizes[i])
a2 <- unlist(Map(function(x) x/n, a1))
abline(h = a2, col = "black", lwd = 4)
dev.off()
