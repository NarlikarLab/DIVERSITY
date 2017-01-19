library(MCMCpack)

source("config.R")

getdst <- function(n){
    y <- Map(function(x) rdirichlet(1, rep(runif(1, 0, 0.1), 4)), 1:n)
    return(y)
}

revdst <- function(dst){
       for(i in 1:length(dst)){
       	     tmp <- dst[[i]][1]
	     dst[[i]][1] <- dst[[i]][4]
	     dst[[i]][4] <- tmp
	     tmp <- dst[[i]][2]
	     dst[[i]][2] <- dst[[i]][3]
	     dst[[i]][3] <- tmp
       }
       return(rev(dst))
}

##pos <- Map(function(x) sample(1:(l - d + 1), size = x, replace = TRUE), ksize)
pos <- Map(function(x) sample((d + 1):(l - 3*d + 1), size = x, replace = TRUE), ksize)


cSize <- ksize
labels <- unlist(Map(function(x) rep(x, times = cSize[x]), 1:kcount))
labels <- sample(labels, size <- n, replace = FALSE, prob = NULL)

distribution <- Map(function(x) getdst(d), 1:kcount)
#print(distribution)

data <- list(c(), c(), c())
j <- 1
for(i in labels){
    start <- head(pos[[i]], 1)
    pos[[i]] <- tail(pos[[i]], -1)
#    if(i == 3) print(pos[[i]])
    # print("Distribution")
    # print(distribution[[i]])
    # print("Reversed")
    # print(revdst(distribution[[i]]))
##################################### Apply rev ###################################################
    strandInfo <- sample(c(0, 1), size = 1)
    if(strandInfo == 0){
    m <- paste(unlist(Map(function(x) sample(c("A", "C", "G", "T"), size = 1, prob = x), distribution[[i]])), collapse = "")
    }
    else{
    m <- paste(unlist(Map(function(x) sample(c("A", "C", "G", "T"), size = 1, prob = x), revdst(distribution[[i]]))), collapse = "")
    }

####################################################################################################
    left <- paste(sample(c("A", "C", "G", "T"), size = start - 1, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25)), collapse = "")
    right <- paste(sample(c("A", "C", "G", "T"), size = l - start - d + 1, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25)), collapse = "")
    seq <- paste(paste(">seq_", j, "\n", sep = ""), paste(left, m, right, sep = ""), sep = "")
    data[[i]] <- c(data[[i]], m)
##    print(seq)
    
    if(j == 1){
        write(seq, "new_data.txt", append = FALSE)
        write(i, "new_labels.txt", append = FALSE)
        write(start - 1, "new_start.txt", append = FALSE)
    }
    else{
        write(seq, "new_data.txt", append = TRUE)
        write(i, "new_labels.txt", append = TRUE)
        write(start - 1, "new_start.txt", append = TRUE)
    }
    j <- j + 1
    ## print(paste(m, collapse = ""))
    ## print(left)
}

write(d, "new_width.txt", append = FALSE)
for(i in 2:kcount) write(d, "new_width.txt", append = TRUE)

print(data[[1]])
print(data[[2]])
print(data[[3]])
