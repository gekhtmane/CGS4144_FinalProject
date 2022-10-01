# First make sure BioInfoFile1 is in your
# data environment before running this script
BioInfoFile1 <- read.csv("C:/Users/surfc/Desktop/BioInformatics/BioInfoFile1.csv", header=TRUE)
x <- BioInfoFile1[,-(1:9)]
#boxplot(log2(x+1))
y <- log2(x+1)
z <- apply(X = y, MARGIN = 1, FUN = median)
plot(density(z))
end()       
 