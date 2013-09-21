library(mygtm)

# Load the `oilflow` dataset.
data(oilflow)
T <- scale(as.matrix(subset(oilflow, select=-c(label))))
stopifnot(ncol(T) == 12)

model <- computeGTM(T,
                    grid = c(50,50),
                    M = 16,
                    sigma = 2/50*4,
                    maxIter = 15,
                    verb=TRUE)

plot.new()
par(mfrow=c(3,1),mar=c(2,2,2,0.2),xaxt="s",yaxt="s")
# First plot
P <- gtmPosteriorMean(model)
plot(P, pch='+', col=oilflow$label, main="Posterior Mean")
# Second plot
P <- gtmPosteriorMode(model)
plot(P, pch='+', col=oilflow$label, main="Posterior Mode")
# Third plot
aux <- svd(T)
P <- T %*% aux$v[,1:2]
plot(P, pch='+', col=oilflow$label, main="PCA")
