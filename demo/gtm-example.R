library(mygtm)

# Load the `oilflow` dataset.
data(oilflow)

grid <- {
  tmp <- seq(0, 1, length.out=5)
  as.matrix(expand.grid(tmp, tmp))
}

T <- scale(as.matrix(subset(oilflow, select=-c(label))))
stopifnot(ncol(T) == 12)
model <- computeGTM(T,
                    grid = c(100,100),
                    M = 16,
                    sigma = 1/4,
                    maxIter = 15,
                    verb=T)

plot.new()
par(mfrow=c(2,1),mar=c(2,2,2,0.2),xaxt="s",yaxt="s")
# First plot
P <- gtmPosteriorMean(model)
plot(P, pch='+', col=oilflow$label, main="Posterior Mean")
# Second plot
P <- gtmPosteriorMode(model)
plot(P, pch='+', col=oilflow$label, main="Posterior Mode")
