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
                    grid = grid,
                    sigma = 1/10,
                    K = nrow(grid)*100,
                    maxIter = 15,
                    callback=function(iter,llh) {
                      cat(sprintf("Cycle %d (%f)\n", iter, llh))
                    })

plot.new()
par(mfrow=c(2,1),mar=c(2,2,0,0)+0.1,xaxt="s",yaxt="s")
# First plot
plot(model$X, pch=".")
grid()
# Second plot
P <- gtmProject(model$X, model$Rin)
plot(P, pch='+', col=oilflow$label)
