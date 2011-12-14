# Computes the GTM from 3-dimensional data-space to 2-dimensional latent-space.
# Data-space points are drawn from the Gaussians. This example is silly because
# of the low dimensionality of the data-space.

library(mygtm)
library(mytools) # for computeNDimensionalGrid

T <- rbind(rmvnorm(n=100, mu=c(0,0,0), sigma=diag(c(1,1,1))),
           rmvnorm(n=100, mu=c(-10,-10,-10), sigma=diag(c(1,3,2))),
           rmvnorm(n=100, mu=c(7,23,0), sigma=diag(c(4,7,1))))
grid <- computeNDimensionalGrid(2, 10^2)
model <- gtm.compute(T, grid, sigma=1/40, K=nrow(grid)*100, epsilon=1e-5, 
                     maxIterations=10,
                     callback=function(iter,llh) {
                       cat(sprintf("Cycle %d (%f)\n", iter, llh))
                     })

P <- gtm.project(model$X, model$Rin)
extX = range(P[,1])
extY = range(P[,2])

plot.new()
par(mfrow=c(2,1),mar=c(2,2,0,0),xaxt="s",yaxt="s")
plot(model$X[,1], model$X[,2], pch=".")
grid()
plot(P[1:100,1], P[1:100,2], pch="+", col="red", xlim=extX, ylim=extY)
par(new=TRUE, xaxt="n", yaxt="n")
plot(P[101:200,1], P[101:200,2], pch="o", col="green", xlim=extX, ylim=extY)
par(new=TRUE, xaxt="n", yaxt="n")
plot(P[201:300,1], P[201:300,2], pch="*", col="blue", xlim=extX, ylim=extY)
grid()
