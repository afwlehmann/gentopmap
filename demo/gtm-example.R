#
# gtm-example.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#
# This file is part of gentopmap.
#
# gentopmap is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gentopmap is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gentopmap.  If not, see <http://www.gnu.org/licenses/>.


library(gentopmap)

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

