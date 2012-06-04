library("rjson")
require(graphics)
library(LPE)

means <- fromJSON(file = "CGmap/lowess_fragment_means.json")
#png('output/batch1_vs_2_fragment_means.png', width = 800, height = 800)


#means <- fromJSON(file = "CGmap/lowess_means.json")
#png('output/batch1_vs_2_site_means.png', width = 800, height = 800)

plot(unlist(means["new"]), unlist(means["old"]), main = "lowess(means)", pch = '.')

#lines(lowess(unlist(means["new"]), unlist(means["old"]), iter=10), col = 2)

fit = lowess(unlist(means["new"]), unlist(means["old"]), f=.0001, iter =10)

lines(fit, col = 'red', lwd = 3)

#plot(unlist(means["new"]), unlist(means["old"]) + fit$x - fit$y, main = "lowess(means)", pch = '.')


lines(c(0,1), c(0,1), col = 'blue', lwd = 3)

#legend(5, 120, c(paste("f = ", c("2/3", ".2"))), lty = 1, col = 2:3)

dev.off()


# PLOT NORMALIZED VALUES
#png('output/lowess_normalized.png', width = 800, height = 800)
#fit = lowess(unlist(means["new"]), unlist(means["old"]), f=.001, iter =10)
#plot(unlist(means["new"]), unlist(means["old"]) + fit$x - fit$y, main = "lowess(means)", pch = '.')
#lines(c(0,1), c(0,1), col = 'blue', lwd = 3)
#legend(5, 120, c(paste("f = ", c("2/3", ".2"))), lty = 1, col = 2:3)
#dev.off()




#sink("CGmap/lowess_output_for_fragment_means.json")
#cat(toJSON(fit))


