library("rjson")
require(graphics)
library(LPE)


means <- fromJSON(file = "/home/pf/hoffman2/twindata/CGmap/lowess_fragment_means.json")
new = unlist(means["new"])
old = unlist(means["old"])

#dat <- cbind(unlist(means["old"]), unlist(means["new"]))

#load(file = 'dat.RData')

#head(dat)
#colnames(dat) <- c('old', 'new')
#dat <- as.data.frame(dat)

new.loess <- loess(new ~ old, data = data.frame(x = old, y = new), span = 0.001)

str(new.loess)
plot(new.loess, pch = '.', xlab = 'batch 1', ylab = 'batch 2')
abline(0, 1, col = 'blue', lwd = 2)

new.pred <- predict(new.loess, data.frame(x = old))

lines(old, new.pred, col = 'red', lwd = 2)

str(new.pred)
str(new.loess)

dat1 <- with(new.loess[c('x', 'y')], {
  new_y <- y + (x - new.pred)
  new_x <- x
  return(cbind(new_x, new_y))
})

plot(dat1, pch = '.', xlab = 'batch1', ylab = 'normalized batch 2')
abline(0,1, col = 'blue', lwd = 3)

colnames(dat1) <- c('old', 'normalized_new')
rownames(dat1) <- NULL
sink("/home/pf/hoffman2/twindata/CGmap/loess_output_for_fragment_means.json")
cat(toJSON(list((dat1[,1]), (dat1[,2]))))
#cat(toJSON((dat1[,2])))

cat("\n")
sink()

