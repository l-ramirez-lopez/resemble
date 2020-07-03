## build logo for resemble

require(resemble)
require(prospectr)

set.seed(200904)
x_orig <- x <- cbind(x = rnorm(5000), y = rnorm(5000))

plot(x)


diss_p <- f_diss(x, t(c(0,0)))

x_2 <- x[diss_p > 0.6 & diss_p < 9,]
x_3 <- x[diss_p > 1 & diss_p < 1.1,]
x <- x[diss_p < 0.4,]

set.seed(820206)
fsample <- sample(1:nrow(x), 30)
x <- x[fsample,]

set.seed(801124)
fsample_2 <- sample(1:nrow(x_2), 30)
x_2 <- x_2[fsample_2,]

set.seed(200904)
fsample_3 <- sample(1:nrow(x_3), 4)
x_3 <- x_3[fsample_3,]

x <- rbind(x, x_2, x_3)
plot(x)


diss_p <- f_diss(x, t(c(0,0)))
x <- x[order(diss_p), ]
diss_p <- sort(diss_p)

csize <- 5
psize <- 50
par(bg = "black",col.axis = NA, col.lab = NA)
plot(x, pch = 16, cex = 1.2* psize, col = NA, axes = FALSE, ann = FALSE)
points(x_orig, pch = 16, cex = csize, col = "black")
points(x_orig, pch = 16, cex = csize, col = rgb(0.8, 0.7, 0, 0.3))


points(x, pch = 16, cex = 1.2* psize, col = "black")

for( i in 1:nrow(x)) {
  lines(c(x[i,1], 0), c(x[i,2], 0), lwd = psize * 0.4, col = rgb(0.5, 0.5, 0.5, 0.6))
}

points(x, 
       pch = 16, 
       cex = psize, 
       col = heat.colors(1.7*length(diss_p), alpha = 0.6)[1:length(diss_p)])
points(c(0,0), c(0,0), pch = 16, cex = psize, col = "black")
points(c(0,0), c(0,0), pch = 16, cex = psize*0.8, col = "white")




