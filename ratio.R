hm <- rnorm(35, 10000000, 5000000)
ht1a <- rnorm(1500, 5000000, 1000000)
ht1 <- ht1a[which(ht1a < 7.5e+06)] #non-recombinant region = 7.5m-12.5m
ht2a <- rnorm(1500, 15000000, 1000000)
ht2 <- ht2a[which(ht2a > 1.25e+07)] #non-recombinant region = 7.5m-12.5m
ht <- c(ht1, ht2)

hmd <- density(hm, from=0, to=18585056)
htd <- density(ht, from=0, to=18585056)
x <- (1:512)*36298
y <- hmd$y/htd$y
plot(x, y)

# I think this will do