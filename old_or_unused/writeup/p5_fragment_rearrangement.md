Fragment rearrangement part 5 
========================================================


```r
hm <- rnorm(70, 1e+07, 2e+06)
ht1 <- rnorm(35, 7e+06, 2e+06)
ht2 <- rnorm(35, 1.3e+07, 2e+06)
ht <- c(ht1, ht2)

density(hm)$y[1]/density(ht)$y[1]  #this gives good values to use
```

```
## [1] 4.188
```

```r

hm <- rnorm(70, 1e+07, 3e+06)
ht1 <- rnorm(350, 5e+06, 1e+06)
ht2 <- rnorm(350, 1.5e+07, 1e+06)
ht <- c(ht1, ht2)
plot(1:512, density(hm)$y/density(ht)$y)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 


Got these distributions, now you got to work out a way to get signal ratio on y axis, with position (bp) on x axis, for kernel density plot.
Actually no, we want signal ratio on the y axis.
