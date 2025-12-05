# Lab 4
Emily Chase (PID: A14656894)
2025-10-10

``` r
# our first script
print("hello world")
```

    [1] "hello world"

``` r
x <- 1:50

plot(x)
```

![](class_04_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
plot(x, col='pink')
```

![](class_04_files/figure-commonmark/unnamed-chunk-1-2.png)

``` r
plot(x, sin(x), col='pink', typ="l", lwd="5")
```

![](class_04_files/figure-commonmark/unnamed-chunk-1-3.png)

``` r
log(10, base=10)
```

    [1] 1

``` r
plot(x, log(x))
```

![](class_04_files/figure-commonmark/unnamed-chunk-1-4.png)

``` r
log10(10)
```

    [1] 1
