# Functions in R
Emily Chase (PID: A14656894)

- [Write a function that adds
  numbers](#write-a-function-that-adds-numbers)
- [Another function](#another-function)
  - [Conditionals interlude](#conditionals-interlude)
  - [Now integrate it into functions](#now-integrate-it-into-functions)
- [Design a protein generating
  function](#design-a-protein-generating-function)
  - [Generate sequences length between 6 and
    12](#generate-sequences-length-between-6-and-12)
  - [Search in BLAST](#search-in-blast)

## Write a function that adds numbers

``` r
add <- function(x, y=1){
  x+y
  
}
```

Call the function

``` r
add(10,100)
```

    [1] 110

## Another function

Write a function to generate random nucleotide sequeces of a user
specified length:

The `sample` function can be helpful here. It samples randomly from an
input vector.

``` r
sample(c("A", "C", "G", "T"), size=5, replace=TRUE) # can only have size>n if replace=TRUE
```

    [1] "C" "T" "T" "A" "C"

`sample` gives us a vector, but what if we want a fasta format instead?
ie a 1 element long character vector

``` r
s <- sample(c("A", "C", "G", "T"), size=5, replace=TRUE) # can only have size>n if replace=TRUE
paste(s, collapse="")
```

    [1] "CACGA"

We can use `paste` (kinda like `join` in python). `paste` has two
arguments that are new to me though:

1.  `collapse` = \_\_\_ ~ determines what will separate **ELEMENTS** in
    a vector

2.  `sep` = \_\_\_ ~ determines what will separate **arguments**. you
    can give `paste` multiple vectors and it will combine them
    pair-wise, and you can tell it what to separate those pairs by.

### Conditionals interlude

``` r
fasta <- TRUE
if (fasta){
  cat("HELLO WORLD!")
} else {
  cat("NOPE")
}
```

    HELLO WORLD!

``` r
# what's cat()

print("HI")
```

    [1] "HI"

``` r
print(c("HI", "BYE"))
```

    [1] "HI"  "BYE"

``` r
cat(c("HI", "BYE"))
```

    HI BYE

### Now integrate it into functions

Add the ability to return a multi element vector or a single element
fasta like vector

``` r
generate_fasta <- function(size=50, fasta=TRUE){
  v <- sample(c("A", "G", "C","T"), size=size, replace=TRUE)
  
  
  if (fasta){
    ans <- paste(v, collapse="")
    return(ans)
  } else {
    return(v)
  }
  
}

generate_fasta(fasta=TRUE)
```

    [1] "TCGAGATGTCTCAATCGTAGAAACGAAGGCTAACTCGCCCACCCATCAGA"

``` r
generate_fasta(fasta=FALSE)
```

     [1] "A" "G" "T" "A" "C" "A" "A" "T" "A" "T" "A" "C" "A" "G" "T" "A" "A" "A" "T"
    [20] "T" "T" "C" "C" "T" "A" "A" "T" "T" "C" "C" "A" "A" "A" "A" "T" "A" "A" "G"
    [39] "A" "C" "G" "C" "G" "G" "C" "C" "C" "A" "A" "C"

# Design a protein generating function

``` r
generate_pfasta <- function(size=50, fasta=TRUE){
  v <- sample(c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), size=size, replace=TRUE)
  
  
  if (fasta){
    ans <- paste(v, collapse="")
    return(ans)
  } else {
    return(v)
  }
  
}

generate_pfasta(fasta=TRUE)
```

    [1] "FAPCACNIPNYYHCSFILIWIEGMPEILTITCEYFHIPTTLDNEPFLCGN"

``` r
generate_pfasta(fasta=FALSE)
```

     [1] "I" "E" "F" "C" "Y" "P" "H" "R" "G" "F" "V" "V" "I" "P" "C" "L" "D" "L" "R"
    [20] "K" "Q" "Q" "G" "D" "F" "V" "R" "D" "A" "L" "Q" "N" "H" "W" "D" "F" "M" "G"
    [39] "Q" "Q" "S" "K" "D" "N" "V" "Y" "S" "A" "Q" "I"

``` r
generate_pfasta(6)
```

    [1] "WSWMRE"

## Generate sequences length between 6 and 12

Use our new `generate_pfasta()` function to make random protein
sequences of length 6 to 12 (ie one length 6, one length 7, up to len
12)

``` r
# brute force
generate_pfasta(6)
```

    [1] "VDFDYP"

``` r
generate_pfasta(7)
```

    [1] "HENWVHE"

``` r
generate_pfasta(8)
```

    [1] "HHLAWLMR"

``` r
# using a for loop

#less efficient
lengths <- 6:12
for(i in lengths){
  cat(">",i, "\n", sep="")
  aa <- generate_pfasta(i)
  cat(aa)
  cat("\n")
}
```

    >6
    RHFMSC
    >7
    QHNIFQD
    >8
    ADAEQELT
    >9
    SASYVKKAL
    >10
    QKTYCWIIHG
    >11
    SFARVPGNLLN
    >12
    RDDMFRWDLQFE

``` r
# more efficient

sapply(lengths,generate_pfasta)
```

    [1] "DGHQKS"       "CDGSETH"      "TCERIHSN"     "MMGISWANK"    "VRWANDYWYE"  
    [6] "HQIFSGIFSDV"  "WKEQHSWPDQYQ"

## Search in BLAST
