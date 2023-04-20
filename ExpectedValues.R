# Node balance function
myf <- function(p) {
  if(p==0 | p==1) return(0)
  return(-p*log2(p) - (1-p)*log2(1-p))
}

# Average node balance, which approaches 1 / (2 log(2)) ≈ myf(0.199709902553977)
sumf <- function(n) {
  j <- 1:(n-1)
  sum(sapply(j/n, myf)) / (n-1)
}

# exact expression for E(I_S) for uniform model, from Elizalde 2023:
f1a <- function(n) n * (4^(n-1) / choose(2*n-2, n-1) - 1)

# exact expression for E(I_S) for uniform model, from Mir et al. 2013:
f2 <- function(n) {
  prod1 <- 1
  for(i in c(1:(n-1))) prod1 <- prod1 * (2*n-2*i)
  prod2 <- 1
  for(i in c(1:(n-1))) prod2 <- prod2 * (2*n-2*i-1)
  return(n * (prod1/prod2 - 1))
}

# good approximation to E(I_S) for uniform model, from Stirling's approximation:
f3 <- function(n) n * sqrt((n-1)*pi) - n

# limit approximation to E(I_S) for uniform model:
f3a <- function(n) n * sqrt(n*pi)

# good approximation to n log2(n) / E(I_S) for uniform model:
f3_inverse <- function(n) log2(n) / (sqrt((n-1)*pi) - 1)

harmonic_fn <- function(n) sum(1/(2:n))

# exact expression for n log2(n) / E(I_S) for Yule process:
g1 <- function(n) log2(n) / (2*harmonic_fn(n))

# good approximation to n log2(n) / E(I_S) for Yule process:
g2 <- function(n, terms = 1) {
  den <- log(n) -digamma(1) - 1 + 1/(2*n)
  if(terms > 1) den <- den - 1/(8*n^2)
  res <- log2(n) / (2*den)
  return(res)
}

n_exact <- 2:10

# Exact values of n log2(n) / E(J^1) and E(I_S) for Yule process:
EJ1 <- c(2, 5, 216/25, 728/57, 1162800/67217, 199806750/9017743, 27.29901, 32.68993, 38.30246)
EIS <- c(2, 5, 26/3, 77/6, 87/5, 223/10, 27.48571, 32.92143, 38.57937)

EJ1_over_2n <- c(1/2, 5/6, 27/25, 364/285, 96900/67217, 28543821/18035486)
EIS_over_2n <- c(1/2, 5/6, 13/12, 77/60, 29/20, 223/140)

EJ1_over_2n_plus_1 <- c(3/2, 11/6, 52/25, 649/285, 164117/67217, 46579307/18035486)
EIS_over_2n_plus_1 <- c(3/2, 11/6, 25/12, 137/60, 49/20, 363/140)
ratios <- c(1, 1, 624/625, 2596/2603, 3282340/3293633, 465793070/467634387)

# Exact values of n log2(n) / E(J^1) and E(I_S) for uniform model:
EJ1_Unif <- c(2, 5, 360/41, 3822/289, 18.25643, 23.81979, 29.87282, 36.38201, 43.31989)
EIS_Unif <- c(2, 5, 44/5, 93/7, 18.38095, 24.0303, 30.19114, 36.82937, 43.91691)

EJ1_Unif_over_n <- c(1, 5/3, 90/41, 3822/1445, NA, NA)
EIS_Unif_over_n <- c(1, 5/3, 11/5, 93/35, NA, NA)

# E(I_S) for uniform model:
n <- 2:50
plot(sapply(n, f2) ~ n) # exact expression
lines(sapply(n, f3) ~ n, col = "blue", lwd = 2) # good approximation
lines(sapply(n, f3a) ~ n, col = "orange") # limit approximation
points(EIS_Unif ~ n_exact, col = "red", pch = 3) # Exact values

# n log2(n) / E(I_S):
#n <- unique(c(2:16, round(2^seq(4, 9, length = 40))))
n <- 2:20
pdf("ExpectedValuesCloseUp.pdf", width = 5, height = 5)
par(mar = c(4,4,1,1))
plot(n*log2(n) / sapply(n, f1a) ~ n, xlim = c(2, 20), ylim = c(0.6, 1),
     xlab = "number of leaves, n", ylab = "index value") # uniform model exact expression
lines(sapply(n, f3_inverse) ~ n, col = "grey", lwd = 2) # uniform model approximation
points(sapply(n, g1) ~ n) # Yule process exact expression
lines(sapply(n, g2) ~ n, col = "grey", lwd = 2) # Yule process approximation
points(n_exact*log2(n_exact) / EJ1 ~ n_exact, col = "red", pch = 3) # Yule process exact values
points(n_exact*log2(n_exact) / EJ1_Unif ~ n_exact, col = "red", pch = 3) # uniform model exact values
#abline(h = 1 / (2 * log(2)), lty = 2)
#points(sapply(n, sumf) ~ n, col = "orange") # balance score of root node
legend(x = 2, y = 0.4, legend = c("n log n / E(I_S)", 
                                "E(J^1)"), 
       col = c("black", "red"), 
       pch = c(1, 3), 
       bty = "n")
legend(x = 2, y = 0.2, legend = c("n log n / [approximation of E(I_S)]"), 
       col = c("grey"), 
       lty = 1, 
       bty = "n")
dev.off()

# E(1 / X) >= 1 / E(X)
# E(J1) = E(log2(n) / I_S) = log2(n) / EJ1 >= log2(n) / E(I_S) = log2(n) / EIS
# EJ1 <= EIS

J1_caterpillar <- function(n) 2*n*log2(n) / ((n-1)*(n+2))

# I_S for sawtooth tree:
I_S_sawtooth <- function(n) {
  if((n - 1) %% 3) stop("n - 1 must be divisible by 3")
  add <- 0
  if(n > 4) {
    maxi <- (n-4)/3
    add <- 6 * sum(1:maxi)
  }
  return(8*(n-1)/3 + add)
}

J1_sawtooth <- function(n) n*log2(n) / I_S_sawtooth(n)

#######

library(Zseq)

# n!!
acdoub <- function(n) {
  if(n %% 2 == 0) {
    m <- n/2 + 1
    return(Factorial.Double(m, odd = FALSE)[m])
  }
  m <- round((n+3)/2)
  return(Factorial.Double(m, odd = TRUE)[m])
}

# OEIS sequence A305578:
A305578 <- function(n, k) {
  s <- 0
  for(k in 0:n) s <- s + choose(n, k) * acdoub(k) * acdoub(n - k)
  return(s)
}

#######

neg_binom <- function(r, p, k) choose(k + r - 1, k)*(1-p)^k*p^r

k <- 0:20
plot(sapply(k, neg_binom, p = 1/2, r = 2) ~ k, col = "#FF0000")
lines(sapply(k, neg_binom, p = 1/2, r = 2) ~ k, col = "#FF0000")
points(sapply(k-1, neg_binom, p = 1/2, r = 3) ~ k, col = "#FF0099")
lines(sapply(k-1, neg_binom, p = 1/2, r = 3) ~ k, col = "#FF0099")

#######

I_S_index <- function(tree) { # Sackin's index
  I_S <- 0 # initialise sum, which will eventually be equal to Sackin's index
  for(i in 1:length(tree)) { # loop over all levels of the tree
    if(i < length(tree)) { # if the current level isn't the bottom level
      leaves_at_level_i <- sum(tree[[i]]) - sum(tree[[i + 1]]) # number of leaves on current level
    } else {
      leaves_at_level_i <- sum(tree[[i]]) # number of leaves on bottom level
    }
    I_S <- I_S + leaves_at_level_i * i
  }
  return(I_S)
}

H_index <- function(tree) {
  sum1 <- 0 # initialise sum, which will eventually be equal to H * I_S
  I_S <- I_S_index(tree) # Sackin's index
  for(i in 1:length(tree)) { # loop over all levels of the tree
    g <- tree[[i]] # vector of all branch weights at current level
    S <- sum(g) # sum of branch weights at current level
    sum1 <- sum1 - sum(g * log(g/S))
  }
  return(sum1 / I_S)
}

# n = 5:

trees <- list(
  list(c(4,1), c(3,1), c(2,1), c(1,1)),
  list(c(4,1), c(2,2), c(1,1,1,1)),
  list(c(3,2), c(2,1,1,1), c(1,1))
)

H <- sapply(trees, H_index)

tree_counts <- c(2, 1, 3)

tree_counts %*% H / sum(tree_counts) # arithmetic mean of H

I_S <- sapply(trees, I_S_index)

tree_counts %*% I_S / sum(tree_counts) # arithmetic mean of I_S
1 / (tree_counts %*% (1 / I_S) / sum(tree_counts)) # harmonic mean of I_S

# n = 7:

tree4 <- list(c(5,2), c(4,1,1,1), c(3,1), c(2,1), c(1,1))
tree5 <- list(c(6,1), c(4,2), c(2,2,1,1), c(1,1,1,1))

identical(I_S_index(tree4), I_S_index(tree5)) # TRUE
identical(H_index(tree4), H_index(tree5)) # FALSE

########

s <- 0
for(i in 1:1e6) {
  X <- rexp(1)
  Y <- rexp(1)
  #s <- s + X / (X + Y)
  s <- s + min(X, Y) / max(X, Y)
}
s / 1e6

########

exp_phylo_ent <- function(k) {
  i <- 2:k
  res <- k / (k - 1) * sum(log(i) / (i * (i - 1)))
  return(res)
}
# limit ≈ 1.257747 ≈ log(3.517488)

n <- floor(2^seq(1, log2(5000), length = 100))
plot(sapply(n, exp_phylo_ent) ~ n, type = "l")
abline(h = 1.257747, lty = 2)

