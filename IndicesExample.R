# Index values for an ultrametric leafy tree with 
# three equally sized leaves and two internal nodes. 
# The tree has height 1 and the non-root internal node is at depth h 
# (hence the four branch lengths are 1, h, 1-h, and 1-h).
# Below is an illustration of the case h = 1/2.

#      /\
#     /  \
#    /   /\
#   /   /  \
#  O   O    O

Ebar <- function(h) { # effective per-node entropy
  E1 <- 1/3*(2*h*(log(3) - log(2)) + log(3))
  if(h > 1/2) E2 <- 2/3*(1-h)*log(2)
  else E2 <- 2/3*(h*log(2) + (1-2*h)*log(3))
  return(E1 + E2)
  }
Mbar <- function(h) { # log effective node degree
  m1 <- 1/3*(h * (3*log(2) - log(3)) + log(3))
  if(h > 1/2) m2 <- 2/3*(1-h)*log(2)
  else m2 <- 2/3*(h*log(2) + (1-2*h)*log(3))
  return(m1 + m2)
}
Jbar <- function(h) h*(log2(3) - 5/3) + 1 # tree balance
Hbar <- function(h) log(3) - 2/3*h*log(2) # phylogenetic entropy

h <- seq(0,1,length=1e3+1)
wid <- 2
pdf("ThreeLeavesExample.pdf",width=4.2, height = 4.6)
plot(sapply(h, Ebar) ~ h, type = "l", xlim = c(0,1), ylim = c(0,1.2), col = "magenta", 
     xaxt = "n", yaxt = "n", ylab = "index value", lwd = wid)
axis(2, at = c(0, log(2), 1, log(3)), labels = c("0", "log 2", "1", "log 3"), las = 1)
axis(1, at = c(0, 1/2, 1), labels = c("0", "0.5", 1))
lines(sapply(h, Hbar) ~ h, col = "black", lty = 3, lwd = wid)
abline(h = 1, lty = 2, col = "grey")
abline(h = log(2), lty = 2, col = "grey")
abline(h = log(3), lty = 2, col = "grey")
abline(v = 0, lty = 2, col = "grey")
abline(v = 1/2, lty = 2, col = "grey")
abline(v = 1, lty = 2, col = "grey")
lines(sapply(h, Mbar) ~ h, col = "darkcyan", lwd = wid)
lines(sapply(h, Jbar) ~ h, col = "orange", lwd = wid)
legend("bottomleft", col = c("magenta", "darkcyan", "orange", "black"), 
       lty = c(1, 1, 1, 3), 
       legend = c("E (effective per-node entropy)", "M (log effective node degree)", "J (tree balance)", expression('H'["p"]*' (phylogenetic entropy)')), 
       bty = "n", lwd = wid)
dev.off()

p <- 9/10
q <- (1-p)/2
points(0, (-2*q*log(q) - p*log(p))/log(3), col = "orange")
points(1, (-q*log(q) - (p+q)*log(p+q))/log(2), col = "orange")
points(0, -2*q*log(q) - p*log(p), col = "magenta")
points(1, -q*log(q) - (p+q)*log(p+q), col = "magenta")
points(0, (-2*q*log(q) - p*log(p))/Mbar(0), col = "blue")
points(1, (-q*log(q) - (p+q)*log(p+q))/Mbar(1), col = "blue")

# The ratio Ebar / Mbar is generally not equal to Jbar 
# but is approximately equal for this particular tree
plot(sapply(h, Ebar) ~ h, type = "l", xlim = c(0,1), ylim = c(0.9,1), col = "white", lwd = wid)
lines(sapply(h, Jbar) ~ h, col = "orange", lwd = wid)
lines(sapply(h, Ebar) / sapply(h, Mbar) ~ h, col = "blue", lwd = wid)

##########

# Generate a random value of h, based on the coalescent model:
stoch_fn <- function() {
  t1 <- rexp(n = 1, rate = 3) # time during which there are exactly three lineages
  t2 <- rexp(n = 1, rate = 1) # time during which there are exactly two lineages
  h <- t1 / (t1 + t2)
  return(h)
}

# vector of random h values:
h_vec <- replicate(1e6, stoch_fn())
mean(h_vec) # ≈ 0.324

# average index values:
mean(sapply(h_vec, Ebar)) # ~0.8907
mean(sapply(h_vec, Mbar)) # ~0.9090
mean(sapply(h_vec, Jbar)) # ~0.9736
mean(sapply(h_vec, Hbar)) # ~0.9490



