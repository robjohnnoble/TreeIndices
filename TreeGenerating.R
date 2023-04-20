library(ape)
library(TreeTools)
library(treebalance)
library(rlist)
library(tictoc)
library(ggplot2)
library(dplyr)

# create a three-node, two-leaf tree:
cherry_tree <- function() {
  tree <- list()
  tree$edge <- matrix(c(3,1, 3,2), nrow = 2, byrow = TRUE)
  tree$tip.label <- c("t1", "t2")
  tree$Nnode <- 1
  class(tree) <- "phylo"
  attr(tree, "order") <- "cladewise"
  tree$count <- 1
  return(tree)
}

# given a set of trees with n leaves, create a set of trees with n+1 leaves:
next_set <- function(tree_list) {
  new_list <- list()
  index <- 1
  for(i in 1:length(tree_list)) {
    tree <- tree_list[[i]]
    ntips <- length(tree$tip.label)
    for(j in 1:ntips) {
      new_tree <- Renumber(AddTip(tree, j, paste0("t", ntips + 1)))
      new_tree$tip.label <- paste0("t", 1:(ntips + 1))
      new_tree$key <- paste(as.numeric(new_tree$edge), collapse = " ")
      new_list[[index]] <- new_tree
      index <- index + 1
    }
  }
  return(new_list)
}

# given a set of trees, replace each set of duplicate trees with a single tree and its count
# this reduces the number of trees (e.g. from 5040 to 429 trees on 8 leaves) 
# at the cost of extra processing time
add_counts <- function(tree_list) {
  all_edges <- unlist(lapply(tree_list, function(x) x$key)) # vector of all keys
  unique_edges <- unique(all_edges) # vector of unique keys
  unique_trees <- list()
  
  match <- function(a, b, i, j) a[i] == b[j]
  
  for(i in 1:length(unique_edges)) {
    edges <- unique_edges[[i]]
    count <- 0
    for(j in 1:length(all_edges)) {
      max_j <- length(all_edges)
      if(unique_edges[i] == all_edges[j]) {
        unique_trees[[i]] <- tree_list[[j]]
        count <- count + tree_list[[j]]$count
        if(j < max_j) {
          if(unique_edges[i] == all_edges[max_j]) count <- count + tree_list[[max_j]]$count
          all_edges[j] <- all_edges[max_j]
          all_edges <- all_edges[-max_j]
          tree_list[[j]] <- tree_list[[max_j]]
          tree_list[[max_j]] <- NULL
        } else {
          all_edges <- all_edges[-max_j]
          tree_list[[max_j]] <- NULL
        }
      }
      if(j >= length(all_edges)) break
    }
    unique_trees[[i]]$count <- count
  }
  print(as.numeric(lapply(unique_trees, function(x) x$count)))
  return(unique_trees)
}

# returns (in order):
# * E(I_S) for Yule process
# * n log2(n) / E(J^1) for Yule process
# * E(I_S) for uniform process
# * n log2(n) / E(J^1) for uniform process
get_index_values <- function(tree_list) {
  S <- as.numeric(lapply(tree_list, sackinI))
  counts <- as.numeric(lapply(tree_list, function(x) x$count))
  
  IS_Yule <- weighted.mean(S, counts)
  J1_Yule <- sum(counts) / sum(counts/S)
  IS_Unif <- mean(S)
  J1_Unif <- 1 / mean(1/S)
  
  return(c(IS_Yule = IS_Yule, J1_Yule = J1_Yule, IS_Unif = IS_Unif, J1_Unif = J1_Unif))
}

# returns (in order):
# * n log2(n) / E(I_S) for Yule process
# * E(J^1) for Yule process
# * n log2(n) / E(I_S) for uniform process
# * E(J^1) for uniform process
get_reciprocal_index_values <- function(tree_list, n) {
  numer <- n * log2(n)
  indices <- get_index_values(tree_list)
  return(numer / indices)
}

# plot histogram of index values
plot_hist <- function(tree_list, n) {
  
  numer <- n * log2(n)
  S <- as.numeric(lapply(tree_list, sackinI))
  if((max(S) == min(S))) stop("Only one index value")
  S <- numer / S
  counts_unif <- rep(1, length(S))
  counts_Yule <- as.numeric(lapply(tree_list, function(x) x$count))
  df <- data.frame(S = rep(S, 2), 
                   counts = c(counts_unif, counts_Yule), 
                   model = c(rep("uniform", length(S)), rep("Yule", length(S)))
                   )
  
  dfs <- group_by(df, S, model) %>% 
    summarise(counts = sum(counts)) %>% 
    ungroup()
  
  dfs <- group_by(dfs, model) %>% 
    mutate(freq = counts / sum(counts)) %>% 
    ungroup()
  
  indices <- get_reciprocal_index_values(tree_list, n)
  
  dodge_width <- (max(S) - min(S)) / 25
  
  ggplot(dfs, aes(x=S, group = model, color = model)) +
    geom_point( aes(y=freq), position = position_dodge(width = dodge_width) ) +
    geom_linerange(aes(ymax=freq), ymin=0, position = position_dodge(width = dodge_width)) + 
    geom_vline(xintercept = indices["IS_Yule"], linetype = "dashed", color = "red") + 
    geom_vline(xintercept = indices["J1_Yule"], linetype = "dotted", color = "red") + 
    geom_vline(xintercept = indices["IS_Unif"], linetype = "dashed", color = "blue") + 
    geom_vline(xintercept = indices["J1_Unif"], linetype = "dotted", color = "blue") + 
    scale_colour_manual(values = c("blue", "red")) + 
    labs(x = "index value", y = "frequency") + 
    theme_classic()
}

# initial list containing the unique 2-leaf tree:
l2 <- list(cherry_tree())

# create lists of trees with n>2 leaves:
l3 <- add_counts(next_set(l2))
l4 <- add_counts(next_set(l3))
l5 <- add_counts(next_set(l4))
l6 <- add_counts(next_set(l5))
l7 <- add_counts(next_set(l6))
l8 <- add_counts(next_set(l7))
l9 <- add_counts(next_set(l8))
tic() # this one takes about 25 seconds
l10 <- add_counts(next_set(l9)) # there are 4862 trees on 10 leaves (takes up 12.8 MB)
toc()
tic() # this one takes about XX seconds
l11 <- next_set(l10)
l11 <- add_counts(next_set(l11)) # there are XXX trees on 10 leaves (takes up XXX MB)
toc()

# index values:
get_index_values(l10)

# histogram:
pdf("HistogramTenLeaves.pdf", width = 5, height = 3)
plot_hist(l10, 10)
dev.off()

