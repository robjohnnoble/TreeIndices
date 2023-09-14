# Checks file extension, reads in tree and converts to phylo object
# Adds equal branch lengths of size 1 if none are given
read_convert <- function(file){

  suppressWarnings({
    tree <- try(ape::read.tree(file), silent = TRUE)
    if ("try-error" %in% class(tree)){
      tree <- try(ape::read.nexus(file), silent = TRUE)
      if ("try-error" %in% class(tree)){
        tree <- try(ape::read.tree(text=file), silent = TRUE)
        if ("try-error" %in% class(tree)){
          if (is.phylo(file)){
            tree <- file
          }else if(!(is.phylo(file))){
            return("Tree must be in Newick or NEXUS format, or be a phylo object.") 
          }
        }
      }
    }
  })  
  
  if (length(tree$edge.length) == 0){ # If no branch lengths, assign all to be one
    tree$edge.length <- rep(1, times = length(tree$edge[,1]))
  }

  return(tree)
}

# Function to calculate total abundance descending from each branch of a phylo object,
abundance_phylo <- function(tree, abundance_data) {
  
  node_labels <- append(tree$tip.label, tree$node.label)
  
  # Initialize a dictionary to store the total abundance for each branch
  total_abundance <- list()
  
  # Define a recursive function to calculate total abundance
  calculate_abundance <- function(node) {

    # Initialize the total abundance for the current node
    
    node_abundance <- abundance_data[which(abundance_data[,1] == node_labels[node]), 2]
    
    # Loop through each child of the current node
    for (child in tree$edge[tree$edge[, 1] == node, 2]) {
      
      # If the child is not a leaf, recursively calculate the total abundance of its subtree
      if (!(child %in% c(1:length(tree$tip.label)))) {
        child_abundance <- calculate_abundance(child)
        
      }
      # If the child is a leaf, use its abundance as the total abundance of its subtree
      else {
        child_abundance <- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
        total_abundance[[as.character(child)]] <<- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
      }
      
      # Add the total abundance of the child subtree to the current node's total abundance
      node_abundance <- node_abundance + child_abundance
    }
    
    # Store the total abundance for the current node
    total_abundance[[as.character(node)]] <<- node_abundance
    
    # Return the total abundance for the current node
    return(node_abundance)
  }
  
  # Start the recursive calculation from the root node
  calculate_abundance(tree$edge[1, 1])
  
  
  # Return the dictionary of total abundances
  return(total_abundance)
}

# Finds distance between two nodes when top_node is an ancestor of bottom_node
distance <- function(tree, top_node, bottom_node){
  curr_node <- bottom_node # Start at bottom node
  dist <- 0 # Initialise sum
  while (curr_node != top_node) {
    dist <- dist + tree$edge.length[which(tree$edge[,2] == curr_node)] # Sum distance
    curr_node <- tree$edge[which(tree$edge[,2] == curr_node),1] # Move to parent node
  }
  return(dist)
}

# Calculates all ancestor nodes of a given node, including itself
getAllAncestors <- function(tree,node){
  getAncestor <- function(tree, node){ # Returns immediate ancestor of given node
    i <- which(tree$edge[, 2] == node)
    return(tree$edge[i, 1])
  }
  root_node <- length(tree$tip.label) + 1 # Number assigned to root node
  ancestors <- c(node) # Include node as an ancestor of itself
  anc_node <- root_node + 1 # Assign so while loop will start
  if (node != root_node){
    while (anc_node > root_node){
      anc_node <- getAncestor(tree, node) # Ancestor of node
      ancestors <- append(ancestors, anc_node) # Append ancestor to list
      node <- anc_node # Move to ancestor
    }
  }
  return(ancestors)
}

# Calculates all descendant branches of a node, recording distance from node or "x"
get_descendant_branches <- function(tree, node){
  library(TreeTools)
  branches <- which(tree$edge[,1] == node) # Row index of direct descendant branches
  if (length(branches) == 0) { # If no descendant branches return empty matrix
      return(data.frame("Start_node" = numeric(), "End_node" = numeric(),
                        "Branch_length" = numeric(), "x" = numeric()))}
  else{
    descendants <- c() # Empty array
    for (i in 1:length(branches)){ # Calculates row index of all branches descendant from node
      descendants <- append(descendants, which(DescendantEdges(edge = branches[i],
                                                               tree$edge[,1], tree$edge[,2])))
    }
    sub_tree <- tree$edge[descendants,] # Select subtree from node
    x <- sapply(sub_tree[,2], distance, top_node=node, tree=tree) # Record distance between node and end of each branch
    df_branch_info <- data.frame("Start_node" = sub_tree[,1], "End_node" = sub_tree[,2],
                       "Branch_length" = tree$edge.length[descendants], 
                       "x" = x)
  }
  
  return(df_branch_info)
}

# Calculates T_i and S_i for a given node i in a phylogenetic tree
compute_T_i_S_i <- function(tree, node, abundances) {

  # Choose branches that are immediate descendants of node and record data in dataframe
  df_descen_info <- data.frame("Start_node" = tree$edge[tree$edge[, 1] == node, 1],
                               "End_node" = tree$edge[tree$edge[, 1] == node, 2],
                               "Branch_lengths" = tree$edge.length[which(tree$edge[, 1] == node)])

  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  
  T_i <- 0 # Initialise T_i sum
  prev_x <- 0 # Keep track of previous value of x
  for (j in 1:length(unique(df_descen_info[,"Branch_lengths"]))) { # For each branch as all could have different branch lengths
    abundance_sum <- 0 # Initialise sum
    if (nrow(df_descen_info) != 0){
      for (i in 1:nrow(df_descen_info)) {
        if (!(is.na(df_descen_info["End_node"][i,]))){
          end_node <- df_descen_info["End_node"][i,] # Select end node of current branch
          abundance_sum <- abundance_sum + abundances[[as.character(end_node)]] # Sum abundances
        }
      }
    }

    if (nrow(df_descen_info) != 0){
      DF_S_i <- rbind(DF_S_i, data.frame("S_i" = abundance_sum, "x" = min(df_descen_info["Branch_lengths"]))) # Append abundance and corresponding value of x
      
      T_i <- T_i + (abundance_sum * (min(df_descen_info["Branch_lengths"]) - prev_x)) # Sum T_i
      prev_x <- min(df_descen_info["Branch_lengths"]) # Update previous x
      
      df_descen_info <- df_descen_info[-(which(df_descen_info["Branch_lengths"] == min(df_descen_info["Branch_lengths"]))),] # Remove branch/s corresponding to x value
    }
  }

  return(list(T_i, DF_S_i))
}

# Calculate S_i_a for a given node i and ancestor a in a phylogenetic tree
calculate_S_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i) {
  
  temp <- function(a,b){return(a[[b]])}
  
  if (node != curr_ancestor){ # If not considering itself as ancestor  
    df_node_info <- get_descendant_branches(tree, curr_ancestor) # Get descendant branches from ancestor
    df_node_info <- df_node_info[df_node_info$x > h,] # Select branches that pass node
  
    curr_index <- which(df_node_info[,"x"] - df_node_info[,"Branch_length"] < h)
    df_node_info[curr_index,"Branch_length"] <- df_node_info[curr_index,"x"] - h  # For any branches that start before the node, correct branch lengths to be measured from node
    df_node_info$x <- df_node_info$x - h # Change distance to be distance from node to end of branch
  }else{ # Considering itself as ancestor
    df_node_info <- get_descendant_branches(tree, node) # All descendant branches from node in dataframe
  }

  df_node_info_ed <- df_node_info
  
  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  abund_list <- list() # Empty dictionary for abundances
  x <- 0 # Keep track of distance from node
  for (j in 1:nrow(df_node_info)) { # For each branch as all could have different branch lengths
    if (x != l_i){
      list_ab <- c() # Empty list to hold branch abundances on each iteration
      index_curr_branches <- which(df_node_info_ed[,"x"] %==% df_node_info_ed[,"Branch_length"]) # Row index of current branches

      if (nrow(df_node_info_ed) != 0){
        abund_list[[as.character(j)]] <- sapply(as.character(end_node <- df_node_info_ed[index_curr_branches,"End_node"]),
                                                temp, a = abundances) # Store all branch abundances present for this value of x
        
        abundance_sum <- sum(abund_list[[as.character(j)]])
        
        DF_S_i <- rbind(DF_S_i, data.frame("S_i" = abundance_sum,
                                           "x" = (min(df_node_info_ed[index_curr_branches,]["x"]) + x))) # Append abundance and corresponding value of x
        
        prev_x <- min(df_node_info_ed[index_curr_branches,]["x"]) # Update previous x
        x <- x + prev_x # Update x
        
        index_to_be_deleted <- which(df_node_info_ed[,"x"] %==% min(df_node_info_ed[index_curr_branches,]["x"])) # Indices of rows to be deleted 
        df_node_info_ed$x <- df_node_info_ed$x - prev_x # "Reset" x = 0 level
        df_node_info_ed[index_curr_branches,]["Branch_length"] <- df_node_info_ed[index_curr_branches,]["Branch_length"] - prev_x # Shorten current branch lengths by previous x
        df_node_info_ed <- df_node_info_ed[-index_to_be_deleted,] # Remove branch/s corresponding to value of x
        
      }
    }
  }

  return(list(DF_S_i, abund_list))
}

# Calculates E/J/M or all given i and a
calculate_EJM_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i, 
                              index_letter, q, individual){
  
  term <- function(a,b,log_base){-(a/b) * log(a/b, base = log_base)}
  
  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h, l_i) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list

  if (individual == FALSE){ # Create empty dataframe/s
    df_E_i_a <- data.frame("E_i_a" = numeric(), "x" = numeric()) 
    df_J_i_a <- data.frame("J_i_a" = numeric(), "x" = numeric())
    df_M_i_a <- data.frame("M_i_a" = numeric(), "x" = numeric())
  }else if (individual == TRUE){
    df_In_i_a <- data.frame("E_i_a" = numeric(), "x" = numeric()) 
  }


  for (k in 1:length(abund_list)){
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    
    if (individual == FALSE){
      m <- log(length(abund_vec)) # Calculate out-degree term
      e <- sum(sapply(abund_vec, term, b = S, log_base=exp(1))) # Calculate diversity term
      if (length(abund_vec) > 1 ){ # Removes case of only one branch in region
        j <- sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) # Calculate balance term
      }else if (length(abund_vec) == 1){
        j <- 1
      }
      # Store values and corresponding x
      df_M_i_a <- rbind(df_M_i_a, data.frame("M_i_a" = m, "x" = df_S_i_a[k, "x"]))
      df_E_i_a <- rbind(df_E_i_a, data.frame("E_i_a" = e, "x" = df_S_i_a[k, "x"]))
      df_J_i_a <- rbind(df_J_i_a, data.frame("J_i_a" = j, "x" = df_S_i_a[k, "x"]))
    }else if (individual == TRUE){
      if ((index_letter == "D")&(q == 0)){
        ind_val <- log(length(abund_vec))
      }else if ((index_letter == "D")&(q == 1)){
        ind_val <- sum(sapply(abund_vec, term, b = S, log_base=exp(1)))
      }else if (index_letter == "J"){
        if (length(abund_vec) > 1 ){ # Removes case of only one branch in region
          ind_val <- sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) # Calculate balance term
        }else if (length(abund_vec) == 1){
          ind_val <- 1
        }
      }
      df_In_i_a <- rbind(df_In_i_a, data.frame("In_i_a" = ind_val, "x" = df_S_i_a[k, "x"]))
    }
  }
  if (individual == FALSE){
    return(list("1DN" = df_E_i_a, "1JN" = df_J_i_a, "0DN" = df_M_i_a))
  }else if (individual == TRUE){
    return(df_In_i_a)
  }
}

# Calculate S_i_a by first attaching all of trees branches to root node/creates star tree
calculate_S_i_a_star <- function(tree, node, abundance_data, l_i) {
  
  library(fpCompare)
  
  temp <- function(a,b){return(a[[b]])}
  
  abundances <- abundance_phylo(tree, abundance_data)
  
  node_labels <- append(tree$tip.label, tree$node.label)
  
  df_node_info <- get_descendant_branches(tree, node) # Get descendant branches from ancestor
  df_node_info$x <- df_node_info$Branch_length
  df_node_info_ed <- df_node_info
  
  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  abund_list <- list() # Empty dictionary for abundances
  for (j in 1:length(unique(df_node_info$x))) { # For each branch as all could have different branch lengths
    if (nrow(df_node_info_ed) != 0){
      end_nodes <- df_node_info_ed$End_node
      abund_list[[as.character(j)]] <- sapply(as.character(end_nodes), temp, a = abundances) # Store all branch abundances present for this value of x
        #c(abundance_data[which(abundance_data[,1] %in% node_labels[end_nodes]), 2]) # Store all branch abundances present for this value of x
      
      abundance_sum <- sum(abund_list[[as.character(j)]])
      
      DF_S_i <- rbind(DF_S_i, data.frame("S_i" = abundance_sum,
                                         "x" = (min(df_node_info_ed$x)))) # Append abundance and corresponding value of x
      
      df_node_info_ed <- df_node_info_ed[-(which(df_node_info_ed$x == min(df_node_info_ed$x))),] # Remove branch/s corresponding to x value
    }
  }
  
  return(list(DF_S_i, abund_list))
}

# Calculates integral for each ancestor
calculate_integral <- function(tree, node, curr_ancestor, S_i, abundances, index_letter,
                               q, individual){

  if (length(S_i["x"] > 1)){
    l_i <- max(S_i["x"]) # Select longest branch length of a immediate descendant
  }else{
    l_i <- S_i["x"] # If only one value, select that
  }
  
  if (l_i == 0){ # h = 0 do not do integral
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }

  h <- distance(tree,curr_ancestor,node) # Distance between current node and ancestor
  
  # Assign distance to parent
  if (curr_ancestor == (length(tree$tip.label) + 1)){ # If ancestor is the root node
    d_parent <- l_i  # By definition its Inf, but integral is only nonzero to l_i
  }else{
    d_parent <- tree$edge.length[which(tree$edge[,2] == curr_ancestor)]
  }

  # Ancestor integral
  int_2 <- 0 # Initialise sum
  prev_x <- h # Keep track of x, h as integral is from h
  if (l_i <= h){ # Nodes furthest away child is closer than distance to ancestor
    int_2 <- 0
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }else{
    current_row <- which(S_i[,"x"] > h)[1] # Start sum over first x that "reaches" ancestor
    while (current_row <= length(S_i[,"x"])) { # Sum over all rows  of S_i
      x_s <- S_i[current_row,2] # Current value of x for S_i
      value_s <- S_i[current_row,1] # Current value of S_i for given x
      if (x_s > d_parent){
        int_2 <- int_2 + (value_s * (d_parent - prev_x)) # Integral
        prev_x <- x_s
        break # Exit while loop
      }else { # x_s < d_parent
        int_2 <- int_2 + (value_s * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row <- current_row + 1 # Move to next row
      }
    }
  }
  
  integral <- function(df){
    prev_x <- 0 # Keep track of x
    current_row_s <- 1 # Start sum over rows of S_i
    current_row_e <- 1 # Start sum over rows of "other value" (E/J/M)
    int_1 <- 0 # Initialise integral sum
    
    # E/J/M integral
    while ((current_row_s <= length(S_i[,"x"]))&(current_row_e <= length(df$x))){ # Sum over all rows of S_i
      x_s <- S_i[current_row_s,2] # Current value of x for S_i
      x_e <- df[current_row_e, 2] # Current value of x for E_i_a
      value_s <- S_i[current_row_s,1] # Current value of S_i for given x
      value_e <- df[current_row_e,1] # Current value of E_i_a for given x
      if (x_s < x_e){
        int_1 <- int_1 + (value_s * value_e * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
      }else if (isTRUE(all.equal(x_s, x_e))){
        int_1 <- int_1 + (value_s * value_e * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
        current_row_e <- current_row_e + 1 # Move to next row
      }else{ # x_e > x_s
        if (x_e < l_i){
          int_1 <- int_1 + (value_s * value_e * (x_e - prev_x)) # Sum integral
          prev_x <- x_e # Update previous x
        }else{ # x_e > l_i
          int_1 <- int_1 + (value_s * value_e * (l_i - prev_x)) # Sum integral
          prev_x <- x_e # Update previous x
        }
        current_row_e <- current_row_e + 1 # Move to next row
      }
    }
    return(int_1) 
  }
  
  if (any(int_2 != 0)){
    # # Selects desired function
    # if (individual == TRUE){
    #   if (index_letter == "E"){
    #     Sum_i_a <- calculate_E_i_a(tree,node,abundances,curr_ancestor, h, l_i) # Run function
    #   }else if(index_letter == "J"){
    #     Sum_i_a <- calculate_J_i_a(tree, node, abundances, curr_ancestor, h, l_i) # Run function 
    #   }else if (index_letter == "M"){
    #     Sum_i_a <- calculate_M_i_a(tree,node,abundances,curr_ancestor, h, l_i) # Run function
    #   }
    # }else{
    #   Sum_i_a <- calculate_EJM_i_a(tree, node, abundances, curr_ancestor, h, l_i)
    # }
    
    Sum_i_a <- calculate_EJM_i_a(tree, node, abundances, curr_ancestor, h, l_i, 
                                 index_letter, q, individual)
    
    if (individual == TRUE){
      int <- integral(Sum_i_a)*int_2
    }else{
      int <- sapply(Sum_i_a, integral)*int_2
    }
  }else if (any(int_2 == 0)){
    if (individual == TRUE){
      int <- 0
    }else{
      int <- c(0,0,0)
    }
  }

  return(int)
}

# Calculates E_i, J_i and M_i 
calculate_EJM_i <- function(tree, node, abundances, index_letter, q, individual){
  
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe

  EJM_i <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                  abundances=abundances, index_letter = index_letter, q = q,
                  individual=individual)

  
  if (individual == FALSE){
    ejm_i <- apply(EJM_i,1,sum)
    if (T_i != 0){ # h > 0
      E_i <- (1/T_i) * ejm_i[1]
      J_i <- (1/T_i) *ejm_i[2]
      M_i <- (1/T_i) *ejm_i[3]
    }else if (T_i == 0){ # h = 0
      E_i <- 0
      J_i <- 1
      M_i <- 0
    }
    v <- c(E_i, J_i, M_i) 
  }else if (individual == TRUE){
    ejm_i <- sum(EJM_i)
    if (T_i != 0){ # h > 0
      v <- (1/T_i) * ejm_i
    }else if ((T_i == 0)&(index_letter == "J")){ # h = 0
      v <- 1
    }else if (T_i == 0){
      v <- 0
    }
  }

  
  return(v)
}

# Calculates longitudinal or star mean for D for q=0 or 1, or J for q =1, 
# or all for one mean type
long_star <- function(file, node_abundances, mean_type, index_letter, q, individual){
  
  term <- function(a,b,log_base){-a * log(a/b, base = log_base)}
  
  index_letter <- toupper(index_letter)
  
  tree <- read_convert(file)
  node <- tree$edge[1,1] # Root node
  if (is.data.frame(node_abundances)){ # If tree has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))){ #if doesn't, assign leaves to be equally abundant, internal nodes size zero
    num_tips <- tree$edge[1,1] - 1
    tree$node.label<- as.character(c(1:num_tips))
    num_nodes <- tree$Nnode
    tree$node.label<- as.character(c((num_tips+1):(num_tips + num_nodes)))
    
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }
  
  # Chooses mean type
  if (toupper(mean_type) == "LONGITUDINAL"){
    S_i_a_res <- calculate_S_i_a(tree, node, abundances, node, 0, sum(tree$edge.length)) # Run function
  }else if (toupper(mean_type) == "STAR"){
    S_i_a_res <- calculate_S_i_a_star(tree, node, node_abundances, max(tree$edge.length)) # Run function
  }

  # Calculate value/s
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  if (individual == TRUE){# Initialise sum/s
    h <- 0
  }else if (individual == FALSE){
    h1 <- 0
    h0 <- 0
    j1 <- 0
  }
  T_S_sum <- 0 # Initialise sum
  prev_x <- 0 # Keep track of x
  for (k in 1:length(abund_list)){
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    x <- df_S_i_a[k, "x"] - prev_x # Assign interval length
    abund_vec <- unlist(abund_list[k], use.names = FALSE) # Branch abundances present in interval
    T_S_sum <- T_S_sum + (S * x) # H_bar

    if (individual == TRUE){ # One q value
      if (index_letter == "J"){ # J for q = 1
        if (length(abund_vec) != 1){
          h <- h + (sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) * x)
        }else if (length(abund_vec) == 1){
          h <- h + (1 * S * x)
        }
      }
      if (q == 1){ # Check order
        h <- h + (sum(sapply(abund_vec, term, b = S, log_base=exp(1))) * x)
      }else if (q == 0){
          h <- h + (S * x * log(length(abund_vec), base=exp(1)))
      }
      
    }else if (individual == FALSE){ # All q values
        h1 <- h1 + (sum(sapply(abund_vec, term, b = S, log_base=exp(1))) * x)
        h0 <- h0 + (S * x * log(length(abund_vec), base=exp(1)))
        if (length(abund_vec) != 1){
          j1 <- j1 + (sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) * x)
        }else if (length(abund_vec) == 1){
          j1 <- j1 + (1 * S * x)
        }
    }

    prev_x <- df_S_i_a[k,"x"] # Update x
  }

  # Normalise index/indices
  if (individual == TRUE){
    if (index_letter == "J"){
      H <- (h/T_S_sum) 
    }else{
      H <- exp((h/T_S_sum)) 
    }
  }else if (individual == FALSE){
      h <- c(h0, h1)
      H <- c(exp((h/T_S_sum)), j1/T_S_sum)
  }
  
  return(H)
}

# Calculates node-wise mean for D for q = 0 or 1, or for J for q = 1, or all
node <- function(file, node_abundances, index_letter, q, individual){
  tree <- read_convert(file)
  index_letter <- toupper(index_letter)
  
  if (is.data.frame(node_abundances)){ # Has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))){ #Doesn't, assign leaves to be equally abundance, internal nodes size zero
    num_tips <- tree$edge[1,1] - 1
    tree$node.label<- as.character(c(1:num_tips))
    num_nodes <- tree$Nnode
    tree$node.label<- as.character(c((num_tips+1):(num_tips + num_nodes)))
    
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }

  # Calculates h_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }

  nodes <- c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)) # All internal node numbers
  
  if (individual == FALSE){ # Calculates all node-averaged indices
    EJM <- sapply(nodes, calculate_EJM_i, tree=tree, abundances = abundances, 
                  individual = individual, q = q, index_letter = "D")
    
    D1N <- (1/T)*sum(EJM[1,])
    J1N <- (1/T)*sum(EJM[2,])
    D0N <- (1/T)*sum(EJM[3,])
    
    List <- list("D1N"= exp(D1N),"J1N" = J1N,"D0N" = exp(D0N))
    return(List)
  }else if (individual == TRUE){ # Calculates one index
    EJM <- sapply(nodes, calculate_EJM_i, tree=tree, abundances = abundances, 
                  individual = individual, index_letter = index_letter, q = q)
    if (index_letter == "J"){
      index <- (1/T)*sum(EJM)
    }else if (!(index_letter == "J")){
      index <- exp((1/T)*sum(EJM))
    }
    
    return(index)
  }
}

# Returns a list of all index values, key is letter for index with first number second
# e.g. 1JN has key J1N
all_indices <- function(file, node_abundances){
  
  node <- node(file, node_abundances, "D", 1, FALSE)
  star <- long_star(file, node_abundances, "Star", "D", 0, FALSE)
  long <- long_star(file, node_abundances, "Longitudinal", "D", 0, FALSE)

  
  values <- list("D0N"= node$D0N,"D1N" = node$D1N,"J1N" = node$J1N, "D0S" = star[1],
                 "D1S" = star[2], "J1S" = star[3], "D0L" = long[1], "D1L" = long[2],
                 "J1L" = long[3])
  return(values)
}
