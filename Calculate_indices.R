# Function to calculate total abundance descending from each branch of a phylo object,
abundance_phylo <- function(tree, abundances) {
  
  # Initialize a dictionary to store the total abundance for each branch
  total_abundance <- list()
  
  # Define a recursive function to calculate total abundance
  calculate_abundance <- function(node) {

    # Initialize the total abundance for the current node
    node_abundance <- abundances[node]
    
    # Loop through each child of the current node
    for (child in tree$edge[tree$edge[, 1] == node, 2]) {
      
      # If the child is not a leaf, recursively calculate the total abundance of its subtree
      if (!(child %in% tree$tip.label)) {
        child_abundance <- calculate_abundance(child)
      }
      
      # If the child is a leaf, use its abundance as the total abundance of its subtree
      else {
        child_abundance <- abundances[child]
        total_abundance[[as.character(child)]] <<- abundances[child]
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
  recursive <- function(t, n) {
    node_descendant_index <- which(t$edge[,1] == n) # Row index of immediate descendant branches of n
    if (length(node_descendant_index) == 0) { # If no descendant branches return empty matrix
      return(data.frame("Start_node" = numeric(),
                        "End_node" = numeric(),
                        "Branch_length" = numeric(),
                        "x" = numeric()))
    }
    df_branch_info <- data.frame("Start_node" = numeric(length(node_descendant_index)),
                                 "End_node" = numeric(length(node_descendant_index)),
                                 "Branch_length" = numeric(length(node_descendant_index)),
                                 "x" = numeric(length(node_descendant_index))) # Create empty dataframe
    for (i in 1:length(node_descendant_index)) {
      descendant_index <- node_descendant_index[i] # Row index of descendant
      df_branch_info[i,"Start_node"] <- t$edge[descendant_index,1] # Start node
      df_branch_info[i,"End_node"] <- t$edge[descendant_index,2] # End node
      df_branch_info[i,"Branch_length"] <- t$edge.length[descendant_index] # Branch length
      df_branch_info[i,"x"] <- distance(t,node,t$edge[descendant_index,2]) # Distance between node and end of descendant branch
      if (length(which(t$edge[,1] == t$edge[descendant_index,2])) != 0) { # If descendant node is not at a leaf 
        child_branches <- recursive(t, t$edge[descendant_index,2]) # Recursively calls function
        df_branch_info<- rbind(df_branch_info, child_branches)
      }
    }
    return(df_branch_info)
  }
  branches <- recursive(tree, node)
  return(branches)
}

# Calculates all descendant branches of an ancestor of a node, recording distance from ancestor or "x+h"
get_descendant_branches_anc <- function(tree, node, ancestor){
  recursive <- function(Tree, Node, Ancestor) {
    ances_descendant_index <- which((Tree$edge[,1] == Ancestor)) # Row index of immediate descendant branches of Ancestor
    if (length(which(ances_descendant_index == which(Tree$edge[,2] == node))) != 0){ # Only non-zero when Ancestor is parent of node
      ances_descendant_index <- ances_descendant_index[-which(ances_descendant_index == which(Tree$edge[,2] == Node))] # Delete branch connecting node to parent
    }
    if (length(ances_descendant_index) == 0) { # If no descendant branches return empty matrix
      return(data.frame("Start_node" = numeric(),
                        "End_node" = numeric(),
                        "Branch_length" = numeric(),
                        "x" = numeric()))
    }
    df_branch_info <- data.frame("Start_node" = numeric(length(ances_descendant_index)),
                                 "End_node" = numeric(length(ances_descendant_index)),
                                 "Branch_length" = numeric(length(ances_descendant_index)),
                                 "x" = numeric(length(ances_descendant_index))) # Create empty dataframe
    for (i in 1:length(ances_descendant_index)) {
      descendant_index <- ances_descendant_index[i] # Row index of descendant
      df_branch_info[i,"Start_node"] <- Tree$edge[descendant_index,1] # Start node
      df_branch_info[i,"End_node"] <- Tree$edge[descendant_index,2] # End node
      df_branch_info[i,"Branch_length"] <- Tree$edge.length[descendant_index] # Branch length
      df_branch_info[i,"x"] <- distance(Tree,ancestor,Tree$edge[descendant_index,2]) # # Distance between ancestor and end of descendant branch
      if (length(which(Tree$edge[,1] == Tree$edge[descendant_index,2])) != 0) { # If descendant node is not at a leaf
        child_branches <- recursive(Tree, Node, Tree$edge[descendant_index,2]) # Recursively calls function
        df_branch_info <- rbind(df_branch_info, child_branches)
      }
    }
    return(df_branch_info)
  }
  branches <- recursive(tree, node, ancestor)
  return(branches)
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
  for (j in 1:nrow(df_descen_info)) { # For each branch as all could have different branch lengths
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
calculate_S_i_a <- function(tree, node, abundances, curr_ancestor, h) {
  
  elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
  
  select_abundance <- function(end_node, abundances){
    abundance <- abundances[[as.character(end_node)]]
    return(abundance)
  }

  df_descen_info <- get_descendant_branches(tree, node) # All descendant branches from node in dataframe
  
  if (node != curr_ancestor){ # If not considering itself as ancestor  
    df_ances_info <- get_descendant_branches_anc(tree, node, curr_ancestor) # Get descendant branches from ancestor
    df_ances_info <- df_ances_info[df_ances_info$x > h,] # Select branches that pass node
    df_ances_info[which(df_ances_info[,"x"] - df_ances_info[,"Branch_length"] < h),"Branch_length"] <- df_ances_info[which(df_ances_info[,"x"] - df_ances_info[,"Branch_length"] < h),"x"] - h  # For any branches that start before the node, correct branch lengths to be measured from node
    df_ances_info$x <- df_ances_info$x - h # Change distance to be distance from node to end of branch
  }else{ # Considering itself as ancestor
    df_ances_info <- data.frame("Start_node" = numeric(),
                                 "End_node" = numeric(),
                                 "x" = numeric())
  }
  if (nrow(df_ances_info) != 0){# Case when there are ancestor branches that pass node
    df_node_info <- rbind(df_descen_info, df_ances_info) # Append dataframes together
    df_node_info_ed <- df_node_info # Create identical dataframe to delete values from
  }else{ # Case when there are no ancestor branches to consider
    df_node_info <- df_descen_info
    df_node_info_ed <- df_node_info # Create identical dataframe to delete values from
  }
  
  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  abund_list <- list() # Empty dictionary for abundances
  x <- 0 # Keep track of distance from node
  for (j in 1:nrow(df_node_info)) { # For each branch as all could have different branch lengths
    list_ab <- c() # Empty list to hold branch abundances on each iteration
    index_curr_branches <- elementwise.all.equal(df_node_info_ed[,"x"], df_node_info_ed[,"Branch_length"]) # Row index of current branches
    abundance_sum <- 0 # Initialise sum
    if (nrow(df_node_info_ed) != 0){
      for (i in 1:nrow(df_node_info_ed[index_curr_branches,])){
        end_node <- df_node_info_ed[index_curr_branches,][i,"End_node"] # Select end node of current branch
        abundance_sum <- abundance_sum + abundances[[as.character(end_node)]] # Sum abundances
      }
    }
    if (nrow(df_node_info_ed) != 0){
      abund_list[[as.character(j)]] <- sapply(df_node_info_ed[index_curr_branches,"End_node"], FUN = select_abundance, abundances=abundances) # Store all branch abundances present for this value of x
      DF_S_i <- rbind(DF_S_i, data.frame("S_i" = abundance_sum, "x" = (min(df_node_info_ed[index_curr_branches,]["x"]) + x))) # Append abundance and corresponding value of x
      
      prev_x <- min(df_node_info_ed[index_curr_branches,]["x"]) # Update previous x
      x <- x + prev_x # Update x
      
      index_to_be_deleted <- elementwise.all.equal(df_node_info_ed[,"x"], min(df_node_info_ed[index_curr_branches,]["x"])) # Boolean list where indices to be deleted are TRUE 
      df_node_info_ed$x <- df_node_info_ed$x - prev_x # "Reset" x = 0 level
      df_node_info_ed[index_curr_branches,]["Branch_length"] <- df_node_info_ed[index_curr_branches,]["Branch_length"] - prev_x # Shorten current branch lengths by previous x
      df_node_info_ed <- df_node_info_ed[!index_to_be_deleted,] # Remove branch/s corresponding to value of x

    }
  }
  
  return(list(DF_S_i, abund_list))
}

# Calculates diversity term given i and a
calculate_E_i_a <- function(tree, node, abundances, curr_ancestor, h){

  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list

  df_E_i_a <- data.frame("E_i_a" = numeric(), "x" = numeric()) # Create empty dataframe
  for (j in 1:length(abund_list)){
    e <- 0 # Initialise sum
    S <- df_S_i_a[j,"S_i"] # Select S_i_a
    abund_vec <- unlist(abund_list[j], use.names = FALSE)
      for (i in 1:length(abund_vec)){ # Sum over all branches in each region of x
        g <- abund_vec[i] # Branch abundance
        e <- e + -(g / S)*log(g/S) # Calculate diversity term
      }
    df_E_i_a <- rbind(df_E_i_a, data.frame("E_i_a" = e, "x" = df_S_i_a[j, "x"])) # Store diversity term and corresponding x
  }
  return(df_E_i_a)
}

# Calculates balance term given i and a
calculate_J_i_a <- function(tree, node, abundances, curr_ancestor, h){
  
  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list

  df_J_i_a <- data.frame("J_i_a" = numeric(), "x" = numeric()) # Create empty dataframe
  for (k in 1:length(abund_list)){
    j <- 0 # Initialise sum
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    for (i in 1:length(abund_vec)){ # Sum over all branches in each region of x
      g <- abund_vec[i] # Branch abundance
      if (length(abund_vec) > 1){ # Avoids NaN in dataframe
        j <- j + -(g / S)*log((g/S), base=length(abund_vec)) # Calculate balance term
      }
    }
    df_J_i_a <- rbind(df_J_i_a, data.frame("J_i_a" = j, "x" = df_S_i_a[k, "x"])) # Store diversity term and corresponding x
  }
  return(df_J_i_a)
}

# Calculates log out-degree term given i and a
calculate_M_i_a <- function(tree, node, abundances, curr_ancestor, h){
  
  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  
  df_M_i_a <- data.frame("M_i_a" = numeric(), "x" = numeric()) # Create empty dataframe
  for (k in 1:length(abund_list)){
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    m <- log(length(abund_vec)) # Calculate out-degree term
    df_M_i_a <- rbind(df_M_i_a, data.frame("M_i_a" = m, "x" = df_S_i_a[k, "x"])) # Store out-degree term and corresponding x
  }
  return(df_M_i_a)
}

# Calculates integral for each ancestor
calculate_integral <- function(tree, node, curr_ancestor, S_i, abundances, integral){

  if (length(S_i["x"] > 1)){
    l_i <- max(S_i["x"]) # Select longest branch length of a immediate descendant
  }else{
    l_i <- S_i["x"] # If only one value, select that
  }
  
  h <- distance(tree,curr_ancestor,node)
  # Selects desired function
  if (integral == "E"){
    Sum_i_a <- calculate_E_i_a(tree,node,abundances,curr_ancestor, h) # Run function
  }else if(integral == "J"){
    Sum_i_a <- calculate_J_i_a(tree, node, abundances, curr_ancestor, h) # Run function 
  }else if (integral == "M"){
    Sum_i_a <- calculate_M_i_a(tree,node,abundances,curr_ancestor,h) # Run function
  }

  if (curr_ancestor == (length(tree$tip.label) + 1)){ # If ancestor is the root node
    d_parent <- l_i  # By definition its Inf, but integral is only nonzero to l_i
  }else{
    d_parent <- tree$edge.length[which(tree$edge[,2] == curr_ancestor)]
  }

  prev_x <- 0 # Keep track of x
  current_row_s <- 1 # Start sum over rows of S_i
  current_row_e <- 1 # Start sum over rows of "other value" (E/J/M)
  int_1 <- 0 # Initialise integral sum
  
  # E/J/M integral
  while (current_row_s <= length(S_i[,"x"])){ # Sum over all rows of S_i
    x_s <- S_i[current_row_s,2] # Current value of x for S_i
    x_e <- Sum_i_a[current_row_e, 2] # Current value of x for E_i_a
    value_s <- S_i[current_row_s,1] # Current value of S_i for given x
    value_e <- Sum_i_a[current_row_e,1] # Current value of E_i_a for given x
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

  int_2 <- 0 # Initialise sum
  prev_x <- h # Keep track of x, h as integral is from h
  if (l_i <= h){ # Nodes furthers away child is closer than distance to ancestor
    int_2 <- 0
  }else{
    current_row <- which(S_i[,"x"] > h)[1] # Start sum over first x that "reaches" ancestor
      while (current_row <= length(S_i[,"x"])) { # Sum over all rows  of S_i
        x_s <- S_i[current_row,2] # Current value of x for S_i
        value_s <- S_i[current_row,1] # Current value of S_i for given x
        if (x_s > d_parent){ # First descendant branch passes ancestors parent
          int_2 <- int_2 + (value_s * (d_parent - h)) # Integral
          prev_x <- x_s
          break # Exit while loop
        }else { # x_s < d_parent
          int_2 <- int_2 + (value_s * (x_s - prev_x)) # Sum integral
          prev_x <- x_s # Update previous x
          current_row <- current_row + 1 # Move to next row
        }
    }
  }
  integral <- int_1*int_2
  return(integral)
}

# Calculates E_i
calculate_E_i <- function(tree, node, abundances){
  
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node

  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe

  all_int <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i, 
                    abundances=abundances, integral = "E") # Vector of integrals for all ancestors

  E_i <- (1/T_i)*(sum(all_int)) # Calculate E_i
  
  return(E_i)
}

# Calculates J_i
calculate_J_i <- function(tree, node, abundances){
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe
  
  all_int <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                    abundances=abundances, integral = "J") # Vector of integrals for all ancestors
  
  J_i <- (1/T_i)*(sum(all_int)) # Calculate J_i
  
  return(J_i)
}

# Calculates M_i
calculate_M_i <- function(tree, node, abundances){
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe

  all_int <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                    abundances=abundances, integral = "M") # Vector of integrals for all ancestors
  
  M_i <- (1/T_i)*(sum(all_int)) # Calculate M_i
  
  return(M_i)
}

# Calculates E_i, J_i and M_i 
calculate_EJM_i <- function(tree, node, abundances){
  
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe
  
  E_i <- (1/T_i)*sum(sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i, 
                    abundances=abundances, integral = "E")) # Vector of integrals for all ancestors
  J_i <- (1/T_i)*sum(sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                    abundances=abundances, integral = "J")) # Vector of integrals for all ancestors
  M_i <- (1/T_i)*sum(sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                    abundances=abundances, integral = "M")) # Vector of integrals for all ancestors
  
  List <- list("E"= E_i,"J" = J_i,"M" = M_i)
  return(List)
}

# Calculates E
E <- function(tree, node_abundances){
  abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  
  # Calculates T_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }

  E_value <- (1/T) * sum(sapply(c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)), 
                              calculate_E_i, tree=tree,abundances=abundances))
  return(E_value)
}

# Calculates J
J <- function(tree, node_abundances){
  abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  
  # Calculates T_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }
  
  J_value <- (1/T) * sum(sapply(c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)), 
                                calculate_J_i, tree=tree,abundances=abundances))
  return(J_value)
}

# Calculates M (average log node out-degree)
M <- function(tree, node_abundances){
  abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  
  # Calculates T_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }
  
  M_value <- (1/T) * sum(sapply(c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)), 
                                calculate_M_i, tree=tree,abundances=abundances))
  return(M_value)
}

# Calculates H (corrected phylogenetic entropy)
H <- function(tree, node_abundances){
  
  node <- tree$edge[1,1] # Root node
  abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  S_i_a_res <- calculate_S_i_a(tree, node, abundances, node, 0) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  
  T_S_sum <- 0 # Initialise sum
  h <- 0 # Initialise sum
  prev_x <- 0 # Keep track of x
  for (k in 1:length(abund_list)){
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    x <- df_S_i_a[k, "x"] - prev_x # Assign interval length
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    T_S_sum <- T_S_sum + (S * x)
    for (i in 1:length(abund_vec)){ # Sum over all branches in each region of x
      g <- abund_vec[i] # Branch abundance
      h <- h + -(g * log(g/S) * x) # Calculate h term
    }
    prev_x <- df_S_i_a[k,"x"] # Update x
  }
  H <- (h/T_S_sum) # Normalise H
  return(H)
}

# Calculates E, J and M
EJM <- function(tree, node_abundances){
  abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  
  # Calculates T_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }
  
  E <- c()
  J <- c()
  M <- c()
  for (i in (length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)){
    ind <- calculate_EJM_i(tree,i,abundances)
    E <- append(E, ind$E)
    J <- append(J, ind$J)
    M <- append(M, ind$M)
  }
  List <- list("E"= (1/T) * sum(E),"J" = (1/T) * sum(J),"M" = (1/T) * sum(M))
  return(List)
}

# Returns a list of all index values, key is letter for index
all_indices <- function(tree, node_abundances){
  ejm <- EJM(tree, node_abundances)
  h <- H(tree, node_abundances)
  values <- list("E"= ejm$E,"J" = ejm$J,"M" = ejm$M,"H" = h)
  return(values)
}