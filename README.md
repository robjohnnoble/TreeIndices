# TreeIndices
 
Tree indices work following [Lemant *et al.* 2022](https://academic.oup.com/sysbio/article/71/5/1210/6567363).

`Calculate_indices`
* Calculates tree shape indices accounting for branch lengths
* Can read various formats of tree data and convert to phylo object

`ExpectedValues.R`
* Calculates various approximations of $E(J^1)$ under the Yule and uniform models
* Plots the expected values and approximations
* Calculates Sackin's index and corrected phylogenetic entropy of bifurcating trees with uniform branch lengths (given trees encoded in terms of branch weights)

`IndicesExample.R`
* Calculates and plots index values $E$, $M$, $J$ and $H$ for an ultrametric leafy tree with three equally sized leaves and two internal nodes, as a function of the depth of the non-root internal node $h$
* Generates random values of $h$, based on the coalescent model, and thus calculates average index values for random three-leaf trees under the coalescent model

`TreeGenerating.R`
* Generates sets of trees on $n$ leaves, with corresponding counts according to the Yule process
* Calculates $E(J^1)$ and approximations under the Yule and uniform processes
* Plots histograms of $J^1$ values for the Yule and uniform processes

Goals:
* Improve commenting and documentation
* Integrate this code with previous code from the `RUtreebalance` repository
