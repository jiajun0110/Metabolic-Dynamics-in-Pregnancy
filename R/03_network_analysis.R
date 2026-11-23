## Intra-stage co-variation network analysis
##
## Input:  dynamics_matrix — numeric matrix [n_individuals x n_metabolites]
## Output: per-metabolite centrality metrics

library(Hmisc)
library(igraph)

build_network <- function(dynamics_matrix) {
  corr <- Hmisc::rcorr(dynamics_matrix, type = "spearman")
  p_adj <- matrix(p.adjust(as.vector(corr$P), method = "BH"),
                  nrow = nrow(corr$P), dimnames = dimnames(corr$P))
  
  adj <- corr$r
  adj[p_adj >= 0.05 | abs(corr$r) < 0.3] <- 0
  diag(adj) <- 0
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  
  data.frame(
    metabolite             = V(g)$name,
    degree                 = degree(g),
    normalized_degree      = degree(g) / (vcount(g) - 1),
    local_clustering_coeff = transitivity(g, type = "local")
  )
}

# Compare local clustering coefficients between transitions
# wilcox.test(net_early$local_clustering_coeff, net_late$local_clustering_coeff)
