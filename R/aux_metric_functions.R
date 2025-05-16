
# a suite of functions to obtain metrics from networks either
# in matrix or edge list format

# -------------------------------------------------------------------------
library(moments)

# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)
# -------------------------------------------------------------------------
# richness functions

richness.matrix <- function(A){
  rich <- NA
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      rich <- nrow(A)
    }else{
      rich <- length(unique(c(colnames(A),rownames(A))))
    }
  }else{
    if(nrow(A)==ncol(A)){
      rich <- nrow(A)
    }else{
      rich <- nrow(A) + ncol(A)
    }
  }# if not null names
  return(rich)
}

richness.edge.list <- function(el){
  rich <- length(unique(c(el$node_from,el$node_to)))
  return(rich)
}

# -------------------------------------------------------------------------
# check whether a network is bipartite

bipartiteness.matrix <- function(A){
  r <- unique(rownames(A))
  co <- unique(colnames(A))
  common.nodes <- dplyr::intersect(r,co)
  bipartiteness <- ifelse(length(common.nodes)>0,FALSE,TRUE)
  return(bipartiteness)
}

bipartiteness.edge.list <- function(el){
  r <- unique(el$node_from)
  co <- unique(el$node_to)
  common.nodes <- dplyr::intersect(r,co)
  bipartiteness <- ifelse(length(common.nodes)>0,FALSE,TRUE)
  return(bipartiteness)
}

# connectance -------------------------------------------------------------
# for now, binary
# connectance.edge.list <- function(el, diag = F){
#   A <- matrix_from_edge_list(el)
#   c <- connectance(interaction.matrix = A,quant = F)
#   return(c)
# }

connectance.matrix <- function(A, diago = F){
  return(MultitrophicFun::connectance(A,quant = F, diago = diago))
}

# degree distribution metrics ---------------------------------------------
degree.metrics.matrix <- function(A, directed = F){
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      degree.vec <- numeric(length = nrow(A))
      
      for(i in 1:nrow(A)){
        if(directed){
          degree.vec[i] <- sum(A[i,] != 0 | A[,i] != 0)
        }else{
          degree.vec[i] <- sum(A[i,] != 0)
        }
      }
      
    }else{
      degree.vec <- numeric(length = (nrow(A) + ncol(A)))
      for(i in 1:nrow(A)){
        degree.vec[i] <- sum(A[i,] != 0)
      }
      for(j in 1:ncol(A)){
        degree.vec[nrow(A)+j] <- sum(A[,j] != 0)
      }
    }
  }else{
    if(nrow(A)==ncol(A)){
      degree.vec <- numeric(length = nrow(A))
      
      for(i in 1:nrow(A)){
        if(directed){
          degree.vec[i] <- sum(A[i,] != 0 | A[,i] != 0)
        }else{
          degree.vec[i] <- sum(A[i,] != 0)
        }      
      }
    }else{
      degree.vec <- numeric(length = (nrow(A) + ncol(A)))
      for(i in 1:nrow(A)){
        degree.vec[i] <- sum(A[i,] != 0)
      }
      for(j in 1:ncol(A)){
        degree.vec[nrow(A)+j] <- sum(A[,j] != 0)
      }
    }
  }# if not null names
  
  return(list(mean.degree = mean(degree.vec),
              sd.degree = sd(degree.vec),
              skewness.degree = moments::skewness(degree.vec),
              kurtosis.degree = moments::kurtosis(degree.vec)
              ))
  
}

# degree distribution -----------------------------------------------------

degree.matrix <- function(A, directed = F){
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      degree.vec <- numeric(length = nrow(A))
      
      for(i in 1:nrow(A)){
        if(directed){
          degree.vec[i] <- sum(A[i,] != 0 | A[,i] != 0)
        }else{
          degree.vec[i] <- sum(A[i,] != 0)
        }
      }
      
    }else{
      degree.vec <- numeric(length = (nrow(A) + ncol(A)))
      for(i in 1:nrow(A)){
        degree.vec[i] <- sum(A[i,] != 0)
      }
      for(j in 1:ncol(A)){
        degree.vec[nrow(A)+j] <- sum(A[,j] != 0)
      }
    }
  }else{
    if(nrow(A)==ncol(A)){
      degree.vec <- numeric(length = nrow(A))
      
      for(i in 1:nrow(A)){
        if(directed){
          degree.vec[i] <- sum(A[i,] != 0 | A[,i] != 0)
        }else{
          degree.vec[i] <- sum(A[i,] != 0)
        }      
      }
    }else{
      degree.vec <- numeric(length = (nrow(A) + ncol(A)))
      for(i in 1:nrow(A)){
        degree.vec[i] <- sum(A[i,] != 0)
      }
      for(j in 1:ncol(A)){
        degree.vec[nrow(A)+j] <- sum(A[,j] != 0)
      }
    }
  }# if not null names
  return(vegan::diversity(degree.vec,index = "shannon"))
}

degree.edge.list <- function(el){
  degree.vec <- as.numeric(table(c(el$node_from,el$node_to)))
  return(vegan::diversity(degree.vec,index = "shannon"))
}

# -------------------------------------------------------------------------
# link density

# defined as "marginal totals of weighted diversity of interactions per species"
# (bipartite documentation of "linkage density")
# this includes diagonal links, and considers in- and out- links
# this naturally allows for directed/undirected links in unipartite matrices

link.density.matrix <- function(A, directed = T){
  
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      links.per.sp <- apply(A,1,FUN = function(x) sum(x!=0,na.rm = T))
    }else{
      links.per.sp.rows <- apply(A,1,FUN = function(x) sum(x!=0,na.rm = T))
      links.per.sp.cols <- apply(A,2,FUN = function(x) sum(x!=0,na.rm = T))
      links.per.sp <- c(links.per.sp.rows,links.per.sp.cols)
      
    }
  }else{
    if(nrow(A)==ncol(A)){

      links.per.sp <- apply(A,1,FUN = function(x) sum(x!=0,na.rm = T))
      
    }else{

      links.per.sp.rows <- apply(A,1,FUN = function(x) sum(x!=0,na.rm = T))
      links.per.sp.cols <- apply(A,2,FUN = function(x) sum(x!=0,na.rm = T))
      links.per.sp <- c(links.per.sp.rows,links.per.sp.cols)
    }
  }# if not null names
  return(mean(links.per.sp))
}

link.density.edge.list <- function(el){
  library(tidyverse)
  # make sure no duplicated links
  temp.el <- el
  temp.el$dup <- FALSE
  temp.el$flag <- FALSE
  for(i in 1:nrow(temp.el)){
    if(!temp.el$flag[i]){
      replicas <- which((temp.el$node_from == temp.el$node_from[i] &
                        temp.el$node_to == temp.el$node_to[i]) | 
                          (temp.el$node_to == temp.el$node_from[i] &
                             temp.el$node_from == temp.el$node_to[i]))
      if(length(replicas)>1){
        temp.el$dup[replicas[2:length(replicas)]] <- TRUE
        temp.el$flag[replicas[2:length(replicas)]] <- TRUE
      }
      temp.el$flag[i] <- TRUE
    }
  }
  
  temp.el <- subset(temp.el, !dup)
  
  links.per.sp.from <- temp.el %>%
    group_by(node_from) %>%
    summarise(links.from = n())
  links.per.sp.to <- temp.el %>%
    group_by(node_to) %>%
    summarise(links.to = n())
  
  all.links <- full_join(links.per.sp.from,links.per.sp.to,by = c("node_from" = "node_to"))
  all.links[is.na(all.links)] <- 0
  links.per.sp <- all.links$links.from + all.links$links.to
  # library(tidyverse)
  # links.per.sp.from <- el %>% 
  #   group_by(node_from) %>%
  #   summarise(num.links = n())
  return(mean(links.per.sp,na.rm=T))
}

# -------------------------------------------------------------------------
# centralization

# eigenvector centrality at the network level (see the help of e.g. centr_eigen,
# and see Delmas et al. 2018 for a brief summary)

centralization.matrix <- function(A, centrality.type = c("degree",
                                                         "betweenness",
                                                         "eigen")){
  # build igraph object depending on whether A is bipartite
  # or, if no info is given in names, guess its bipartiteness
  # by assuming that square matrices are unipartite.
  # Risk: bipartite nets with ncol = nrow will be taken as unipartite.
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      g <- graph_from_adjacency_matrix(A,mode = "undirected")
    }else{
      g <- graph_from_incidence_matrix(A)
    }
  }else{
    if(nrow(A)==ncol(A)){
      g <- graph_from_adjacency_matrix(A,mode = "undirected")
    }else{
      g <- graph_from_incidence_matrix(A)
    }
  }# if not null names
  
  if(centrality.type == "eigen"){
    # my.cent <- igraph::centr_eigen(g, directed = TRUE)$centralization
    eig.cent <- eigen_centrality(g,directed = F)$vector
    my.cent <- DescTools::Gini(eig.cent)
  }else if(centrality.type == "degree"){
    my.cent <- degree.centralization.graph(g)
  }else if(centrality.type == "betweenness"){
    my.cent <- igraph::centr_betw(g, directed = TRUE)$centralization
  }
  return(my.cent)
}


# -------------------------------------------------------------------------
# auxiliary functions for degree centralization

calc_degree_centrality <- function(g) {
  components <- igraph::components(g)$membership
  centralities <- list()
  
  for(comp in unique(components)) {
    subgraph <- induced_subgraph(g, which(components == comp))
    centrality <- igraph::degree(subgraph, mode = "all")
    centralities[[as.character(comp)]] <- centrality
  }
  
  return(centralities)
}

# this accounts for networks with more than one component,
# in that case, the centralization of the different components are calculated
# separatedly and aggregated for the final metric
# the formula for centralization comes from Butts 2006, Social Networks 28:283-296
degree.centralization.graph <- function(g){
  total_nodes <- length(V(g))
  centralities <- calc_degree_centrality(g)
  
  all_values <- unlist(centralities)
  max_centrality <- max(all_values)
  centralization <- sum(max_centrality - all_values)
  
  norm.centralization <- centralization / ((total_nodes - 1) * (total_nodes - 2))
  return(norm.centralization)
}

centralization.edge.list <- function(el, centrality.type = c("degree",
                                                             "betweenness",
                                                             "eigen")){
  g <- graph_from_data_frame(el[,c("node_from","node_to")])
  if(centrality.type == "eigen"){
    my.cent <- igraph::centr_eigen(g, directed = TRUE)$centralization
  }else if(centrality.type == "degree"){
    my.cent <- degree.centralization.graph(g)
  }else if(centrality.type == "betweenness"){
    my.cent <- igraph::centr_betw(g, directed = TRUE)$centralization
  }
  return(my.cent)
}

# -------------------------------------------------------------------------
# nestedness

# using maxnodf 
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13545

# or using bipartite - conceptually nestedness can also be calculated for a 
# unipartite network

nestedness.matrix <- function(A){
  # maxnodf::NODFc(A)
  bipartite::nested(A,method = "NODF")
}

nestedness.edge.list <- function(el){
  A <- matrix_from_edge_list(el)
  bipartite::nested(A,method = "NODF")
  # maxnodf::NODFc(A)
}

# -------------------------------------------------------------------------
# interaction overlap

interaction.overlap.matrix <- function(A, method = "horn"){
  
  dis.matrix <- as.matrix(vegan::vegdist(A,method = method))
  overlap.matrix <- 1 - dis.matrix
  
  sp.names <- as.character(1:nrow(A))
  if(!is.null(rownames(A))){
    sp.names <- rownames(A)
  }
  
  overlap.pairs <- list()
  
  for(i in 1:nrow(overlap.matrix)){
    for(j in 1:ncol(overlap.matrix)){
      overlap.pairs[[length(overlap.pairs)+1]] <-
        data.frame(sp1 = sp.names[i],sp2 = sp.names[j], overlap = overlap.matrix[i,j])
    }
  }
  overlap.df <- dplyr::bind_rows(overlap.pairs)
  return(overlap.df)
}

interaction.overlap.edge.list <- function(el, method = "horn"){
  A <- matrix_from_edge_list(el)
  return(interaction.overlap.matrix(A,method))
}

# -------------------------------------------------------------------------
# reshuffle interactions
reshuffle.edge.list <- function(el){
  A <- matrix_from_edge_list(el)
  shuff.A <- reshuffle_matrix(A)
  
  shuff.list <- list()
  for(i.row in 1:nrow(shuff.A)){
    for(i.col in 1:ncol(shuff.A)){
      if(shuff.A[i.row,i.col] != 0){
        my.link <- data.frame(node_from = colnames(A)[i.col], 
                              node_to = rownames(A)[i.row])
        shuff.list[[length(shuff.list)+1]] <- my.link
      }
    }# for i.col
  }# for i.row
  
  shuff.el <- bind_rows(shuff.list)
  return(shuff.el)
  # my.sp <- unique(c(el$node_from,el$node_to))
  # n.links <- nrow(el)
  # shuff.el <- el
  # shuff.el$node_from <- sample(my.sp,size = nrow(el),replace = T)
  # shuff.el$node_to <- sample(my.sp,size = nrow(el),replace = T)
  return(shuff.el)
}

# -------------------------------------------------------------------------
# igraph modularity

modularity.matrix <- function(A, modularity.type = c("edge_betweenness",
                                                     # "leading_eigen",
                                                     "infomap")){
  
  # build igraph object depending on whether A is bipartite
  # or, if no info is given in names, guess its bipartiteness
  # by assuming that square matrices are unipartite.
  # Risk: bipartite nets with ncol = nrow will be taken as unipartite.
  if(!is.null(colnames(A)) & !is.null(rownames(A))){
    if(identical(rownames(A),colnames(A))){
      my.graph <- graph_from_adjacency_matrix(A)
    }else{
      my.graph <- graph_from_incidence_matrix(A)
    }
  }else{
    if(nrow(A)==ncol(A)){
      my.graph <- graph_from_adjacency_matrix(A)
    }else{
      my.graph <- graph_from_incidence_matrix(A)
    }
  }# if not null names
  
  if(modularity.type == "edge_betweenness"){
    my.communities <- igraph::cluster_edge_betweenness(my.graph,directed = FALSE)
  # }else if(modularity.type == "leading_eigen"){
  #   my.communities <- igraph::cluster_leading_eigen(my.graph)
  }else if(modularity.type == "infomap"){
    my.communities <- igraph::cluster_infomap(my.graph)
  }
  
  my.mod <- igraph::modularity(my.communities)
  return(my.mod)
}

modularity.edge.list <- function(el, modularity.type = c("edge_betweenness",
                                                         # "spinglass",
                                                     # "leading_eigen",
                                                     "infomap")){
  my.graph <- igraph::graph_from_data_frame(d = el[,c("node_from","node_to")])
  
  if(modularity.type == "edge_betweenness"){
    my.communities <- igraph::cluster_edge_betweenness(my.graph,directed = FALSE)
  # }else if(modularity.type == "leading_eigen"){
  #   my.communities <- igraph::cluster_leading_eigen(my.graph)
  }else if(modularity.type == "infomap"){
    my.communities <- igraph::cluster_infomap(my.graph)
  # }else if(modularity.type == "spinglass"){
  #   my.communities <- igraph::cluster_spinglass(my.graph)
  }
  
  my.mod <- igraph::modularity(my.communities)
  return(my.mod)
}

# -------------------------------------------------------------------------

infomap.modularity.matrix <- function(A,
                               allowed.interactions.mat = NULL,
                               damping = .85, 
                               pr.algo = "prpack",
                               temp.dir = NULL){
  g  <- graph.adjacency(A,weighted=TRUE)
  
  # To run infomap it seems that weights should be non-negative
  my.edge.list <- get.data.frame(g) %>% mutate(weight=abs(weight))
  
  nodes <- my.edge.list$from %>% unique()
  
  nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
  
  # Preparing edge.list and nodes.ID to run infomap
  
  my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
    left_join(nodes.ID,by="species") %>% dplyr::select(-species) %>%
    rename(from=node_id,species=to) %>%
    left_join(nodes.ID,by="species") %>%
    dplyr::select(-species) %>% rename(to=node_id) %>% 
    dplyr::select(from,to,weight)
  
  nodes.ID2 <- nodes.ID %>% rename(node_name=species)
  
  infomap_mono <- create_monolayer_object(x = my.edge.list,
                                          directed = T,
                                          bipartite = F,
                                          node_metadata = nodes.ID2[,c(2,1)])
  
  # infomap_mono2 <- infomap_mono
  # infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
  
  # run Infomap
  if(!is.null(temp.dir)){
    # temp_dir <- paste("R/temp/d",i.id,sep="")
    # temp_dir <- "R/temp/d1"
    dir.create(temp.dir)
    modules_relax_rate <- run_infomap_monolayer3(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 temp_dir = temp.dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
    unlink(temp.dir, recursive = TRUE)
  }else{
    modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 # temp_dir = temp_dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
  }
  
  
  my.modularity <- modules_relax_rate[["L"]] # modularity in bits
  
  # Extract module information
  modules <- modules_relax_rate$modules %>%
    dplyr::select(node_id,module_level1) %>%
    rename(module=module_level1) %>%
    left_join(nodes.ID,by="node_id")
  
  
  # modules.aux <- tibble(year = my.year,
  #                       plot = my.plot,
  #                       guild = my.guild,
  #                       species = modules$species,
  #                       module = modules$module)
  # 
  # module_members <- bind_rows(module_members,modules.aux)
  
  # If we dont want to use bits, I guess that we could translate the previous partition
  # to other units by uning the alternative definition of modularity
  # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).
  
  # Note 2: infomap do not optimise this generalized modularity function. Thus,
  # I dont expect high values.
  
  # According to results for plot 1, year 2019 and guild == "floral visitors",
  # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
  # For each one of those nodes, we will add a new module. We denote such partition
  # as "corrected partition"
  
  # CORRECTED partition from infomap
  module_max <- max(modules$module,na.rm = T)
  
  for(i in 1:nrow(modules)){
    
    if(is.na(modules$module[i])){
      module_max <- module_max + 1
      modules$module[i] <- module_max
    }
  }
  
  # linkrank modularity -----------------------------------------------------
  # Note 1: I will use non-negative weights to be consistent with
  # the inputs that we used to feed the infomap algorithm.
  
  g  <- graph.adjacency(abs(A),weighted=TRUE)
  
  # this returns a "mask", a matrix of the same dimensions as my.matrix
  # with zeros if an interaction is forbidden, 1 if it is allowed.
  # allowed.interactions.mat <- allowed.interactions(my.matrix.null, year.plot.matrix,
  #                                                  my.guild,
  #                                                  sp.names, my.year)
  
  if(is.null(allowed.interactions.mat)){
    all.int <- A
    all.int[] <- 1
  }else{
    all.int <- allowed.interactions.mat
  }
  
  linkrank_modularity <- linkRankModularity(g,
                                            partition=modules$module,
                                            allowed.interactions.mat = all.int,
                                            damping = damping,
                                            pr.algo = pr.algo)
  return(linkrank_modularity)
}

# unweighted version
infomap.modularity.edge.list <- function(el,
                                      allowed.interactions.mat = NULL,
                                      damping = .85, 
                                      pr.algo = "prpack",
                                      temp.dir = NULL){
  g  <- graph_from_data_frame(el[,c("node_from","node_to")],directed=TRUE)
  
  # To run infomap it seems that weights should be non-negative
  # TODO adapt to weighted and non-weighted data

  my.edge.list <- el[c("node_from","node_to")]
  names(my.edge.list) <- c("from","to")
  my.edge.list$weight <- 1
  nodes <- c(el$node_from,el$node_to) %>% unique()
  
  nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
  
  # Preparing edge.list and nodes.ID to run infomap
  
  my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
    left_join(nodes.ID,by="species") %>% dplyr::select(-species) %>%
    rename(from=node_id,species=to) %>%
    left_join(nodes.ID,by="species") %>%
    dplyr::select(-species) %>% rename(to=node_id) %>% 
    dplyr::select(from,to,weight)
  
  nodes.ID2 <- nodes.ID %>% rename(node_name=species)
  
  infomap_mono <- create_monolayer_object(x = my.edge.list,
                                          directed = T,
                                          bipartite = F,
                                          node_metadata = nodes.ID2[,c(2,1)])
  
  # infomap_mono2 <- infomap_mono
  # infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
  
  # run Infomap
  if(!is.null(temp.dir)){
    # temp_dir <- paste("R/temp/d",i.id,sep="")
    # temp_dir <- "R/temp/d1"
    dir.create(temp.dir)
    modules_relax_rate <- run_infomap_monolayer3(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 temp_dir = temp.dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
    unlink(temp.dir, recursive = TRUE)
  }else{
    modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 # temp_dir = temp_dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
  }
  
  
  my.modularity <- modules_relax_rate[["L"]] # modularity in bits
  
  # Extract module information
  modules <- modules_relax_rate$modules %>%
    dplyr::select(node_id,module_level1) %>%
    rename(module=module_level1) %>%
    left_join(nodes.ID,by="node_id")
  
  
  # modules.aux <- tibble(year = my.year,
  #                       plot = my.plot,
  #                       guild = my.guild,
  #                       species = modules$species,
  #                       module = modules$module)
  # 
  # module_members <- bind_rows(module_members,modules.aux)
  
  # If we dont want to use bits, I guess that we could translate the previous partition
  # to other units by uning the alternative definition of modularity
  # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).
  
  # Note 2: infomap do not optimise this generalized modularity function. Thus,
  # I dont expect high values.
  
  # According to results for plot 1, year 2019 and guild == "floral visitors",
  # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
  # For each one of those nodes, we will add a new module. We denote such partition
  # as "corrected partition"
  
  # CORRECTED partition from infomap
  module_max <- max(modules$module,na.rm = T)
  
  for(i in 1:nrow(modules)){
    
    if(is.na(modules$module[i])){
      module_max <- module_max + 1
      modules$module[i] <- module_max
    }
  }
  
  # linkrank modularity -----------------------------------------------------
  # Note 1: I will use non-negative weights to be consistent with
  # the inputs that we used to feed the infomap algorithm.
  
  # g  <- graph.adjacency(abs(A),weighted=TRUE)
  
  # this returns a "mask", a matrix of the same dimensions as my.matrix
  # with zeros if an interaction is forbidden, 1 if it is allowed.
  # allowed.interactions.mat <- allowed.interactions(my.matrix.null, year.plot.matrix,
  #                                                  my.guild,
  #                                                  sp.names, my.year)
  
  if(is.null(allowed.interactions.mat)){
    all.int <- as.matrix(as_adjacency_matrix(g,sparse = T))
    all.int[] <- 1
  }else{
    all.int <- allowed.interactions.mat
  }
  
  linkrank_modularity <- linkRankModularity(g,
                                            partition=modules$module,
                                            allowed.interactions.mat = all.int,
                                            damping = damping,
                                            pr.algo = pr.algo)
  return(linkrank_modularity)
}


