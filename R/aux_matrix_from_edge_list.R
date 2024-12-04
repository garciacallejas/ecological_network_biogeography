
#' Generate adjacency matrix from edge list
#'
#' given an edge list with columns "node_from" and "node_to", returns
#' an adjacency matrix, either binary or weighted. The matrix can be 
#' expanded as a square matrix, with all node names in rows and columns,
#' so it can be treated as a unipartite matrix.
#'
#' @param edge.list dataframe with columns "node_from", "node_to", 
#' and optionally "link_value"
#' @param weighted if TRUE, the adjacency matrix will be weighted
#' @param unipartite if TRUE, it will return a square matrix with all names
#' from "node_to" and "node_from" in rows and columns
#' @param directed only for unipartite matrices. If TRUE, the returning matrix
#' will be directed, i.e. potentially asymmetric. If FALSE, the returning matrix
#' will be undirected, i.e. symmetric.
#' @param diagonal only for unipartite matrices. If TRUE, diagonal values will be 
#' set to 1.
#'
#' @return matrix
#' @export
#'
#' @examples
matrix_from_edge_list <- function(edge.list, weighted = FALSE, 
                                  unipartite = FALSE, directed = FALSE, diagonal = FALSE){
  
  sp.from <- unique(edge.list$node_from)
  sp.to <- unique(edge.list$node_to)
  
  if(unipartite){
    all.sp <- unique(c(sp.from,sp.to))
    A <- matrix(0,nrow = length(all.sp),ncol = length(all.sp), dimnames = list(all.sp,all.sp))
  }else{
    A <- matrix(0,nrow = length(sp.to),ncol = length(sp.from),dimnames = list(sp.to,sp.from))
  }
  
  if(weighted){
    for(i in 1:nrow(edge.list)){
      A[edge.list[i,"node_to"],edge.list[i,"node_from"]] <- edge.list$link_value[i]
      
      if(unipartite & !directed){
        A[edge.list[i,"node_from"],edge.list[i,"node_to"]] <- edge.list$link_value[i]
      }
      
    }
  }else{
    for(i in 1:nrow(edge.list)){
      A[edge.list$node_to[i],edge.list$node_from[i]] <- 1
      
      if(unipartite & !directed){
        A[edge.list$node_from[i],edge.list$node_to[i]] <- 1
      }
    }
  }
  
  if(unipartite & diagonal){
    diag(A) <- 1
  }
  
  # A[edge.list[,c("node_to","node_from")]] <- edge.list[,"link_value"]
  return(A)
  
}
