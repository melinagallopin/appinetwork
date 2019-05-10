#' clustering
#'
#' PPI network clustering to identify proteins of a biological process
#'
#' @param network File or data frame containing the network (full or not)
#' @param method Clustering method (default: "TFit")
#' @return A list of clusters (using labels of vertices)
#'
#' @export
clustering <- function(network, method="TFit") #TODO: any other methods ?
{
  # Obtain input data
  if (is.character(network))
    network <- load_network(network)
  if (ncol(network > 2))
    network <- network[,c(1,3)] #extract edges information (only)
  if (!is.matrix(network))
    network <- as.matrix(network)
  edges <- network

  # Convert input into an adjacency matrix (TODO: what about sparse graphs ?)
  nbEdges <- nrow(edges)
  labels <- unique(c(as.character(edges[,1]), as.character(edges[,2])))
  nbVertices <- length(labels)
  integerLabels <- 1:nbVertices
  # NOTE: following inspired by https://stackoverflow.com/a/46251794
  # There should be an easier solution (?!)
  names(integerLabels) <- labels
  integerLabelsList <- split(unname(integerLabels), names(integerLabels))
  adjacencyMatrix <- matrix(FALSE, nbVertices, nbVertices)
  for (i in 1:nbEdges)
  {
    integerEdge1 <- integerLabelsList[[ edges[i,1] ]]
    integerEdge2 <- integerLabelsList[[ edges[i,2] ]]
    adjacencyMatrix[integerEdge1,integerEdge2] <- TRUE
    adjacencyMatrix[integerEdge2,integerEdge1] <- TRUE #symmetric graph
  }

  # Call the method in arguments:
  if (method == "TFit")
    partition <- .Call("tfit", adjacencyMatrix, PACKAGE = "appinetwork")
#  else if (method == "...")
#    partition <- ...

  res <- list()
  K <- max(partition)
  for (k in 1:K)
    res[[k]] = labels[partition==k]
  res
}
