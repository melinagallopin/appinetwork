#' Cluster PPI network
#'
#' PPI network clustering to identify proteins of a biological process
#'
#' @param network File or data frame containing the network (full or not)
#' @param method Clustering method (default: "TFit")
#' @return A list of clusters (using labels of vertices)
#'
#' @examples
#' filename <- system.file("extdata", "network.example", package = "appinetwork")
#' cl <- clusterNetwork(filename)
#' cl
#'
#' @export
clusterNetwork <- function(network, method="TFit") #TODO: any other methods ?
{
  # Obtain input data
  if (is.character(network))
    network <- load_network(network)
  if (ncol(network) > 2)
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
    # NOTE: as.character(...) because if labels are integers
    # R casts edges[i,k] into an int (and reach wrong cell)
    integerEdge1 <- integerLabelsList[[ as.character(edges[i,1]) ]]
    integerEdge2 <- integerLabelsList[[ as.character(edges[i,2]) ]]
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

#' Cluster (real) dataset
#'
#' Clustering of data based on a distances matrix,
#' used to first build a neighborhood graph
#'
#' @param D n x n distances or dissimilarity matrix (between data rows)
#' @param k Number of neighbors for each vertex. Default: 0 for adaptive
#'
#' @examples
#' cl <- clustData(iris[,1:4], 12)
#' plot(iris[,1], iris[,3], col=cl)
#'
#' @export
clusterData <- function(D, k)
{
  n <- nrow(D)
  neighbs <- list()
  if (k <= 0)
  {
    # TODO: Which strategy ?
    # ...
  }
  else
  {
    for (i in 1:n)
      neighbs[[i]] <- sort(D[i,], index.return=T)$ix[2:(k+1)]
  }
  A <- matrix(nrow=0, ncol=2)
  for (i in 1:n)
    A <- rbind(A, cbind(i, neighbs[[i]]))
  cl <- clustering(A)
  K <- length(cl)
  colors <- rep(0,n)
  for (i in 1:K)
    colors[as.integer(cl[[i]])] <- i
  colors
}
