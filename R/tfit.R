#' tfit
#'
#' Description de tfit
#'
#' @param X Fichier network d'entrée
#' @return Partition du graphe + labels
#'
#' @export
tfit <- function(X, out)
{
  edges <- read.table(X, headers=TRUE)
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
  partition <- .Call("tfit", adjacencyMatrix, PACKAGE = "appinetwork")
  list("partition"=partition, "labels"=labels)
}
