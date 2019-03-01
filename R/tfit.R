#' tfit
#'
#' Description de tfit
#'
#' @param X un nom de fichier ?
#'
#' @return Bonne question ? Il y a un peu de boulot à ce niveau puisque pour l'instant ça ne fait qu'afficher les résultats + stocker dans un fichier.
#'   Option simple : tout faire via des fichiers, puis lire le fichier depuis R
#'   Option moins simple : récupérer une sortie structurée sans passer par un fichier .clas
#'
#' @export
tfit <- function(X, out)
{
  out = .Call("tfit", X, out, PACKAGE = "appinetwork")
}
