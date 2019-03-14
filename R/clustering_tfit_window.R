clustering_tfit <- function(organism, network) {
  
  mainpath <- getwd()
  
  cat("\n\n>CLUSTERING : TFIT\n")
  
  # Lecture des donnees
  tab <- load_network(network)
  network <- as.matrix(tab)
  int_number <- dim(network)[1]
  
  # Creation des dossiers de sauvegarde
  setwd(mainpath)
  organism <- gsub(" ", "-", organism)
  organism.path <- paste(mainpath, organism, sep = '/')
  dir.create(organism.path, showWarnings = FALSE)  
  analyze.path <- paste(organism.path, "Analyzes", sep = '/')
  dir.create(analyze.path, showWarnings = FALSE)  
  saveTFit.path <- paste(analyze.path, '/TFit', sep = '')
  dir.create(saveTFit.path, showWarnings = FALSE)
  
  # Generation d'un fichier contenant les paires de proteines en interaction dans le reseau analyse
  write.table(int_number, file = paste(saveTFit.path, '/graph.gr', sep = ''), row.names = F, quote = F, col.names = F, sep = " ")
  write.table(network[,c(1,3)], file = paste(saveTFit.path, '/graph.gr', sep = ''), row.names = F, quote = F, col.names = F, sep = " ", append = TRUE)
  File.clas <- paste(saveTFit.path, '/class.clas', sep = '')
  File.gr <-  paste(saveTFit.path, '/graph.gr', sep = '')
  setwd(mainpath)
  
  # run TFitW
  cat('\n>Running TFit...')
  
  # Appel de la fonction tfit
  #tfitc("graph\nY\nall") ##################################################
  
  cat("OK\n\n>Clustering TFit is done.\n")
  cat(paste('\n>Results are stocked in : ', saveTFit.path))
  setwd(mainpath)
  
}

clustering_tfit_window <- function(f_pos, mainpanel, panb, mainpath) {

	network <- c()
	return.parameter <- c()

	panel_para <- gwindow("TFit : help(clustering_tfit)", parent = f_pos, visible = T)
	pp <- gvbox(container = panel_para)
	pp$set_borderwidth(10L)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Selection des parametres
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)

	# Selection du fichier
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Network", container = chdb, handler = function(...) {
		network <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(network) == T) {
			gmessage('Selected file is null', icon = 'error')
		}
		if (is.null(network) == F) {
			bouton1$set_value(paste(length(network), 'network file selected'))
			cat(paste('\n>File selected : ', network, sep = ''))
		}
	})

	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Clustering", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# Memorisation des parametres de recherche selectionnes
		organism <- as.character(return.parameter[1])


		# Execution de la fonction lorsque tous les parametres sont fournis
		if (is.null(network) == F) {
			######################################################################################
			########################### Execution de la fonction #################################
  
			mainpath <- getwd()

			cat("\n\n>CLUSTERING : TFIT\n")

			if (organism == "Other") {
			  panelorganism <- gwindow("Organism description", parent = f_pos, visible = F, expand = T)
			  
			  pc <- ggroup(container = panelorganism, horizontal = F, use.scrollwindow = T)
			  pcsb <- ggroup(container = pc, horizontal = F, use.scrollwindow = F)
			  lg <- gvbox(container = pcsb)
			  fllg <- gformlayout(container = lg)
			  
			  organismName <- gedit(initial.msg = 'Organism Name', label = "NAME", container = fllg)
			  
			  visible(panelorganism)  <- T
			  
			  bpc <- ggroup(container = pc); addSpring(bpc)
			  
			  # 1 : Autre organism
			  bouton1 <- gbutton("OK", handler = function(h,...) {
			    # Rassemblement des modifications a apporter
			    org.name <- cbind(names(svalue(fllg)), svalue(fllg))
			    org.name <- data.frame(org.name, row.names = NULL)
			    
			    org <- cbind(org.name)
			    colnames(org) <-  c('Organism_name')
			    
			    organism <- as.character(org[1,2])
			    setwd(mainpath)
			    visible(mainpanel) <<- T
			    clustering_tfit(organism, network)
					dispose(bpc)		
				}, container = bpc)

			}

			# Parametres : organism,network
			else{
			  setwd(mainpath)
			  visible(mainpanel) <<- T
			  clustering_tfit(organism, network)
			}
			dispose(panel_para)
			dispose(ppb)

		}

		# Tests de la presence de tous les elements necessaires a l'execution de la fonction
		else if (is.null(network) == T) {
			gmessage('Network file selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else {
			gmessage('Error : Unable to start clustering', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}

	}, container=ppb)

	# Retour
	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(panb) <<- T
	}, container = ppb)
	
	visible(panel_para) <- T

}



