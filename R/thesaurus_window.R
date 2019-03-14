thesaurus_othero <- function(organism, uniprot, sequence, mainpath) {
  
  cat(paste("\n\nOrganism", organism, sep = " : "))
  cat(paste("\nUniprot file", uniprot, sep = " : "))
  
  cat("\n\n>THESAURUS CONSTRUCTION ... ")
  
  # Creation du dossier correspondant a l'organisme
  organism <- gsub(" ", "-", organism)
  organisme <- gsub("-", "+", organism)
  organisme <- tolower(organisme)
  
  organism.path <- paste(mainpath, organism, sep = '/')
  dir.create(organism.path, showWarnings = FALSE)
  
  sequences <- paste(sequence, "sequences", sep = '-')
  nomsortie <- paste("Thesaurus", "_", organism, "_", sequences, ".txt", sep = '')
  sortie <- paste(organism.path, nomsortie, sep = '/')
  
  setwd(organism.path)
  
  # Construction du thesaurus
  path <- paste(system.file(package = "appinetwork"), "python/thesaurusPy.py", sep = "/")
  file.copy(path, mainpath)
  path.copy <- paste(mainpath, "thesaurusPy.py", sep = "/")
  
  setwd(mainpath)
  setwd("inst")
  command <- paste("python", "thesaurusPy.py", uniprot, sortie, organisme, sep = " ")
  system(command)
  setwd(mainpath)
  
  file.remove(path.copy)
  
  cat("OK\n>Thesaurus construction is done.")
  cat(paste("\nThesaurus file (output)", sortie, sep = " : "))
  cat(paste("\n\n"))
  
  setwd(mainpath)
  
  visible(mainpanel) <- T
  
}


thesaurus_window <- function(f_pos, mainpanel, mainpath) {

	uniprot <- c()
	return.parameter <- c()

	panel_para <- gwindow("Construct an ID correspondence file", parent = f_pos, visible = T)
	pp <- gvbox(container = panel_para)
	pp$set_borderwidth(10L)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Selection of parametres (organism and type of Uniprot file (all or proteome))
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)
	gradio(c("all", "proteome"), selected = 1, horizontal = TRUE, label = "Organism sequences in uniprot file", container = flyt)

	# Select the Uniprot file
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Uniprot file", container = chdb, handler = function(...) {
		uniprot <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(uniprot) == T) {
			gmessage('Selected file is null', icon = 'error')
		}
		if (is.null(uniprot) == F) {
			bouton1$set_value(paste(length(uniprot), 'uniprot file selected'))
			#cat(paste('\n>File selected : ', uniprot, sep = ''))
		}
	})

	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Construct file", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# Memorization of the selected parameters (organism and type of Uniprot file)
		organism <- as.character(return.parameter[1])
		sequence <- as.character(return.parameter[2])


		# Execution of thesaurusPy.py script when all the parameters are provided
		if (is.null(uniprot) == F) {
			mainpath <- getwd()

			if (organism == "Other") {
				panelorganism <- gwindow("Organism description", parent = f_pos, visible = F, expand = T)

				pc <- ggroup(container = panelorganism, horizontal = F, use.scrollwindow = T)
				pcsb <- ggroup(container = pc, horizontal = F, use.scrollwindow = F)
				lg <- gvbox(container = pcsb)
				fllg <- gformlayout(container = lg)

				organismName <- gedit(initial.msg = 'Organism Name', label = "NAME", container = fllg)

				visible(panelorganism)  <- T

				bpc <- ggroup(container = pc); addSpring(bpc)

				# 1 : Other organism
				bouton1 <- gbutton("OK", handler = function(h,...) {
					# Group the different modifications to do
					org.name <- cbind(names(svalue(fllg)), svalue(fllg))
					org.name <- data.frame(org.name, row.names = NULL)

					org <- cbind(org.name)
					colnames(org) <-  c('Organism_name')

					organism <- as.character(org[1,2])

					thesaurus_othero(organism, uniprot, sequence, mainpath)

					dispose(bpc)
				}, container = bpc)

			}

			else {

				cat(paste("\n\nOrganism", organism, sep = " : "))
				cat(paste("\nUniprot file", uniprot, sep = " : "))

				cat("\n\n>BUILDING OF THE ID CORRESPONDENCE FILE  ... ")

				# Create a directory having the name of the organism
				organism <- gsub(" ", "-", organism)
				organisme <- gsub("-", "+", organism)
				organisme <- tolower(organisme)

				organism.path <- paste(mainpath, organism, sep = '/')
				dir.create(organism.path, showWarnings = FALSE)

				sequences <- paste(sequence, "sequences", sep = '-')
				nomsortie <- paste("Thesaurus", "_", organism, "_", sequences, ".txt", sep = '')
				sortie <- paste(organism.path, nomsortie, sep = '/')

				setwd(organism.path)

				# execute thesaurusPy.py (python script)  to create the ID correspondence file
				path <- paste(system.file(package = "appinetwork"), "thesaurusPy.py", sep = "/")
				file.copy(path, mainpath)
				path.copy <- paste(mainpath, "thesaurusPy.py", sep = "/")

				setwd(mainpath)
				cat("\n CAT 1 :",mainpath)
				command <- paste("python", "thesaurusPy.py", uniprot, sortie, organisme, sep = " ")
        cat("\n CAT 2 :",command)
				system(command)

				file.remove(path.copy)

				cat("OK\n>the ID correspondence file is done.")
				cat(paste("\nID correspondence file (output)", sortie, sep = " : "))
				cat(paste("\n\n"))

				setwd(mainpath)

				visible(mainpanel) <<- T

			}

			dispose(panel_para)
			dispose(ppb)

		}

		# Test the presence of all elements necessary to execute thesaurusPy.py
		else if (is.null(uniprot) == T) {
			gmessage('File selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else {
			gmessage('Error : Unable to construct the ID correspondence file ', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}

	}, container=ppb)

	# Return
	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(mainpanel) <<- T
	}, container = ppb)

	#visible(panel_para) <- T

}

