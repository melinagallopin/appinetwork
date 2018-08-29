irefindex_window <- function(f_pos, mainpanel, pana, mainpath) {

	file <- c()
	return.parameter <- c()

	panel_para <- gwindow("IrefIndex database file formatting", parent = f_pos, visible = T)
	pp <- gvbox(container = panel_para)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Parametres selection
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)
	gcombobox(c(0:20), label = "Maximum number of proteins in a complex", selected = 1, container = flyt)

	# File section
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Irefindex file", container = chdb, handler = function(...) {
		file <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(file) == T) {
			gmessage('Selected file is null', icon = 'error')
		}
		if (is.null(file) == F) {
			bouton1$set_value(paste(length(file), 'irefindex file selected'))
			#cat(paste('\n>File selected : ', file, sep = ''))
		}
	})

	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Formate file", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# Storage of the parameters
		organism <- as.character(return.parameter[1])
		complex <- as.character(return.parameter[2])


		# Parse and formate
		if (is.null(file) == F) {

			mainpath <- getwd()
            # Other organism
			if (organism == "Other") {
				panelorganism <- gwindow("Organism description", parent = f_pos, visible = F, expand = T)
	
				pc <- ggroup(container = panelorganism, horizontal = F, use.scrollwindow = T)
				pcsb <- ggroup(container = pc, horizontal = F, use.scrollwindow = F)
				lg <- gvbox(container = pcsb)
				fllg <- gformlayout(container = lg)

				organismName <- gedit(initial.msg = 'Organism Name', label = "NAME", container = fllg)

				visible(panelorganism)  <- T
 
				bpc <- ggroup(container = pc); addSpring(bpc)

				bouton1 <- gbutton("OK", handler = function(h,...) {
					# Rassemblement des modifications a apporter
					org.name <- cbind(names(svalue(fllg)), svalue(fllg))
					org.name <- data.frame(org.name, row.names = NULL)

					org <- cbind(org.name)
					colnames(org) <-  c('Organism_name')

					organism <- as.character(org[1,2])

					irefindex_othero(organism, complex, file, mainpath)

					dispose(bpc)		
				}, container = bpc)

			}
            # Organism proposed by the dialog window
			else {

				cat("\n>FORMATING IREFINDEX DATABASE\n")

				if (organism == "Saccharomyces cerevisiae") {
					organisme <- "Saccharomyces cerevisiae S288C"
				}
				else {
					organisme <- organism
				}

				cat(paste("\nOrganism", organisme, sep = " : "))
				cat(paste("\nIrefIndex file", file, sep = " : "))
				cat(paste("\n"))

				while (!file.exists(file)) {
					cat("\n\t-> incorrect database name ! try again\n")
					file <- readline("\t> Type the name of the IrefIndex data set file that you want to parse : ")
				}
	
				irefindexfile <- file

				setwd(mainpath)

				cat("\n>Loading database ... ")
				Base <- read.delim(file = irefindexfile, header = T, sep = "\t")
				cat("OK")

				if (dim(Base)[2] != 54) {
					cat(paste("\nIncorrect data dimension, read help(irefindex) for more informations"))
					stop()
				}

				colnames(Base)[1] <- "uidA"
				Base.f <- Base[,c(1:2,5:13,15,54)]
				Interactions <- Base.f

				seuil <- as.integer(complex)

				cat("\n>Parsing interactions ... ")

				# Add PPIs between proteins of the complexes having less proteins than the chosen Maximum number of proteins in a complex and delete lines containing complexes
				lignecomplexe <- c(grep("complex:", Interactions[,1]), grep("complex:", Interactions[,2]))
				Interactions.noComplex <- Interactions[-lignecomplexe,]
				Interactions.Complex <- Interactions[lignecomplexe,]
				Interactions.Complex <- Interactions.Complex[as.integer(Interactions.Complex[,13]) <= seuil,]
				nbcomplex <- 0
				l <- dim(Interactions.Complex)[1]
				i <- 1
				while (i < l) {
					ListProtofComplex <- unique(Interactions.Complex[Interactions.Complex[,1] == Interactions.Complex[i,1],2])
					ListProtofComplex <- as.vector(ListProtofComplex)
					ListAliasofComplex <- unique(Interactions.Complex[Interactions.Complex[,1] == Interactions.Complex[i,1],4])
					ListAliasofComplex <- as.vector(ListAliasofComplex)
					n <- length(ListProtofComplex)
					if (n > 1) {
						nbcomplex <- nbcomplex + 1
						ListPPIofComplex <- as.data.frame(matrix(NA, nrow = (n * (n - 1) / 2), ncol = dim(Interactions.Complex)[2]))
						j <- 1
						for (k in 1:(n - 1)) {
							for (l in (k + 1):n) {
								ListPPIofComplex[j,1] <- ListProtofComplex[k]
								ListPPIofComplex[j,2] <- ListProtofComplex[l]
								ListPPIofComplex[j,3] <- ListAliasofComplex[k]
								ListPPIofComplex[j,4] <- ListAliasofComplex[l]
								ListPPIofComplex[j,5] <- as.vector(Interactions.Complex[i,5])
								ListPPIofComplex[j,6] <- as.vector(Interactions.Complex[i,6])
								ListPPIofComplex[j,7] <- as.vector(Interactions.Complex[i,7])
								ListPPIofComplex[j,8] <- as.vector(Interactions.Complex[i,8])
								ListPPIofComplex[j,9] <- as.vector(Interactions.Complex[i,9])
								ListPPIofComplex[j,10] <- as.vector(Interactions.Complex[i,10])
								ListPPIofComplex[j,11] <- as.vector(Interactions.Complex[i,11])
								ListPPIofComplex[j,12] <- as.vector(Interactions.Complex[i,12])
								ListPPIofComplex[j,13] <- n
								j <- j + 1
							}
						}
		
						if (nbcomplex == 1) {
							tmp <- ListPPIofComplex
						}
						if (nbcomplex > 1) {
							tmp <- rbind(tmp, ListPPIofComplex)
						}
					}
					Interactions.Complex <- Interactions.Complex[-c(1:dim(Interactions.Complex)[1])[Interactions.Complex[,1] == Interactions.Complex[i,1]],]
					l <- dim(Interactions.Complex)[1]
				}
				colnames(tmp) <- colnames(Interactions)
				Interactions <- rbind(Interactions.noComplex, tmp)

				cat("OK")

				# All lines that don't contain "uniprotkb" and "refseq:NP_" in the columns uidA and uidB are deleted.
				cat(paste("\nChecking uniprot-IDs\n"))
				prsbar <- txtProgressBar(min = 1, max = dim(Interactions)[1], style = 3)
				SUP <- rep(TRUE, dim(Interactions)[1])
				for (i in 1:(dim(Interactions)[1])) {
					if ((length(grep("uniprotkb", Interactions[i,1], value = FALSE)) != 1) && (length(grep("refseq:NP_", Interactions[i,1], value = FALSE)) != 1)) {
						SUP[i] <- FALSE
					}
					else {
						if ((length(grep("uniprotkb", Interactions[i,2], value = FALSE)) != 1) && (length(grep("refseq:NP_", Interactions[i,2], value = FALSE)) != 1)) {
							SUP[i] <- FALSE
						}
					}
					setTxtProgressBar(prsbar, i)
				}
				Interactions <- Interactions[SUP == TRUE,]

				# All lines that don't contain the specie's name in the columns taxA and taxB are deleted.
                
				cat(paste("\nChecking taxons\n"))
				prsbar <- txtProgressBar(min = 1, max = dim(Interactions)[1], style = 3)
				SUP <- rep(FALSE, dim(Interactions)[1])
				for (i in 1:(dim(Interactions)[1])) {
					if (length(grep(organisme, Interactions[i,8], value = FALSE)) == 1 && length(grep(organisme, Interactions[i,9], value = FALSE)) == 1) {
						SUP[i] <- TRUE
					}
					else if (length(grep(organisme, Interactions[i,8], value = FALSE)) == 1 && length(grep("-", Interactions[i,9], value = FALSE)) == 1) {
						SUP[i] <- TRUE
					}
					else if (length(grep("-", Interactions[i,8], value = FALSE)) == 1 && length(grep(organisme, Interactions[i,9], value = FALSE)) == 1) {
						SUP[i] <- TRUE
					}
					setTxtProgressBar(prsbar, i)
				}
				Interactions <- Interactions[SUP == TRUE,]
				# Interactions<-Interactions[as.character(Interactions[,8])==as.character(Interactions[,9]),]

				Interactions <- as.matrix(Interactions)

				# Troncate uids
				Interactions[,1] <- apply(as.matrix(apply(as.matrix(Interactions[,1]), 1, strsplit, split = ":")), 1, unlist)[2,]
				Interactions[,2] <- apply(as.matrix(apply(as.matrix(Interactions[,2]), 1, strsplit, split = ":")), 1, unlist)[2,]

				# Give the geneNames
				GeneNameA <- Interactions[,1]
				GeneNameB <- Interactions[,2]

				cat(paste("\nSearching gene name\n"))
				prsbar <- txtProgressBar(min = 1, max = dim(Interactions)[1], style = 3)
				for (i in 1:(dim(Interactions)[1])) {
					listA <- unlist(strsplit(as.character(Interactions[i,3]), "\\|"))
					for (j in 1:length(listA)) {
						if (length(grep("locuslink:", listA[j])) == 1) {
							GeneNameA[i] <- unlist(strsplit(as.character(listA[j]), "locuslink:"))[2]
						}
						else if (length(grep("hgnc:", listA[j])) == 1) {
							GeneNameA[i] <- unlist(strsplit(as.character(listA[j]), "hgnc:"))[2]
						}
					}
					if (GeneNameA[i] == "NA") {
						GeneNameA[i] <- Interactions[i,1]
					}

					listB <- unlist(strsplit(as.character(Interactions[i,4]), "\\|"))
					for (j in 1:length(listB)) {
						if (length(grep("locuslink:", listB[j])) == 1) {
							GeneNameB[i] <- unlist(strsplit(as.character(listB[j]), "locuslink:"))[2]
						}
						else if (length(grep("hgnc:", listB[j])) == 1) {
							GeneNameB[i] <- unlist(strsplit(as.character(listB[j]), "hgnc:"))[2]
						}
					}
					if (GeneNameB[i] == "NA") {
						GeneNameB[i] <- Interactions[i,2]
					}
		
					Interactions[i,3] <- unlist(strsplit(as.character(Interactions[i,3]), "\\|"))[1]
					Interactions[i,3] <- unlist(strsplit(as.character(Interactions[i,3]), ":"))[2]
					Interactions[i,4] <- unlist(strsplit(as.character(Interactions[i,4]), "\\|"))[1]
					Interactions[i,4] <- unlist(strsplit(as.character(Interactions[i,4]), ":"))[2]
					Interactions[i,11] <- unlist(strsplit(as.character(Interactions[i,11]), "\\("))[2]
					Interactions[i,11] <- unlist(strsplit(as.character(Interactions[i,11]), "\\)"))[1]
					setTxtProgressBar(prsbar, i)
				}

				Interactions <- cbind(Interactions, GeneNameA, GeneNameB)
				Interactions[Interactions[,11] == "grid",11] <- "BioGrid"

				# Give the interactions type
				Interactions[grep("MI:0254", Interactions[,5], value = FALSE),10] <- "genetic"
				Interactions[grep("MI:0001", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0008", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0024", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0045", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0063", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0064", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0363", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0364", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0403", Interactions[,5], value = FALSE),10] <- "other"
				Interactions[grep("MI:0686", Interactions[,5], value = FALSE),10] <- "other"
		
				Interactions[(Interactions[,10] != "other") & (Interactions[,10] != "genetic"),10] <- "physical"

				cat("\n>Saving database(s) ... ")

				List.database <- unlist(unique(Interactions[,11]))
	
				organism <- gsub(" ", "-", organism)
				organism.path <- paste(mainpath, organism, sep = '/')
				dir.create(organism.path, showWarnings = FALSE)
				database.path <- paste(organism.path, "Databases", sep = '/')
				dir.create(database.path, showWarnings = FALSE)
	
				setwd(database.path)

				for(i in 1:(length(List.database))) {
					database.output.name <- paste(organism, "irefindex", List.database[i], "threshold", complex, sep = "_")
					database.output <- Interactions[Interactions[,11] == List.database[i],]
					write.table(database.output, file = paste(database.output.name, "txt", sep = "."), row.names = F, col.names = T, quote = F, sep = "\t")
				}
	
				cat(paste("OK\n>Formating irefindex database is done.\n\nDatabase files are saved in", database.path, sep = " : "))
				cat(paste("\n\n"))
		
				setwd(mainpath)

				visible(mainpanel) <<- T

			}

			dispose(panel_para)
			dispose(ppb)

		}

		# Test the presence of all the elements necessary to execute the function
		else if (is.null(file) == T) {        
			gmessage('File selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else {
			gmessage('Error : Unable to start formatting', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}

	}, container=ppb)

	# Return
	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(pana) <<- T
	}, container = ppb)
	
	visible(panel_para) <- T

}

