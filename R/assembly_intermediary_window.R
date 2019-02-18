
assembly_intermediary_window <- function(f_pos, mainpanel) {

	network <- c()
	inputlist <- c()
	return.parameter <- c()

	panel_para <- gwindow("Assembly intermediaries : help(assembly_intermediary)", parent = f_pos, width = '500', height = '200', visible = T)
	pp <- gvbox(container = panel_para)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Selection des parametres
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)
	gcombobox(c('Chi2', 'dice', 'jaccard', 'liddell', 'ms', 'zscore'), label = "Score", selected = 3, container = flyt)

	# Selection du fichier
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Network", container = chdb, handler = function(...) {
		network <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(network) == T) {
			gmessage('Selected network file is null', icon = 'error')
		}
		if (is.null(network) == F) {
			bouton1$set_value(paste(length(network), 'network file selected'))
			cat(paste('\n>Network selected : ', network, sep = ''))
		}
	})
	bouton2 <- gbutton("Input list", container = chdb, handler = function(...) {
		inputlist <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(inputlist) == T) {
			gmessage('Selected input list file is null', icon = 'error')
		}
		if (is.null(inputlist) == F) {
			bouton2$set_value(paste(length(inputlist), 'input list file selected'))
			cat(paste('\n>Input list selected : ', inputlist, sep = ''))
		}
	})

	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Modeling", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# Memorisation des parametres de recherche selectionnes
		organism <- as.character(return.parameter[1])
		score <- as.character(return.parameter[2])


		# Execution de la fonction lorsque tous les parametres sont fournis
		if (is.null(network) == F & is.null(inputlist) == F) {

			######################################################################################
			########################### Execution de la fonction #################################

			pathracine <- getwd()

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

					assembly_intermediary_othero(organism, score, network, inputlist, pathracine)

					dispose(bpc)
				}, container = bpc)

			}

			else {

				cat("\n\n>FINDING POTENTIAL ASSEMBLY INTERMEDIARIES\n")

				# Recuperation des donnees
				network2cluter <- load_network(network)
				inputListFileName <- load_inputlist(inputlist)
				inputName <- inputListFileName[,1]
				inputID<- inputListFileName[,2]
				listProt <- unique(rbind(as.matrix(network2cluter[,4]), as.matrix(network2cluter[,5])))
				unredundant.direct <- network2cluter
				
				# Making the proximity score matrix
				cat("\n\n>Organism : ", organism)
				cat("\n>SCORE : ", score)
				cat("\n\n>Making a distance matrix...")

				mat <- matrix(nrow = length(inputID), ncol = length(inputID), dimnames = list(inputID, inputID))
				N <- length(inputID)
				allProt <- union(unredundant.direct[,4], unredundant.direct[,5])
				cat("inputID[N] : ")
				for (i in 1:N) {
					listProti <- unique(c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[i]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[i]),4])))
					for (j in 1:N) {
					  if(i!=j){
					    listProtj <- unique(c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[j]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[j]),4])))
					    mat[i,j] <- proximity_score(listProti, listProtj, inputID, allProt, score)
					  }
					}
				}
        
				# Warning message if a protein has no interaction with the others
				for (i in 1:N) {
					if (colSums(mat, na.rm = TRUE)[i] == 0) {
						cat(paste("\n\t> Warning message: ", colnames(mat)[i], " has 0 interaction with others sub-units, please remove it from the Input List and try again", "\n "), sep = "")
						stop()
					}
				}

				# Transform the matrix into values between 0 and 1
				if (score == "liddell" || score == "zscore" || score == "Chi2") {
					mat <- normalize_mat(mat)
				}

				# Making the Jaccard distance matrix
				mat.JaccardDistance <- 1 - mat

				cat("OK")

				setwd(pathracine)
				organism <- gsub(" ", "-", organism)
				dir.create(organism, showWarnings = FALSE)
				setwd(organism)
				dir.create("Analyzes", showWarnings = FALSE)
				setwd("Analyzes")
				dir.create("Assembly_intermediaries", showWarnings = FALSE)
				setwd("Assembly_intermediaries")
				dir.create(score, showWarnings = FALSE)
				setwd(score)
				colnames(mat) <- inputName
				colnames(mat.JaccardDistance) <- inputName
				rownames(mat) <- inputName
				rownames(mat.JaccardDistance) <- inputName

				# Finding Sub Complex
				cat("\n>Finding sub complex...")

				# Building of the intermediary matrix
				write.table(round(mat.JaccardDistance, 3), file = "DIS.txt")

				mat.JaccardDistance.current <- mat.JaccardDistance
				e <- c()

				while (dim(mat.JaccardDistance.current)[1] > 2) {
					indices <- indices_min_jaccard(mat.JaccardDistance.current)
					namesProtMerged <- colnames(mat.JaccardDistance.current)[indices]
					part1 <- strsplit(namesProtMerged[1], split = ",")[[1]][1]
					part2 <- strsplit(namesProtMerged[2], split = ",")[[1]][1]
					d <- mat.JaccardDistance.current[indices[1],indices[2]]
					e <- construct(part1, part2, e, mat, inputID, inputName)
					mat.JaccardDistance.current <- mat.JaccardDistance.current[-indices[2],]
					mat.JaccardDistance.current <- mat.JaccardDistance.current[,-indices[2]]

					colnames(mat.JaccardDistance.current)[indices[1]] <- paste(namesProtMerged[1], namesProtMerged[2], sep = ",")
					rownames(mat.JaccardDistance.current)[indices[1]] <- paste(namesProtMerged[1], namesProtMerged[2], sep = ",")

					ListIntersectName <- paste("Proteins", paste(namesProtMerged[1], namesProtMerged[2], sep = ","), sep = "_")

					# list of the proteins interacting with the cluster proteins
					namesProtMerged <- unlist(c(strsplit(namesProtMerged[1], ","), strsplit(namesProtMerged[2], ",")))
					indiceUnion <- NA
					L <- length(namesProtMerged)
					for (j in 1:L) {
						indiceUnion <- c(indiceUnion, (1:length(inputID))[inputName == namesProtMerged[j]])
					}
					indiceUnion <- indiceUnion[-1]
					v <- c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[indiceUnion[1]]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[indiceUnion[1]]),4]))

					ListProtIntersect <- v
					M <- length(indiceUnion)
					for (j in 2:M) {
						v <- c(v, c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[indiceUnion[j]]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[indiceUnion[j]]),4])))
						ListProtIntersect <- intersect(ListProtIntersect, c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[indiceUnion[j]]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[indiceUnion[j]]),4])))
					}
					v <- unique(v)
					unionListProtMerged <- v
					ListProtIntersect <- unique(ListProtIntersect)


					ListProtNameIntersect<-c()
					for(k in 1:length(ListProtIntersect)){
					  if(length(unredundant.direct[unredundant.direct[,4]==ListProtIntersect[k],11])>0){
					    ListProtNameIntersect<-c(ListProtNameIntersect,unique(unredundant.direct[unredundant.direct[,4]==ListProtIntersect[k],11]))
					  }
					  else{
					    ListProtNameIntersect<-c(ListProtNameIntersect,unique(unredundant.direct[unredundant.direct[,5]==ListProtIntersect[k],12]))
					  }
					}
					print(ListProtNameIntersect)

					write.table(ListProtNameIntersect, file = paste(ListIntersectName, "txt", sep = "."), row.names = FALSE, col.names = FALSE)

					distProtMerged <- rep(NA, dim(mat.JaccardDistance.current)[1])

					for (j in 1:dim(mat.JaccardDistance.current)[1]) {
						if (j != indices[1]) {
							# list of the proteins interacting with the protein(s)of the j column
							namesProtMergedj <- unlist(strsplit(colnames(mat.JaccardDistance.current)[j], ","))
							indiceUnionj <- NA
							for (k in 1:length(namesProtMergedj)) {
								indiceUnionj <- c(indiceUnionj, (1:length(inputName))[inputName == namesProtMergedj[k]])
							}
							indiceUnionj <- indiceUnionj[-1]
							v <- c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[indiceUnionj[1]]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[indiceUnionj[1]]),4]))
							for (k in 2:length(indiceUnionj)) {
								v <- c(v, c(as.vector(unredundant.direct[unredundant.direct[,4] == as.vector(inputID[indiceUnionj[k]]),5]), as.vector(unredundant.direct[unredundant.direct[,5] == as.vector(inputID[indiceUnionj[k]]),4])))
							}
							v <- unique(v)
							unionListProtMergedj <- v
							distProtMerged[j] <- 1 - proximity_score(unionListProtMerged, unionListProtMergedj, inputID, allProt, score)
						}
					}
					mat.JaccardDistance.current[indices[1],] <- distProtMerged
					mat.JaccardDistance.current[,indices[1]] <- distProtMerged
					if (score == "liddell" || score == "zscore" || score == "Chi2") {
						mat.JaccardDistance.current <- 1 - normalize_mat(1 - mat.JaccardDistance.current)
					}
					write.table(round(mat.JaccardDistance.current, 3), file = "DIS.txt", append = TRUE)
				}
				indices <- indices_min_jaccard(mat.JaccardDistance.current)
				namesProtMerged <- colnames(mat.JaccardDistance.current)
				part1 <- strsplit(namesProtMerged[1], split = ",")[[1]][1]
				part2 <- strsplit(namesProtMerged[2], split = ",")[[1]][1]
				d <- mat.JaccardDistance.current[indices[1],indices[2]]
				e <- construct(part1, part2, e, mat, inputID, inputName)
				e <- paste(e, ";", sep = "")
				write(e, file = "HC.txt")
				tr <- read.tree("HC.txt")
				jpeg("Tree.jpeg")
				plot(tr, type = "cladogram", root.edge = TRUE)
				dev.off()

				setwd(pathracine)

				cat("OK\n\n>Potential assembly intermediaries modelization is done.\n")
				cat("\n\nResults are saved in : ", pathracine, "/", organism, "/Analyzes/Assembly_intermediaries\n\n")

				visible(mainpanel) <<- T

			}

			########################### Execution de la fonction #################################
			######################################################################################

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
		else if (is.null(inputlist) == T) {
			gmessage('Input list file selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else {
			gmessage('Error : Unable to start modeling', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}

	}, container = ppb)

	# Retour
	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(mainpanel) <<- T
	}, container = ppb)

	visible(panel_para) <- T

}

