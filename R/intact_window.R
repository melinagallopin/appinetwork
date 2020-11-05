# intact_othero <- function(organism, organismID, file, pathracine) {
#   
#   cat("\n\n>FORMATING INTACT DATABASE")
#   
#   intactFileName <- file
#   
#   cat("\n>Loading database ... ")
#   
#   intact <- read.delim(file = intactFileName, header = T, sep = "\t")
#   
#   if (dim(intact)[2] != 15) {
#     cat(paste("\nIncorrect data dimension, read userguide for more informations"))
#     stop()
#   }
#   
#   intact <- intact[,-c(3,4,14)]
#   intact <- intact[grep("uniprotkb:", intact[,1]),]
#   intact <- intact[grep("uniprotkb:", intact[,2]),]
#   
#   if (organism == "Saccharomyces cerevisiae") {
#     organismID <- "559292"
#   }
#   if (organism == "Homo-sapiens") {
#     organismID <- "9606"
#   }
#   if (organism == "Caenorhabditis elegans") {
#     organismID <- "6239"
#   }
#   if (organism == "Drosophila melanogaster") {
#     organismID <- "7227"
#   }
#   if (organism == "Escherichia coli") {
#     organismID <- "562"
#   }
#   if (organism == "Mus musculus") {
#     organismID <- "10090"
#   }
#   if (organism == "Rattus norvegicus") {
#     organismID <- "10116"
#   }
#   
#   rtaxon <- paste("taxid:", organismID, "(", organism, ")", sep = '')
#   taxon <- organismID
#   
#   intact <- intact[grep(taxon, intact[,8]),]
#   intact <- intact[grep(taxon, intact[,9]),]
#   
#   intact[,8] <- rtaxon
#   intact[,9] <- rtaxon
#   
#   if (length(grep("colocalization", intact[,10])) > 0) {
#     intact <- intact[-grep("colocalization", intact[,10]),]
#   }
#   irefindex <- intact
#   irefindex[,10] <- rep("physical", dim(irefindex)[1])
#   irefindex[,1] <- apply(as.matrix(apply(as.matrix(irefindex[,1]), 1, strsplit, split = ":")), 1, unlist)[2,]
#   irefindex[,2] <- apply(as.matrix(apply(as.matrix(irefindex[,2]), 1, strsplit, split = ":")), 1, unlist)[2,]
#   if (length(grep("-", irefindex[,1])) > 0) {
#     irefindex[grep("-", irefindex[,1]),1] <- apply(as.matrix(apply(as.matrix(irefindex[grep("-", irefindex[,1]),1]), 1, strsplit, split = "-")), 1, unlist)[1,]
#   }
#   if (length(grep("-", irefindex[,2])) > 0) {
#     irefindex[grep("-", irefindex[,2]),2] <- apply(as.matrix(apply(as.matrix(irefindex[grep("-", irefindex[,2]),2]), 1, strsplit, split = "-")), 1, unlist)[1,]
#   }
#   
#   tmpA1 <- rep("NA", length(irefindex[,3]))
#   tmpB1 <- rep("NA", length(irefindex[,3]))
#   tmpA2 <- rep("NA", length(irefindex[,3]))
#   tmpB2 <- rep("NA", length(irefindex[,3]))
#   
#   tmpA <- apply(as.matrix(apply(as.matrix(irefindex[,3]), 1, strsplit, split = "\\|")), 1, unlist)
#   tmpB <- apply(as.matrix(apply(as.matrix(irefindex[,4]), 1, strsplit, split = "\\|")), 1, unlist)
#   
#   cat(paste("OK\n>Parsing interactions\n"))
#   prsbar <- txtProgressBar(min = 1, max = length(irefindex[,3]), style = 3)
#   for (i in 1:length(irefindex[,3])) {
#     tmpA1[i] <- tmpA[[i]][1]
#     if (length(grep("\\(gene name\\)", tmpA[[i]])) > 0) {
#       tmpA2[i] <- tmpA[[i]][grep("\\(gene name\\)", tmpA[[i]])]
#     }
#     else {
#       if (length(grep("\\(locus name\\)", tmpA[[i]])) > 0) {
#         tmpA2[i] <- tmpA[[i]][grep("\\(locus name\\)", tmpA[[i]])]
#       }
#       else {
#         if (length(grep("\\(orf name\\)", tmpA[[i]])) > 0) {
#           tmpA2[i] <- tmpA[[i]][grep("\\(orf name\\)", tmpA[[i]])]
#         }
#         else {
#           if (length(grep("\\(display_long\\)", tmpA[[i]])) > 0) {
#             tmpA2[i] <- tmpA[[i]][grep("\\(display_long\\)", tmpA[[i]])]
#           }
#         }
#       }
#     }
#     tmpB1[i] <- tmpB[[i]][1]
#     if (length(grep("\\(gene name\\)", tmpB[[i]])) > 0) {
#       tmpB2[i] <- tmpB[[i]][grep("\\(gene name\\)", tmpB[[i]])]
#     }
#     else {
#       if (length(grep("\\(locus name\\)", tmpB[[i]])) > 0) {
#         tmpB2[i] <- tmpB[[i]][grep("\\(locus name\\)", tmpB[[i]])]
#       }
#       else {
#         if (length(grep("\\(orf name\\)", tmpB[[i]])) > 0) {
#           tmpB2[i] <- tmpB[[i]][grep("\\(orf name\\)", tmpB[[i]])]
#         }
#         else {
#           if (length(grep("\\(display_long\\)", tmpB[[i]])) > 0) {
#             tmpB2[i] <- tmpB[[i]][grep("\\((display_long)\\)", tmpB[[i]])]
#           }
#         }
#       }
#     }
#     setTxtProgressBar(prsbar, i)
#   }
#   
#   GeneNameA <- sub("\\(gene name\\)", '', sub("uniprotkb:", '', tmpA2))
#   GeneNameA <- sub("\\(locus name\\)", '', GeneNameA)
#   GeneNameA <- sub("\\(orf name\\)", '', GeneNameA)
#   GeneNameA <- sub("\\(gene name synonym\\)", '', GeneNameA)
#   GeneNameA <- sub("\\(display_long\\)", '', GeneNameA)
#   GeneNameA[grep("CELE_", GeneNameA)] <- sub("CELE_", "", GeneNameA[grep("CELE_", GeneNameA)])
#   GeneNameA[grep("psi-mi:", GeneNameA)] <- sub("_caeel", "", sub("psi-mi:", "", GeneNameA[grep("psi-mi:", GeneNameA)]))
#   
#   GeneNameB <- sub("\\(gene name\\)", '', sub("uniprotkb:", '', tmpB2))
#   GeneNameB <- sub("\\(locus name\\)", '', GeneNameB)
#   GeneNameB <- sub("\\(orf name\\)", '', GeneNameB)
#   GeneNameB <- sub("\\(gene name synonym\\)", '', GeneNameB)
#   GeneNameB <- sub("\\(display_long\\)", '', GeneNameB)
#   GeneNameB[grep("CELE_", GeneNameB)] <- sub("CELE_", "", GeneNameB[grep("CELE_", GeneNameB)])
#   GeneNameB[grep("psi-mi:", GeneNameB)] <- sub("_caeel", "", sub("psi-mi:", "", GeneNameB[grep("psi-mi:", GeneNameB)]))
#   
#   irefindex[,3] <- sub("\\(display_long\\)", '', sub("psi-mi:", '', tmpA1))
#   irefindex[,4] <- sub("\\(display_long\\)", '', sub("psi-mi:", '', tmpB1))
#   
#   numParticipants <- rep(2, dim(irefindex)[1])
#   
#   sourcedb <- rep("Intact", dim(irefindex)[1])
#   irefindex[,11] <- sourcedb
#   
#   cat("OK\n>Saving database(s) ... ")
#   
#   setwd(pathracine)
#   
#   organism <- gsub(" ", "-", organism)
#   organism.path <- paste(pathracine, organism, sep = '/')
#   dir.create(organism.path, showWarnings = FALSE)
#   database.path <- paste(organism.path, "Databases", sep = '/')
#   dir.create(database.path, showWarnings = FALSE)
#   setwd(database.path)
#   
#   irefindex <- irefindex[,c(1:12)]
#   
#   irefindex <- cbind(irefindex, numParticipants, GeneNameA, GeneNameB)
#   colnames(irefindex) <- c("uidA", "uidB", "aliasA", "aliasB", "method", "author", "pmids", "taxa", "taxb", "interactionType", "sourcedb", "confidence", "numParticipants", "GeneNameA", "GeneNameB")
#   outputfilename <- paste(organism, "_intact.txt", sep = "")
#   write.table(irefindex, fi le = outputfilename, row.names = F, col.names = T, quote = F, sep = "\t")
#   
#   cat(paste("OK\n\n>Formating intact database is done.\n\nDatabase file is saved in", database.path, sep = " : "))
#   cat(paste("\n\n"))
#   setwd(pathracine)
#   
#   visible(mainpanel) <- T
#   
# }

format_intact = function(organism, organismID, file, pathracine, other) {
  
  
  cat("\n\n>FORMATING INTACT DATABASE")
  
  intactFileName <- file
  
  cat("\n>Loading database ... ")
  
  intact <- read.delim(file = intactFileName, header = T, sep = "\t")
  
  intact <- intact[,-c(3,4,14)]
  intact <- intact[grep("uniprotkb:", intact[,1]),]
  intact <- intact[grep("uniprotkb:", intact[,2]),]
  
  
  if(other != "Other") {
    # Si l'organisme est dans les choix prédéfinis :
    list.organism<-cbind(c("Saccharomyces cerevisiae","Homo sapiens","Caenorhabditis elegans","Drosophila melanogaster","Escherichia coli","Mus musculus","Rattus norvegicus","Arabidopsis thaliana"),c("559292","9606","6239","7227","562","10090","10116","3702"))
    organismID<-list.organism[list.organism[,1]==organism,2]
  }
  
  rtaxon <- paste("taxid:", organismID, "(", organism, ")", sep = '')
  taxon <- organismID
  
  intact <- intact[grep(taxon, intact[,8]),]
  intact <- intact[grep(taxon, intact[,9]),]
  
  intact[,8] <- rtaxon
  intact[,9] <- rtaxon
  
  if (length(grep("colocalization", intact[,10])) > 0) {
    intact <- intact[-grep("colocalization", intact[,10]),]
  }
  irefindex <- intact
  irefindex[,10] <- rep("physical", dim(irefindex)[1])
  irefindex[,1] <- apply(as.matrix(apply(as.matrix(irefindex[,1]), 1, strsplit, split = ":")), 1, unlist)[2,]
  irefindex[,2] <- apply(as.matrix(apply(as.matrix(irefindex[,2]), 1, strsplit, split = ":")), 1, unlist)[2,]
  if (length(grep("-", irefindex[,1])) > 0) {
    irefindex[grep("-", irefindex[,1]),1] <- apply(as.matrix(apply(as.matrix(irefindex[grep("-", irefindex[,1]),1]), 1, strsplit, split = "-")), 1, unlist)[1,]
  }
  if (length(grep("-", irefindex[,2])) > 0) {
    irefindex[grep("-", irefindex[,2]),2] <- apply(as.matrix(apply(as.matrix(irefindex[grep("-", irefindex[,2]),2]), 1, strsplit, split = "-")), 1, unlist)[1,]
  }
  
  tmpA1 <- rep("NA", length(irefindex[,3]))
  tmpB1 <- rep("NA", length(irefindex[,3]))
  tmpA2 <- rep("NA", length(irefindex[,3]))
  tmpB2 <- rep("NA", length(irefindex[,3]))
  
  tmpA <- apply(as.matrix(apply(as.matrix(irefindex[,3]), 1, strsplit, split = "\\|")), 1, unlist)
  tmpB <- apply(as.matrix(apply(as.matrix(irefindex[,4]), 1, strsplit, split = "\\|")), 1, unlist)
  
  cat(paste("OK\n>Parsing interactions\n"))
  prsbar <- txtProgressBar(min = 1, max = length(irefindex[,3]), style = 3)
  for (i in 1:length(irefindex[,3])) {
    tmpA1[i] <- tmpA[[i]][1]
    if (length(grep("\\(gene name\\)", tmpA[[i]])) > 0) {
      tmpA2[i] <- tmpA[[i]][grep("\\(gene name\\)", tmpA[[i]])]
    }
    else {
      if (length(grep("\\(locus name\\)", tmpA[[i]])) > 0) {
        tmpA2[i] <- tmpA[[i]][grep("\\(locus name\\)", tmpA[[i]])]
      }
      else {
        if (length(grep("\\(orf name\\)", tmpA[[i]])) > 0) {
          tmpA2[i] <- tmpA[[i]][grep("\\(orf name\\)", tmpA[[i]])]
        }
        else {
          if (length(grep("\\(display_long\\)", tmpA[[i]])) > 0) {
            tmpA2[i] <- tmpA[[i]][grep("\\(display_long\\)", tmpA[[i]])]
          }
        }
      }
    }
    tmpB1[i] <- tmpB[[i]][1]
    if (length(grep("\\(gene name\\)", tmpB[[i]])) > 0) {
      tmpB2[i] <- tmpB[[i]][grep("\\(gene name\\)", tmpB[[i]])]
    }
    else {
      if (length(grep("\\(locus name\\)", tmpB[[i]])) > 0) {
        tmpB2[i] <- tmpB[[i]][grep("\\(locus name\\)", tmpB[[i]])]
      }
      else {
        if (length(grep("\\(orf name\\)", tmpB[[i]])) > 0) {
          tmpB2[i] <- tmpB[[i]][grep("\\(orf name\\)", tmpB[[i]])]
        }
        else {
          if (length(grep("\\(display_long\\)", tmpB[[i]])) > 0) {
            tmpB2[i] <- tmpB[[i]][grep("\\((display_long)\\)", tmpB[[i]])]
          }
        }
      }
    }
    setTxtProgressBar(prsbar, i)
  }
  
  GeneNameA <- sub("\\(gene name\\)", '', sub("uniprotkb:", '', tmpA2))
  GeneNameA <- sub("\\(locus name\\)", '', GeneNameA)
  GeneNameA <- sub("\\(orf name\\)", '', GeneNameA)
  GeneNameA <- sub("\\(gene name synonym\\)", '', GeneNameA)
  GeneNameA <- sub("\\(display_long\\)", '', GeneNameA)
  GeneNameA[grep("CELE_", GeneNameA)] <- sub("CELE_", "", GeneNameA[grep("CELE_", GeneNameA)])
  GeneNameA[grep("psi-mi:", GeneNameA)] <- sub("_caeel", "", sub("psi-mi:", "", GeneNameA[grep("psi-mi:", GeneNameA)]))
  
  GeneNameB <- sub("\\(gene name\\)", '', sub("uniprotkb:", '', tmpB2))
  GeneNameB <- sub("\\(locus name\\)", '', GeneNameB)
  GeneNameB <- sub("\\(orf name\\)", '', GeneNameB)
  GeneNameB <- sub("\\(gene name synonym\\)", '', GeneNameB)
  GeneNameB <- sub("\\(display_long\\)", '', GeneNameB)
  GeneNameB[grep("CELE_", GeneNameB)] <- sub("CELE_", "", GeneNameB[grep("CELE_", GeneNameB)])
  GeneNameB[grep("psi-mi:", GeneNameB)] <- sub("_caeel", "", sub("psi-mi:", "", GeneNameB[grep("psi-mi:", GeneNameB)]))
  
  irefindex[,3] <- sub("\\(display_long\\)", '', sub("psi-mi:", '', tmpA1))
  irefindex[,4] <- sub("\\(display_long\\)", '', sub("psi-mi:", '', tmpB1))
  
  numParticipants <- rep(2, dim(irefindex)[1])
  
  sourcedb <- rep("Intact", dim(irefindex)[1])
  irefindex[,11] <- sourcedb
  
  cat("Ok\n>Saving database(s) ... ")
  
  setwd(pathracine)
  
  organism <- gsub(" ", "-", organism)
  organism.path <- paste(pathracine, organism, sep = '/')
  dir.create(organism.path, showWarnings = FALSE)
  database.path <- paste(organism.path, "Databases", sep = '/')
  dir.create(database.path, showWarnings = FALSE)
  setwd(database.path)
  
  irefindex <- irefindex[,c(1:12)]
  
  irefindex <- cbind(irefindex, numParticipants, GeneNameA, GeneNameB)
  colnames(irefindex) <- c("uidA", "uidB", "aliasA", "aliasB", "method", "author", "pmids", "taxa", "taxb", "interactionType", "sourcedb", "confidence", "numParticipants", "GeneNameA", "GeneNameB")
  outputfilename <- paste(organism, "_intact.txt", sep = "")
  write.table(irefindex, file = outputfilename, row.names = F, col.names = T, quote = F, sep = "\t")
  
  cat(paste("OK\n\n>Formating intact database is done.\n\nDatabase file is saved in", database.path, sep = " : "))
  cat(paste("\n\n"))
  setwd(pathracine)
  
  visible(mainpanel) <<- T
}

intact_window <- function(f_pos, mainpanel, pana, mainpath) {

	file <- c()
	return.parameter <- c()

	panel_para <- gwindow("Intact database file formatting : ", parent = f_pos, visible = T)
	pp <- gvbox(container = panel_para)
	pp$set_borderwidth(10L)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Parametres selection
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae','Arabidopsis thaliana', 'Other'), label = "Organism", selected = 7, container = flyt)

	# File selection
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Intact file", container = chdb, handler = function(...) {
		file <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
		if (is.null(file) == T) {
			gmessage('Selected file is null', icon = 'error')
		}
		if (is.null(file) == F) {
			bouton1$set_value(paste(length(file), 'intact file selected'))
			cat(paste('\n>File selected : ', file, sep = ''))
		}
	})

	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Format file", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# selected parametres storage
		organism <- as.character(return.parameter[1])
		other = organism


		# When all paramaters are Ok
		if (is.null(file) == F) {
			######################################################################################
			########################### FORMATING AND PARSING #################################

			pathracine <- getwd()

			if (organism == "Other") {
				panelorganism <- gwindow("Organism description", parent = f_pos, visible = F, expand = T)
	
				pc <- ggroup(container = panelorganism, horizontal = F, use.scrollwindow = T)
				pcsb <- ggroup(container = pc, horizontal = F, use.scrollwindow = F)
				lg <- gvbox(container = pcsb)
				pg <- gvbox(container = pcsb)
				fllg <- gformlayout(container = lg)
				flpg <- gformlayout(container = pg)

				organismName <- gedit(initial.msg = 'Organism Name', label = "NAME", container = fllg)
				organismId <- gedit(initial.msg = 'Organism ID', label = "NCBI ID", container = flpg)

				visible(panelorganism)  <- T
 
				bpc <- ggroup(container = pc); addSpring(bpc)

				# 1 : Other organism
				bouton1 <- gbutton("OK", handler = function(h,...) {
					# Rassemblement des modifications a apporter
					org.name <- cbind(names(svalue(fllg)), svalue(fllg))
					org.name <- data.frame(org.name, row.names = NULL)
					org.id <- cbind(names(svalue(flpg)), svalue(flpg))
					org.id <- data.frame(org.id, row.names = NULL)

					org <- cbind(org.name, org.id)
					colnames(org) <-  c('Organism_name', 'Organism_id')

					organism <- as.character(org[1,2])
					organismID <- as.character(org[1,4])

					# intact_othero(organism, organismID, file, pathracine)
					
					format_intact(organism, organismID, file, pathracine, other)

					dispose(bpc)		
				}, container = bpc)

			}

			else {
			  
			  format_intact(organism, organismID, file, pathracine, other)
			  
			}

			dispose(panel_para)
			dispose(ppb)

		}

		# Tests all elements necessary to the function execution
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

