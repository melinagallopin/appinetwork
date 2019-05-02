build_network_othero <- function(organism, db, ipl, th, method, degree, remove.sdi, remove.ul, remove.si, update, mainpath, f_pos) {
  
  cat("\n\n>BUILD NETWORK")
  
  # Parametres de la fonction
  Os <- organism
  os <- tolower(organism)
  ex.type <- tolower(method)
  inter.type <- tolower(degree)
  r.i <- tolower(remove.sdi)
  r.s.i <- tolower(remove.si)
  r.u.l <- tolower(remove.ul)
  IPL <- ipl
  DB <- db
  Th <- th
  UpDate <- tolower(update)
  
  # Lecture des fichiers input
  selected.database <- load_data(db)
  input.uniprotID <- read.table(file = ipl, header = FALSE, sep = '\t', stringsAsFactors = F)
  Thesaurus <- read.table(file = th, header = TRUE, quote = "", sep = "\t", stringsAsFactors = F)
  
  ### Tous les parametres et fichiers fournis peuvent etre utilises ###
  
  
  # Completer la liste d'ID a rechercher dans les bases : ajout des anciens ID et isoformes ID
  listsupp <- c()
  # Parcours de tout le fichier thesaurus pour trouver tous les elements de la liste input et ajouter tous les autres uniprot-ID associe a la meme proteine
  for (i in 1:length(Thesaurus[,1])) {
    if (Thesaurus[i,1] %in% input.uniprotID[,1]) {
      oldNames <- Thesaurus[i,7]
      oldNamesList <- unlist(strsplit(oldNames, ";"))
      if (is.na(Thesaurus[i,7]) == FALSE) {
        # Ajout des anciens uniprot-ID a la liste
        for (j in 1:length(oldNamesList[])) {
          listsupp <- rbind(listsupp, oldNamesList[j])
        }
      }
      Nisoformes <- Thesaurus[i,8]
      if (Nisoformes > 0) {
        # Ajout des uniprot-ID des isoformes associes au meme gene
        for (k in 1 : Nisoformes) {
          newID <- paste(Thesaurus[i,1], k, sep = "-")
          listsupp <- rbind(listsupp, paste(Thesaurus[i,1], k, sep = "-"))
        }
      }
    }
  }
  # Apres avoir recupere tous les uniprot-ID associes a ceux de la liste input, les uniprot-ID de la liste input sont ajoutes
  input <- unlist(input.uniprotID)
  for (i in 1:length(input)) {
    uniprotID <- input[i]
    listsupp <- rbind(listsupp, uniprotID)
  }
  
  input.uniprotID <- listsupp
  
  ### La liste des uniprot-ID a rechercher dans les bases est complete ###
  
  
  # Creation d'un dossier pour les fichiers resultats en sortie
  organismdir <- gsub(" ", "-", Os)
  organism.path <- paste(mainpath, organismdir, sep = '/')
  dir.create(organism.path, showWarnings = FALSE)
  network.path <- paste(organism.path, '/Network', sep = '')
  dir.create(network.path, showWarnings = FALSE)
  
  
  # Affichage des parametres selectionnes pour la recherche d'interactions
  cat(paste('\n\n>Species selected :', os))
  cat(paste('\n>Experimental type :', ex.type))
  cat(paste('\n>Interaction type :', inter.type))
  
  # Filtres des resultats sur l'organisme recherche et la methode souhaitee, dans les bases de donnees
  selected.database2 <- selected.database[(length(grep(os, selected.database[,8], ignore.case=TRUE)) > 0 && length(grep(os, selected.database[,9], ignore.case = TRUE)) > 0),]
  if (method == "genetic") {
    selected.database3 <- selected.database2[(length(grep(ex.type, selected.database2[,10], ignore.case = TRUE)) > 0 || length(grep(ex.type, selected.database2[,10], ignore.case = TRUE)) > 0),]
  }
  else {
    selected.database3 <- selected.database2[selected.database2[,10] == ex.type,]
  }
  # Verification de la presence d'interactions dans la base
  if (dim(selected.database3)[1] == 0) {
    cat("\n\nWARNING : No result found in database(s)\nProcesuss stop")
    stop()
  }
  
  # Rassemblement de toutes les databases fournies
  selected.database4 <- as.matrix(data.frame(selected.database3[,1:15], stringAsFactors = F))
  
  ### Les bases sont filtrees sur la recherche souhaitee ###
  
  
  # Recherche des interactions directes pour les proteines de l'input list
  cat('\n\n>Search for first degree interactions')
  PPI.Direct <- recup_ppi(input.uniprotID, selected.database3)
  
  # Verification de la presence d'interactions
  if (dim(PPI.Direct)[1] == 0) {
    cat("\n\nWARNING : No interaction found\nProcesuss stop")
    stop()
  }
  
  ### Les interactions avec les proteines de l'input list sont extraites des bases ###
  
  
  # S'il y a bien des interactions identifiees dans les bases
  
  # Recherches du nombre d'article par interaction
  unredundant <- pubmed_id(PPI.Direct, os,1)
  PPI.Direct <- unredundant
  
  cat('\n')
  
  # Si le reseau demande est de degre 2 on relance une recherche d'interactions
  if (inter.type == "second-degree") {
    # Recuperation des uniprot-ID de toutes les proteines du reseau de degre 1
    Interactions.matrix <- PPI.Direct
    # Recuperation des proteines d'interet et de leurs interactants directes
    col1 <- as.matrix(Interactions.matrix[,1])
    col2 <- as.matrix(Interactions.matrix[,2])
    # Rassemblement des listes pour recuperer les uniprot-ID en une seule copie
    listProt <- rbind(col1, col2)
    listProt <- unique(listProt)
    ListProt.UniprotID <- as.matrix(listProt)
    
    # Seconde recherche d'interactions, pour le degre 2
    cat('\n>Search for second degree interactions')
    PPI.Direct2 <- recup_ppi(ListProt.UniprotID, selected.database3)
    
    # Suppression des interactions de degre 2 avec les proteines qui n'ont qu'une seule interaction dans le reseau de degre 1 : si demande
    if (r.i == 'yes') {
      PPI.Indirect <- remove_ppi(ListProt.UniprotID, PPI.Direct2, input.uniprotID)
    }
    else {
      PPI.Indirect <- PPI.Direct2
    }
    
    # Recherches du nombre d'article par interaction
    unredundant <- pubmed_id(PPI.Indirect, os,2)
    PPI.Indirect<-unredundant
    
    # Rassemblement des resultats d'interactions de degre 1 et degre 2
    network <- rbind(PPI.Direct[,c(3,5,4,1,2,6:13)], PPI.Indirect[,c(3,5,4,1,2,6:13)])
    
    ### Les interactions de degre 2 sont ajoutees au reseau ###
    
  }
  
  else {
    # On ne garde que les resultats de la premiere recherche d'interactions si on souhaite un reseau de degre 1
    network <- rbind(PPI.Direct[,c(3,5,4,1,2,6:13)])
  }
  
  # Verification de la presence de resultats
  if (is.null(network) == T) {
    cat("WARNING :\nThe network is null!!!\nProcessus stoped")
    stop()
  }
  
  cat('\n')
  
  # Premier trie des interactions du reseau en fonction des parametre de construction selectionnes
  # Remove redundants
  cat('\n>Removing redundants : first sort \n')
  selected.database3 <- remove_redundants(network)
  network <<- as.matrix(selected.database3)
  # Remove self-interactant
  if (r.s.i == 'yes'){
    cat('\n>Removing self-interactants : first sort ... ')
    network <<- network[as.vector(network[,4]) != as.vector(network[,5]),]
    cat('OK')
  }
  else {
    cat('\n>Proteins which interact with itselves are kept')
  }
  
  ### Le reseau est complet et un premier trie a ete effectue ###
  
  
  # Correction automatique du reseau avec le thesaurus : on met l'uniprot-ID principal pour chaque proteine (a la place des anciens)
  cat ( '\n>Autocorrection' )
  cat ( '\n' )
  data <- as.matrix(network)
  thesaurus <- read.table(file = th, header = TRUE, quote = "", sep = "\t" , stringsAsFactors = F)
  not_founds <- c()
  # Visualisation de la progression des corrections
  pb1 <- txtProgressBar(min = 0, max = dim(data)[1], style = 3)
  
  # Parcours de toutes les interactions du reseau
  for (i in 1:dim(data)[1]) {
    # Recherche de l'uniprot ID de reference et association du bon nom de proteine et de gene
    # Colonne A
    resultatA <- search_id(data[i,4], thesaurus)
    if (length(resultatA) == 3) {
      ID <- resultatA[1]
      Proteine <- resultatA[2]
      Gene <- resultatA[3]
      data[i,4] <- ID
      data[i,1] <- Proteine
      data[i,11] <- Gene
    }
    else {
      # On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
      not_found <- data[i,c(4,1,11)]
      not_founds <- rbind(not_founds, not_found)
    }
    # Colonne B
    resultatB <- search_id(data[i,5], thesaurus)
    if (length(resultatB) == 3) {
      ID <- resultatB[1]
      Proteine <- resultatB[2]
      Gene <- resultatB[3]
      data[i,5] <- ID
      data[i,3] <- Proteine
      data[i,12] <- Gene
    }
    else {
      # On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
      not_found <- data[i,c(5,3,12)]
      not_founds <- rbind(not_founds, not_found)
    }
    setTxtProgressBar(pb1, i)
  }
  
  # Memorisation du nouveau reseau (correction avec le thesaurus)
  network2 <- data
  
  ### Toutes les corrections automatiques sont faites sur le reseau ###
  
  
  # On recupere les identifiants que le thesaurus n'a pas trouves
  not_founds <<- unique(not_founds)
  if (is.null(dim(not_founds)[1]) == F) {
    cat(paste('\n\n>Search finished,', dim(as.matrix(not_founds))[1], 'IDs are not matched with thesaurus'))
    colnames (not_founds) <- c('uID', 'proteinname', 'genename')
    # Correction manuelle pour les identifiants que le thesaurus n'a pas trouves
    NF <- not_founds
    not.found.ID <- as.matrix(NF)
    
    # Mise en place de la fenetre d'affichage des corrections
    panelcorrection <- gwindow("Manual correction panel", parent = f_pos, visible = F, expand = T)
    
    pc <- ggroup(container = panelcorrection, horizontal = F, use.scrollwindow = T)
    pcsb <- ggroup(container = pc, horizontal = T, use.scrollwindow = F)
    lg <- gvbox(container = pcsb)
    pg <- gvbox(container = pcsb)
    rg <- gvbox(container = pcsb)
    fllg <- gformlayout(container = lg)
    flpg <- gformlayout(container = pg)
    flrg <- gformlayout(container = rg)
    
    # Affichage des identifiants a modifier
    for (i in 1:dim(not.found.ID)[1]) {
      uniprotID <- gedit(initial.msg = 'New UniprotID', label = as.character(not.found.ID[i,1]), container = fllg)
      PROTEINname <- gedit(initial.msg = 'New Protein name', label = as.character(not.found.ID[i,2]), container = flpg)
      GENEname <- gedit(initial.msg = 'New Gene name', label = as.character(not.found.ID[i,3]), container = flrg)
    }
    visible(panelcorrection)  <- T
    
    # Informations sur la correction manuelle
    Info <- '\"Correct Network\" : This button will save your manual corrections in the network and in updated databases.\n\n\"Ignore\" : This step will be ignored and the interactions with the uncorrected proteins will be conserved.'
    
    bpc <- ggroup(container = pc); addSpring(bpc)
    
    # 1 : Informations sur la correction manuelle
    bouton1 <- gbutton("Info", handler = function(h,...) {
      winfo <- gwindow("Info..", parent = f_pos)
      g <- gvbox(container = winfo); g$set_borderwidth(10L)
      glabel(Info, container = g)
      gseparator(container = g)
      bg <- ggroup(container = g); addSpring(bg)
      gbutton("Return", container = bg, handler = function(...) dispose(winfo))
    }, container = bpc)
    
    # 2 : Ignorer la correction manuelle, conservation des identifiants
    bouton2 <- gbutton('Ignore', handler = function(h,...) {
      # Recuperation des resultats finaux et sauvegardes
      finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
      dispose(bpc)
    }, container = bpc)
    
    # 3 : Correction manuelle du reseau
    bouton3 <- gbutton("Correct network", handler = function(h,...) {
      # Rassemblement des modifications a apporter
      cor.uniprotID <- cbind(names(svalue(fllg)), svalue(fllg))
      cor.uniprotID <- data.frame(cor.uniprotID, row.names = NULL)
      cor.proteinname <- cbind(names(svalue(flpg)), svalue(flpg))
      cor.proteinname <- data.frame(cor.proteinname, row.names = NULL)
      cor.genename <- cbind(names(svalue(flrg)), svalue(flrg))
      cor.genename <- data.frame(cor.genename, row.names = NULL)
      cor.manuel <- cbind(cor.uniprotID, cor.proteinname, cor.genename)
      colnames(cor.manuel) <-  c('old_uid', 'corrected_uid', 'old_proteinname', 'corrected_proteinname', 'old_genename', 'corrected_genename')
      # Corrections du reseau et ajout a la liste des corrections
      for (j in 1:dim(network2)[1]) {
        for (i in 1:dim(cor.manuel)[1]) {
          # Correction proteine A
          if (network2[j,4] == cor.manuel[i,2]) {
            network2[j,4] <- cor.manuel[[i,4]]
            network2[j,1] <- cor.manuel[[i,8]]
            network2[j,11] <- cor.manuel[[i,12]]
          }
          # Correction proteine B
          if (network2[j,5] == cor.manuel[i,2]) {
            network2[j,5] <- cor.manuel[[i,4]]
            network2[j,3] <- cor.manuel[[i,8]]
            network2[j,12] <- cor.manuel[[i,12]]
          }
        }
      }
      # Recuperation des resultats finaux et sauvegardes
      finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
      dispose(bpc)		
    }, container = bpc)
  }
  
  else if (is.null(not_founds)[1] == T) {
    cat(paste('\n\n>Search finished, all IDs are matched with thesaurus '))
    # Recuperation des resultats finaux et sauvegardes
    finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
  }
  else {
    cat(paste("\n\n>Search finished, 1 IDs isn't matched with thesaurus "))
    # Recuperation des resultats finaux et sauvegardes
    finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
  }
  
  setwd(mainpath)
  visible(mainpanel) <- T
  
}


build_network_window <- function(f_pos, mainpanel, mainpath) {

	db <- c()
	return.parameter <- c()
	ipl <- c()
	th <- c()

	panel_para <- gwindow("Build network : ", parent = f_pos, visible = T)
	pp <- gvbox(container = panel_para)
	pp$set_borderwidth(10L)

	flyt <- gformlayout(container = pp, expand = TRUE)

	# Selection des options de construction du reseau
	gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)
	gradio(c("Physical", "Genetic"), selected = 1, horizontal = T, label = "Experimental method", container = flyt)
	gradio(c("First-degree", "Second-degree"), selected = 1, horizontal = TRUE, label = "Interaction type", container = flyt)
	gradio(c("Yes", "No"), selected = 1, horizontal = TRUE, label = "Remove second degree unique links", container = flyt)
	gradio(c("Yes", "No"), selected = 2, horizontal = TRUE, label = "Remove all unique links", container = flyt)
	gradio(c("Yes", "No"), selected = 1, horizontal = TRUE, label = "Remove self-interactants", container = flyt)
	gradio(c("Yes", "No"), selected = 1, horizontal = TRUE, label = "Update databases if necessary", container = flyt)

	# Selection des bases de donnees
	chdb <- ggroup(container = pp, horizontale = T)
	addSpring(chdb)
	bouton1 <- gbutton("Select database", container = chdb, handler = function(...) {
		db <<- gfile(text = "Select database", type = "open", multi = T, container = chdb)
		if (is.null(db) == T) {
			gmessage('Selected database is null', icon = 'error')
		}
		if (is.null(db) == F) {
			bouton1$set_value(paste(length(db), 'databases selected'))
			cat(paste('\n>Database selected : ', db, sep = ''))
		}
	})
	# Selection du fichier contenant les noms des proteines d'interet
	bouton2 <- gbutton("Select input list", container = chdb, handler = function(...) {
		ipl <<- gfile(text = "Select input list", type = "open", multi = T,container = chdb)
		if (is.null(ipl) == T) {
			gmessage('Selected input list is null', icon = 'error')
		}
		if (is.null(ipl) == F) {
			bouton2$set_value(paste(length(ipl), 'input list selected'))
			cat(paste('\n>Inputlist selected : ', ipl, sep = ''))
		}
	})
	# Selection du fichier thesaurus
	bouton3 <- gbutton("Input thesaurus", container = chdb, handler = function(...) {
		th <<- gfile(text = "Select a thesaurus", type = "open", multi = F, container = chdb)
		if (is.null(th) == T) {
			gmessage('Selected thesaurus is null', icon = 'error')
		}
		if (is.null(th) == F) {
			bouton3$set_value(paste(length(th), 'thesaurus selected'))
			cat(paste('\n>Thesaurus selected : ', th, sep = ''))
		}
	})



	ppb <- ggroup(container = pp)
	addSpring(ppb)
	gbutton("Build the network", handler = function(h,...) {

		return.parameter <<- svalue(flyt)
		visible(panel_para) <- F
		# Memorisation des parametres de recherche selectionnes
		organism <- as.character(return.parameter[1])
		method <- return.parameter[2]
		degree <- as.character(return.parameter[3])
		remove.sdi <- as.character(return.parameter[4])
		remove.ul <- as.character(return.parameter[5])
		remove.si <- as.character(return.parameter[6])
		update <- as.character(return.parameter[7])

		# Verification de la presence de tous les elements necessaires a la construction du reseau
		# Lancement de la construction du reseau avec les parametres et fichiers donnes
		if (is.null(db) == F && is.null(ipl) == F && is.null(th) == F) {

			#############################################################################################################
			############################################## build_network() ##############################################

			### Construction du reseau

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

					build_network_othero(organism, db, ipl, th, method, degree, remove.sdi, remove.ul, remove.si, update, mainpath, f_pos)

					dispose(bpc)
				}, container = bpc)

			}

			else {

				cat("\n\n>BUILD NETWORK")

				# Parametres de la fonction
				Os <- organism
				os <- tolower(organism)
				ex.type <- tolower(method)
				inter.type <- tolower(degree)
				r.i <- tolower(remove.sdi)
				r.s.i <- tolower(remove.si)
				r.u.l <- tolower(remove.ul)
				IPL <- ipl
				DB <- db
				Th <- th
				UpDate <- tolower(update)

				# Lecture des fichiers input
				selected.database <- load_data(db)
				input.uniprotID <- read.table(file = ipl, header = FALSE, sep = '\t', stringsAsFactors = F)
				Thesaurus <- read.table(file = th, header = TRUE, quote = "", sep = "\t", stringsAsFactors = F)

				### Tous les parametres et fichiers fournis peuvent etre utilises ###


				# Completer la liste d'ID a rechercher dans les bases : ajout des anciens ID et isoformes ID
				listsupp <- c()
				# Parcours de tout le fichier thesaurus pour trouver tous les elements de la liste input et ajouter tous les autres uniprot-ID associe a la meme proteine
				for (i in 1:length(Thesaurus[,1])) {
					if (Thesaurus[i,1] %in% input.uniprotID[,1]) {
						oldNames <- Thesaurus[i,7]
						oldNamesList <- unlist(strsplit(oldNames, ";"))
						if (is.na(Thesaurus[i,7]) == FALSE) {
							# Ajout des anciens uniprot-ID a la liste
							for (j in 1:length(oldNamesList[])) {
								listsupp <- rbind(listsupp, oldNamesList[j])
							}
						}
						Nisoformes <- Thesaurus[i,8]
						if (Nisoformes > 0) {
							# Ajout des uniprot-ID des isoformes associes au meme gene
							for (k in 1 : Nisoformes) {
								newID <- paste(Thesaurus[i,1], k, sep = "-")
								listsupp <- rbind(listsupp, paste(Thesaurus[i,1], k, sep = "-"))
							}
						}
					}
				}
				# Apres avoir recupere tous les uniprot-ID associes a ceux de la liste input, les uniprot-ID de la liste input sont ajoutes
				input <- unlist(input.uniprotID)
				for (i in 1:length(input)) {
					uniprotID <- input[i]
					listsupp <- rbind(listsupp, uniprotID)
				}

				input.uniprotID <- listsupp

				### La liste des uniprot-ID a rechercher dans les bases est complete ###


				# Creation d'un dossier pour les fichiers resultats en sortie
				organismdir <- gsub(" ", "-", Os)
				organism.path <- paste(mainpath, organismdir, sep = '/')
				dir.create(organism.path, showWarnings = FALSE)
				network.path <- paste(organism.path, '/Network', sep = '')
				dir.create(network.path, showWarnings = FALSE)


				# Affichage des parametres selectionnes pour la recherche d'interactions
				cat(paste('\n\n>Species selected :', os))
				cat(paste('\n>Experimental type :', ex.type))
				cat(paste('\n>Interaction type :', inter.type))

				# Filtres des resultats sur l'organisme recherche et la methode souhaitee, dans les bases de donnees
				selected.database2 <- selected.database[(length(grep(os, selected.database[,8], ignore.case=TRUE)) > 0 && length(grep(os, selected.database[,9], ignore.case = TRUE)) > 0),]
				if (method == "genetic") {
					selected.database3 <- selected.database2[(length(grep(ex.type, selected.database2[,10], ignore.case = TRUE)) > 0 || length(grep(ex.type, selected.database2[,10], ignore.case = TRUE)) > 0),]
				}
				else {
					selected.database3 <- selected.database2[selected.database2[,10] == ex.type,]
				}
 				# Verification de la presence d'interactions dans la base
				if (dim(selected.database3)[1] == 0) {
					cat("\n\nWARNING : No result found in database(s)\nProcesuss stop")
					stop()
				}

				# Rassemblement de toutes les databases fournies
				selected.database4 <- as.matrix(data.frame(selected.database3[,1:15], stringAsFactors = F))

				### Les bases sont filtrees sur la recherche souhaitee ###


				# Recherche des interactions directes pour les proteines de l'input list
				cat('\n\n>Search for first degree interactions')
				PPI.Direct <- recup_ppi(input.uniprotID, selected.database3)
				PPI.Direct[!duplicated(PPI.Direct[,c(1:4,6,11:12)]),]
				# Verification de la presence d'interactions
				if (dim(PPI.Direct)[1] == 0) {
					cat("\n\nWARNING : No interaction found\nProcesuss stop")
					stop()
				}

				### Les interactions avec les proteines de l'input list sont extraites des bases ###


				# S'il y a bien des interactions identifiees dans les bases

				# Recherches du nombre d'article par interaction
				unredundant <- pubmed_id(PPI.Direct, os,1)
				PPI.Direct <- unredundant
				# Si le reseau demande est de degre 2 on relance une recherche d'interactions
				if (inter.type == "second-degree") {
					# Recuperation des uniprot-ID de toutes les proteines du reseau de degre 1
					Interactions.matrix <- PPI.Direct
					# Recuperation des proteines d'interet et de leurs interactants directes
					col1 <- as.matrix(Interactions.matrix[,1])
					col2 <- as.matrix(Interactions.matrix[,2])
					# Rassemblement des listes pour recuperer les uniprot-ID en une seule copie
					listProt <- rbind(col1, col2)
					listProt <- unique(listProt)
					ListProt.UniprotID <- as.matrix(listProt)

					# Seconde recherche d'interactions, pour le degre 2
					cat('\n>Search for second degree interactions')
					PPI.Direct2 <- recup_ppi(ListProt.UniprotID, selected.database3)
	
					# Suppression des interactions de degre 2 avec les proteines qui n'ont qu'une seule interaction dans le reseau de degre 1 : si demande
					if (r.i == 'yes') {
						PPI.Indirect <- remove_ppi(ListProt.UniprotID, PPI.Direct2, input.uniprotID)
					}
					else {
						PPI.Indirect <- PPI.Direct2
					}
					PPI.Indirect<-unique(PPI.Indirect)
					PPI.Indirect[!duplicated(PPI.Indirect[,c(1:4,6,11:12)]),]
					cat('\n>Search for PUBMED IDS for each interaction')
					
					# Recherches du nombre d'article par interaction
					unredundant <- pubmed_id(PPI.Indirect, os,2)
					PPI.Indirect<-unredundant

					# Rassemblement des resultats d'interactions de degre 1 et degre 2
					network <- rbind(PPI.Direct[,c(3,5,4,1,2,6:13)], PPI.Indirect[,c(3,5,4,1,2,6:13)])

					### Les interactions de degre 2 sont ajoutees au reseau ###

				}

				else {
					# On ne garde que les resultats de la premiere recherche d'interactions si on souhaite un reseau de degre 1
					network <- rbind(PPI.Direct[,c(3,5,4,1,2,6:13)])
				}

				# Verification de la presence de resultats
				if (is.null(network) == T) {
					cat("WARNING :\nThe network is null!!!\nProcessus stoped")
					stop()
				}

				cat('\n')

				# Premier trie des interactions du reseau en fonction des parametre de construction selectionnes
				# Remove redundants
				cat('\n>Removing redundants ...\n')
				network<-DataBases(network)
				selected.database3 <- remove_redundants(network)
				network <- as.matrix(selected.database3)

				# Remove self-interactant
				if (r.s.i == 'yes'){
					cat('\n>Removing self-interactants  ...\n')
					network <- network[network[,4]!=network[,5],]
					network <<- network
				  cat(' OK')
				}
				else {
					cat('\n>Proteins which interact with itselves are kept')
				}

				### Le reseau est complet et un premier trie a ete effectue ###


				# Correction automatique du reseau avec le thesaurus : on met l'uniprot-ID principal pour chaque proteine (a la place des anciens)
				data <- as.matrix(network)
				thesaurus <- read.table(file = th, header = TRUE, quote = "", sep = "\t" , stringsAsFactors = F)
				# Visualisation de la progression des corrections
				cat ( '\n>Autocorrection' )
				cat ( '\n' )
				pb1 <<- txtProgressBar(min = 0, max = dim(data)[1], style = 3)
				nbpassage<<-0
				nf<<-c()
				apply(data,1,cherche_uniprotID,data=data,thesaurus=thesaurus)
				not_founds<-nf
				  # Parcours de toutes les interactions du reseau
#				for (i in 1:dim(data)[1]) {
				# Recherche de l'uniprot ID de reference et association du bon nom de proteine et de gene
					# Colonne A
#				  resultatA <- search_id(data[i,4], thesaurus)
#					if (length(resultatA) == 3) {
#					  ID <- resultatA[1]
#						Proteine <- resultatA[2]
#						Gene <- resultatA[3]
#						data[i,4] <- ID
#						data[i,1] <- Proteine
#						data[i,11] <- Gene
#					}
#				  v1<-c(grep(data[i,11],data[,11]))
#				  v2<-c(grep(data[i,11],data[,12]))
#				  if((length(unique(c(data[v1,4],data[v2,5])))>1)||(length(resultatA) != 3))
#					{
					# On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
#					  not_found <- data[i,c(4,1,11)]
#					  not_founds <- rbind(not_founds, not_found)
#					}
					# Colonne B
#				  resultatB <- search_id(data[i,5], thesaurus)
#					if (length(resultatB) == 3) {
#					  ID <- resultatB[1]
#						Proteine <- resultatB[2]
#						Gene <- resultatB[3]
#						
#						data[i,5] <- ID
#						data[i,3] <- Proteine
#						data[i,12] <- Gene
#					}
#				  v3<-c(grep(data[i,12],data[,11]))
#				  v4<-c(grep(data[i,12],data[,12]))
#				  if((length(unique(c(data[v3,4],data[v4,5])))>1)||(length(resultatB) != 3))
#					{
					# On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
#					  not_found <- data[i,c(5,3,12)]
#						not_founds <- rbind(not_founds, not_found)
#					}
#					setTxtProgressBar(pb1, i)
#				}#end for
				
				
				# Memorisation du nouveau reseau (correction avec le thesaurus)
				network2 <- data

				### Toutes les corrections automatiques sont faites sur le reseau ###

        
				# On recupere les identifiants que le thesaurus n'a pas trouves
				not_founds <<- unique(not_founds)
				
				if (is.null(dim(not_founds)[1]) == F) {
				  cat(paste('\n\n>Search finished,', dim(as.matrix(not_founds))[1], 'IDs are not matched with thesaurus'))
					colnames (not_founds) <- c('uID', 'proteinname', 'genename')
					# Correction manuelle pour les identifiants que le thesaurus n'a pas trouves
					NF <- not_founds
				  not.found.ID <- as.matrix(NF)
	
  				# Mise en place de la fenetre d'affichage des corrections
	  			panelcorrection <- gwindow("Manual correction panel", parent = f_pos, visible = F, expand = T)
	
		  		pc <- ggroup(container = panelcorrection, horizontal = F, use.scrollwindow = T)
			  	pcsb <- ggroup(container = pc, horizontal = T, use.scrollwindow = F)
				  lg <- gvbox(container = pcsb)
					pg <- gvbox(container = pcsb)
					rg <- gvbox(container = pcsb)
				  fllg <- gformlayout(container = lg)
					flpg <- gformlayout(container = pg)
					flrg <- gformlayout(container = rg)
					not.found.ID<-unique(not.found.ID)
					# Affichage des identifiants a modifier
					for (i in 1:dim(not.found.ID)[1]) {
					  uniprotID <- gedit(initial.msg = 'New UniprotID', label = as.character(not.found.ID[i,1]), container = fllg)
					  PROTEINname <- gedit(initial.msg = 'New Protein name', label = as.character(not.found.ID[i,2]), container = flpg)
						GENEname <- gedit(initial.msg = 'New Gene name', label = as.character(not.found.ID[i,3]), container = flrg)
				  }
					visible(panelcorrection)  <- T
	
  			  # Informations sur la correction manuelle
	  			Info <- '\"Correct Network\" : This button will save your manual corrections in the network and in updated databases.\n\n\"Ignore\" : This step will be ignored and the interactions with the uncorrected proteins will be conserved.'
  
		  		bpc <- ggroup(container = pc); addSpring(bpc)
	
			  	# 1 : Informations sur la correction manuelle
				  bouton1 <- gbutton("Info", handler = function(h,...) {
					  winfo <- gwindow("Info..", parent = f_pos)
						g <- gvbox(container = winfo); g$set_borderwidth(10L)
						glabel(Info, container = g)
						gseparator(container = g)
						bg <- ggroup(container = g); addSpring(bg)
						gbutton("Return", container = bg, handler = function(...) dispose(winfo))
				  }, container = bpc)
		 
					# 2 : Ignorer la correction manuelle, conservation des identifiants
					bouton2 <- gbutton('Ignore', handler = function(h,...) {
					  # Recuperation des resultats finaux et sauvegardes
					  finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
						dispose(bpc)
				  }, container = bpc)
				  #   3 : Correction manuelle du reseau
					bouton3 <- gbutton("Correct network", handler = function(h,...) {
					  # Rassemblement des modifications a apporter
					  if(svalue(fllg)!=""){
						  cor.uniprotID <- cbind(names(svalue(fllg)), svalue(fllg))
					  }
					  else{
					    cor.uniprotID <- cbind(names(svalue(fllg)), names(svalue(fllg)))
					  }
					  if(svalue(flpg)!=""){
						  cor.proteinname <- cbind(names(svalue(flpg)), svalue(flpg))
					  }
					  else{
					    cor.proteinname <- cbind(names(svalue(flpg)), names(svalue(flpg)))
					  }
					  if(svalue(flrg)!=""){
						  cor.genename <- cbind(names(svalue(flrg)), svalue(flrg))
					  }
					  else{
					    cor.genename <- cbind(names(svalue(flrg)), names(svalue(flrg)))
					  }
					  cor.uniprotID <- data.frame(cor.uniprotID, row.names = NULL)
					  cor.proteinname <- data.frame(cor.proteinname, row.names = NULL)
					  cor.genename <- data.frame(cor.genename, row.names = NULL)
						cor.manuel <- cbind(cor.uniprotID, cor.proteinname, cor.genename)
						cor.manuel <- data.frame(cor.manuel, row.names = NULL)
						#
						#
						colnames(cor.manuel) <-  c('old_uid', 'corrected_uid', 'old_proteinname', 'corrected_proteinname', 'old_genename', 'corrected_genename')
						
						# Corrections du reseau et ajout a la liste des corrections
						for (j in 1:dim(network2)[1]) {
						  for (i in 1:dim(cor.manuel)[1]) {
							  # Correction proteine A
								if (network2[j,4] == cor.manuel[i,1]) {
								  network2[j,4] <- gsub("\n", "", cor.manuel[[i,2]][1])
								}
						    if (network2[j,1] == cor.manuel[i,4]) {
						      network2[j,1] <- gsub("\n", "", cor.manuel[[i,4]][1])
						    }
						    if (network2[j,11] == cor.manuel[i,6]) {
						      network2[j,11] <-gsub("\n", "", cor.manuel[[i,6]][1])
						    }
						    
							  # Correction proteine B
								if (network2[j,5] == cor.manuel[i,1]) {
								  network2[j,5] <- gsub("\n", "", cor.manuel[[i,2]][1])
								  }
						    if (network2[j,3] == cor.manuel[i,4]) {
									network2[j,3] <-  gsub("\n", "", cor.manuel[[i,4]][1])
						    }
						    if (network2[j,12] == cor.manuel[i,6]) {
									network2[j,12] <- gsub("\n", "", cor.manuel[[i,6]][1])
						    }	
						  }
						}
						# Recuperation des resultats finaux et sauvegardes
						finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
						  dispose(bpc)		
					}, container = bpc)
				}

				else if (is.null(not_founds)[1] == T) {
				  cat(paste('\n\n>Search finished, all IDs are matched with thesaurus '))
				  # Recuperation des resultats finaux et sauvegardes
			    finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
				}
				else {
					cat(paste("\n\n>Search finished, 1 IDs isn't matched with thesaurus "))
				  # Recuperation des resultats finaux et sauvegardes
					finish(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path)
				}
				
				setwd(mainpath)
			  visible(mainpanel) <<- T

			}

			dispose(panel_para)
			dispose(ppb)

		}
	
		# Affichage d'un message d'erreur s'il manque un fichier
		else if (is.null(db) == T) {
			gmessage('Database selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		
		else if (is.null(ipl) == T) {
			gmessage('Input list selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else if (is.null(th) == T) {
			gmessage('Thesaurus selected is null', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}
		else {
			gmessage('Error : Unable to start search', icon = 'error')

			dispose(ppb)
			visible(panel_para) <- F
			visible(mainpanel) <<- T

		}

	}, container = ppb)

	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(mainpanel) <<- T
	}, container = ppb)
	
	visible(panel_para) <- T

}

