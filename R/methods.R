### function construct
### function branch_length

#############################################################

construct <- function(part1, part2, e, mat, inputID, inputSGD) {
	ind1 <- grep(part1, e)
	ind2 <- grep(part2, e)
	if ((length(ind1) > 0) && (length(ind2) > 0) && (ind1 < ind2)) {	
		ind2 <- ind2-1
	}
	if (length(ind1) > 0) {
		part1 <- e[ind1]
		e <- e[-ind1]
	}
	if (length(ind2) > 0) {
		part2 <- e[ind2]
		e <- e[-ind2]
	}	
	if (length(part1) == length(setdiff(inputID, inputID[inputSGD == part2]))) {
		d1 <- round(branch_length(mat, part1, setdiff(inputID, inputID[inputSGD == part1])), 3) / 2.0
		d2 <- round(branch_length(mat, part2, setdiff(inputID, inputID[inputSGD == part2])), 3) / 2.0
	} 
	else {
		d1 <- round(branch_length(mat, part1, setdiff(inputID, inputID[inputSGD == part1])), 3)
		d2 <- round(branch_length(mat, part2, setdiff(inputID, inputID[inputSGD == part2])), 3)
	}
	
	e <- c(e, paste("(", part1, ":", d1, ",", part2, ":", d2, ")", sep = ""))
	return(e)    

}

#############################################################

branch_length <- function(matrix, set1, set2) {

	goodTriples <- 0
	badTriples <- 0
	# computing triplets with one element in set2 and 2 elements in set1
	if (length(set1) > 1) {
		for (i in 1:length(set1)) {
			for (j in i:length(set1)) {
				for (k in 1:length(set2)) {
					if (k != i && k != j && i != j) {	
						if (matrix[i,j] < max(matrix[i,k], matrix[j,k])) {
							badTriples <- badTriples + 1
						}
						else {
							goodTriples <- goodTriples + 1
						}
					}
				}
			}
		}
	}
	# computing triplets with one element in set1 and 2 elements in set2
	if (length(set2) > 1) {
		for (i in 1:length(set2)) {
			for (j in i:length(set2)) {
				for (k in 1:length(set1)) {
					if (k != i && k != j && i != j) {	
						if (matrix[i,j] < max(matrix[i,k], matrix[j,k])) {
							badTriples <- badTriples + 1
						}
						else {
							goodTriples <- goodTriples + 1
						}
					}
				}
			}
		}
	}

	return (goodTriples / (goodTriples + badTriples));

}

# Search all the databases having the same PPI coming from one or several publications 

DataBases<-function(Final.List.Redondant){
  
  Inter.red<-as.matrix(Final.List.Redondant)
  
  indices<-rep("NA",dim(Inter.red)[1])

  for(i in 1:dim(Inter.red)[1])
  {
    indi<-intersect(c(grep(Inter.red[i,12],Inter.red[,11]),grep(Inter.red[i,12],Inter.red[,12])),c(grep(Inter.red[i,11],Inter.red[,11]),grep(Inter.red[i,11],Inter.red[,12])))
    indices[i]<-toString(unique(Inter.red[indi,10]))
  }
  
  Inter.red[,10]<-indices
  
  return(Inter.red)
}

search_id <- function(cible, thesaurus) {
  
  # Recherche de la cible dans tout le thesaurus pour associer le bon UniprotID avec le nom de la proteine et le nom du gene
  
  # On regarde si la cible peut etre un ID-isoforme pour cherche l'ID seul dans le thesaurus : UniprotID ou OldID
  CIBLE <- unlist(strsplit(cible, "-"))
  ciblePrincipale <- CIBLE[1]
  if (length(CIBLE[]) > 1) {
    complementCible <- paste('-', CIBLE[2], sep = '')
  }
  else {
    complementCible <- ""
  }
  
  # On cherche la cible dans la liste des uniprotID du thesaurus : colonne 1
  ligneThesaurus <- 1 # on commence à la premiere ligne du thesaurus
  ligneMax <- dim(thesaurus)[1]
  resultat <- list()
  recherche <- 'non'
  while (recherche == 'non' & ligneThesaurus <= ligneMax) {
    if (ciblePrincipale == thesaurus[ligneThesaurus,1]) {
      resultat <- c(paste(thesaurus[ligneThesaurus,1], complementCible, sep = ''), paste(thesaurus[ligneThesaurus,5], complementCible, sep = ''), thesaurus[ligneThesaurus,3])
      recherche <- 'oui'
    }
    else {
      ligneThesaurus <- ligneThesaurus + 1
    }
  }
  
  # Si la cible n'a pas ete trouvee dans les uniprotID on cherche dans les anciens ID : colonne 7
  if (length(resultat) == 0) {
    ligneThesaurus <- 1
    resultat <- list()
    recherche <- 'non'
    while (recherche == 'non' & ligneThesaurus <= ligneMax) {
      test <- grep(ciblePrincipale, thesaurus[ligneThesaurus,7])
      if (length(test) != 0) {
        resultat <- c(paste(thesaurus[ligneThesaurus,1], complementCible, sep = ''), paste(thesaurus[ligneThesaurus,5], complementCible, sep = ''), thesaurus[ligneThesaurus,3])
        recherche <- 'oui'
      }
      else {
        ligneThesaurus <- ligneThesaurus + 1
      }
    }
  }
  
  # Si la cible n'a pas ete trouvee on cherche dans les Gene ID : colonne 6
  if (length(resultat) == 0) {
    ligneThesaurus <- 1
    resultat <- list()
    recherche <- 'non'
    while (recherche == 'non' & ligneThesaurus <= ligneMax) {
      if (is.na(thesaurus[ligneThesaurus,6]) == FALSE & cible == thesaurus[ligneThesaurus,6]) {
        resultat <- c(thesaurus[ligneThesaurus,1], thesaurus[ligneThesaurus,5], thesaurus[ligneThesaurus,3])
        recherche <- 'oui'
      }
      else {
        ligneThesaurus <- ligneThesaurus + 1
      }
    }
  }
  
  # Si la cible n'a pas ete trouvee on cherche dans les RefSeq : colonne 4
  if (length(resultat) == 0) {
    ligneThesaurus <- 1
    resultat <- list()
    recherche <- 'non'
    while (recherche == 'non' & ligneThesaurus <= ligneMax) {
      if (is.na(thesaurus[ligneThesaurus,4]) == FALSE & cible == thesaurus[ligneThesaurus,4]) {
        resultat <- c(thesaurus[ligneThesaurus,1], thesaurus[ligneThesaurus,5], thesaurus[ligneThesaurus,3])
        recherche <- 'oui'
      }
      else {
        ligneThesaurus <- ligneThesaurus + 1
      }
    }
  }
  
  # Si la cible n'a pas ete trouvee on cherche dans les Protein Name : colonne 5
  if (length(resultat) == 0) {
    ligneThesaurus <- 1
    resultat <- list()
    recherche <- 'non'
    while (recherche == 'non' & ligneThesaurus <= ligneMax) {
      if (cible == thesaurus[ligneThesaurus,5]) {
        resultat <- c(thesaurus[ligneThesaurus,1], thesaurus[ligneThesaurus,5], thesaurus[ligneThesaurus,3])
        recherche <- 'oui'
      }
      else {
        ligneThesaurus <- ligneThesaurus + 1
      }
    }
  }
  
  # Si la cible n'a pas ete trouvee on cherche dans les Gene Name : colonne 3
  if (length(resultat) == 0) {
    ligneThesaurus <- 1
    resultat <- list()
    recherche <- 'non'
    while (recherche == 'non' & ligneThesaurus <= ligneMax) {
      if (cible == thesaurus[ligneThesaurus,3]) {
        resultat <- c(thesaurus[ligneThesaurus,1], thesaurus[ligneThesaurus,5], thesaurus[ligneThesaurus,3])
        recherche <- 'oui'
      }
      else {
        ligneThesaurus <- ligneThesaurus + 1
      }
    }
  }
  
  return (resultat)
  
}

############################################################################################

recup_ppi <- function(inputListFile, Base) {
  
  # Recuperation des colonnes des bases qui sont utilisees dans le reseau
  Base.f <- as.matrix(Base[,c(1:5,7:11,14:15)])
  ind <- unique(c(as.vector(grep(inputListFile[1], Base.f[,1])), as.vector(grep(inputListFile[1], Base.f[,2]))))
  
  cat('\n>Searching interactions...')
  
  # Recherche des interactions dans les bases, contenant au moins une proteine du fichier de recherche input list
  for (i in 1:length(inputListFile)) {
    ind <- unique(c(ind, c(as.vector(grep(inputListFile[i], Base.f[,1])), as.vector(grep(inputListFile[i], Base.f[,2])))))
  }
  interaction.dir <- Base.f[ind,]
  
  cat('OK')
  
  return(interaction.dir)
  
}

############################################################################################

pubmed_id <- function(Final.List.Redondant, Name, run) {
  
  cat('\n>Searching pubmed IDs ... ')
  print(dim(Final.List.Redondant))
  cat("\n")
  
  # Recuperation des donnees du reseau
  Inter.red <- as.matrix(Final.List.Redondant)
  Inter.partialRed <- as.matrix(Final.List.Redondant[1,])
  if (dim(Inter.partialRed)[2] == 1) {
    Inter.partialRed <- t(Inter.partialRed)
  }
  if(run==2)
    pb <- txtProgressBar(min = 0, max = dim(Inter.red)[1],style = 3)
  # Parcours du reseau pour trouver des interactions redondantes et identifier le nombres d'articles associes (nombre de pubmed ID)
  l <- dim(Inter.red)[1]
  i <- 1
  if (Name != "nofile") {
    unredundant <<- as.matrix(t(c(Inter.red[1,], "NA")))
    colnames(unredundant) <<- c(colnames(Inter.red), "NbPmids")
    nbl <- 1
    while (i <= l) {
      indi <- intersect(c(grep(Inter.red[i,2], Inter.red[,1]), grep(Inter.red[i,2], Inter.red[,2])), c(grep(Inter.red[i,1], Inter.red[,1]), grep(Inter.red[i,1], Inter.red[,2])))
      j <- 1
      # Recherche des redondances d'interactions
      while (j <= length(indi)) {
        if ((((Inter.red[indi[j],1] == Inter.red[i,1]) && (Inter.red[indi[j],2] == Inter.red[i,2])) || ((Inter.red[indi[j],1] == Inter.red[i,2]) && (Inter.red[indi[j],2] == Inter.red[i,1]))) == FALSE) {
          indi <- indi[-j]
          j <- j - 1
        }
        j <- j + 1
      }
      if (nbl > 1) {
        unredundant <<- rbind(unredundant, as.matrix(t(c(Inter.red[i,], "NA"))))
      }
      # Recuperation des differents pubmed-ID
      if (length(indi) > 1) {
        pmids <- unique(Inter.red[indi,6])
        pmids <- as.matrix(pmids)
        pmids2 <- c()
        for (n in 1:length(pmids)) {
          a <- strsplit(pmids[n,1], "\\|")
          a1 <- length(a[[1]])
          for (p in 1:a1) {
            pmids2 <- c(pmids2, a[[1]][p])
          }
        }
        pmids <- unique(as.vector(pmids2))
        # Rassemblement des different pubmed-ID de chaque interaction
        pmids3 <- c()
        b <- length (pmids)
        if (b > 1) {
          for (o in 1:length(pmids)) {
            pmids3 <- paste(pmids3, pmids[o], sep = "|")
          }
          pmids3 <- substr(pmids3, 2, nchar(pmids3))
        }
        if (b == 1) {
          pmids3 <- paste(pmids3, pmids)
        }
        tmp <- Inter.red[indi,]
        for (k in 1:length(pmids)) {
          Inter.partialRed <- rbind(Inter.partialRed, t(as.matrix(tmp[(1:dim(tmp)[1])[tmp[,6] == pmids[k]][1],])))
        }
        unredundant[nbl,6] <<- pmids3
        # On complete la derniere colonne du reseau avec le nombre de pubmed-ID trouves pour chaque interaction
        unredundant[nbl,13] <<- length(pmids)
        Inter.red <- Inter.red[-indi,]
        l = l - length(indi)
        i = i - 1
      }
      if (length(indi) == 1) {
        unredundant[nbl,13] <<- 1
        Inter.partialRed <- rbind(Inter.partialRed, t(as.matrix(Inter.red[indi,])))
      }
      i = i + 1
      nbl <- nbl + 1
      
      cat("i :")
      print(i)
      cat("\n")
      if(run==2)
        setTxtProgressBar(pb, i)
    }
    Inter.partialRed <- Inter.partialRed[-1,]
    if(run==2)
      close(pb)
  }
  cat('OK')
  
  return(unredundant)
  
}

############################################################################################

load_data <- function(db) {
  
  # Recuperation des donnes contenues dans les bases de donnees selectionnees
  cat('\n\n>Loading database...')
  data.name.list <- db
  
  # Verification de la presence de donnees dans la base
  if (length(data.name.list) == 0) {
    cat('ERROR : no databases selected')
    stop()
  }
  
  # Verification de la compatibilite du format des bases de donnees
  multiple.database <- c()
  for (i in 1:length(data.name.list)) {
    one.multiple.database <- read.delim2(data.name.list[i], header = T, sep = '\t')
    tryCatch({
      colnames(one.multiple.database) <- c( "uidA", "uidB", "aliasA", "aliasB", "method", "author", "pmid", "taxA", "taxB", "interactionType", "sourceBD", "confidence", "numParticipants", "GeneNameA", "GeneNameB")
      multiple.database <- rbind(multiple.database, one.multiple.database)
    }, error = function(err) {
      message(err)
      message('Database(s) dimension different')
    })
  }
  
  cat('OK dim :')
  cat(dim(multiple.database))
  # Recuperation des bases si le format correspond
  return(multiple.database)
  
}

############################################################################################

load_network <- function(nw) {
  
  # Recuperation des donnes contenues dans le reseau selectionne
  cat('\n>Loading network ... ')
  network.name.list <- nw
  
  # Verification de la presence d'interactions dans le reseau
  if (length(network.name.list) == 0) {
    cat('ERROR : no network selected')
    stop()
  }
  
  # Verification de la compatibilite de format du reseau
  multiple.network <- c()
  for (i in 1:length(network.name.list)) {
    one.multiple.network <- read.delim2(network.name.list[i], header = T, sep = '\t', stringsAsFactors = F)
    tryCatch({
      colnames(one.multiple.network)  <- c("aliasA", "method", "aliasB", "uidA", "uidB", "pmid", "taxA", "taxB", "interactionType", "sourceBD", "GeneNameA", "GeneNameB","NbPmids" )
      multiple.network <- rbind(multiple.network, one.multiple.network)
    }, error = function(err) {
      message(err)
      message('Network dimension uncorrect')
    })
  }
  
  cat('OK dim :')
  cat(dim(unique(multiple.network)))
  # Recuperation du reseau si le format correspond
  return(unique(multiple.network))
  
}

############################################################################################

proximity_score <- function (listProt1, listProt2, inputUniprotID, allProt, formula) {
  
  score <- (-1.0)
  listinter <- intersect(listProt1, listProt2)
  ninter <- length(listinter)
  listunion <- union(listProt1, listProt2)
  nunion <- length(listunion)
  nUniprotID <- length(inputUniprotID)
  for (i in 1:nUniprotID) {
    ninter <- ninter - length(grep(inputUniprotID[i], listinter))
    nunion <- nunion - length(grep(inputUniprotID[i], listunion))
  }
  
  o1 <- length(listProt1)
  o2 <- length(listProt2)
  O1 <- as.double(o1)
  O2 <- as.double(o2)
  OO <- as.double(O1 + O2)
  # proteins not in the complex, in the neighborhood of subcomplex 1 and subcomplex 2
  O11 <- as.double(ninter)
  O <- as.double(nunion)
  # proteins not in the complex, in the neighborhood of subcomplex 1 but not of subcomplex 2
  O12 <- as.double(O1 - O11)
  # proteins not in the complex, in the neighborhood of subcomplex 2 but not of subcomplex 1
  O21 <- as.double(O2 - O11)
  # O21<-length(listProt2)-length(intersect(listProt2,inputUniprotID))-length(ninter)
  # proteins neither in the complex, nor in the neighborhood of subcomplex 1 nor of subcomplex 2
  N <- length(allProt) - length(inputUniprotID)
  N <- as.double(N)
  O22 <- as.double(N - O)	
  S1 <- as.double(O12 + O22)
  S2 <- as.double(O21 + O22)
  E11 <- as.double(O1 * O2) / N
  E12 <- as.double(O1 * S1) / N
  E21 <- as.double(S2 * O2) / N
  E22 <- as.double(S2 * S1) / N
  
  if (formula == "jaccard") {
    valeur <- (O11 / O)
  }
  if (formula == "liddell") {
    valeur <- (O11 * O22 - O12 * O21) / (O2 * S1)
  }
  if (formula == "dice") {
    valeur <- 2.0 * O11 / (O1 + O2)
  }
  if (formula == "zscore") {
    valeur <- (O11 - E11) ^ 2.0 / sqrt(E11)
  }  
  if (formula == "ms") {
    valeur <- min(O11 / O1, O11 / O2)
  }
  if (formula == "Chi2") {
    valeur <- N * (O11 - E11) ^ 2.0 / (E11 * E22)
  }
  
  score <- round(valeur, 3)
  if (score != 'NA') {
    return (score)	
  }
  else {
    cat('Error : Unable to generate a score, choose an other one')
    stop()
  }
  
}

############################################################################################

normalize_mat <- function(mat) {
  
  matMin = min(mat, na.rm = TRUE)
  matMax = max(mat, na.rm = TRUE)
  if (matMin != matMax) {
    for (i in 1:(length(mat[1,]))) {
      for (j in 1:(length(mat[1,]))) {
        if (matMin < 0) {
          # affine transformation
          mat[i,j] <- 0.99 * (mat[i,j] - matMin) / (matMax - matMin)
        } 
        else {
          # linear transformation
          mat[i,j] <- 0.99 * mat[i,j] / matMax
        }
      }
    }
  }
  else {
    for (i in 1:(length(mat[1,]))) {
      for (j in 1:(length(mat[1,]))) {
        mat[i,j] <- 0 * mat[i,j]
      } 
    }	
  }
  
  return(mat)
  
}

############################################################################################

load_inputlist <- function(data) {
  
  # Recuperation des noms des proteines d'interet
  
  cat('\n>Loading inputlist ... ')
  
  data.name.list <- data
  if (length(data.name.list) == 0) {
    cat('Error no file selected')
    stop()
  }
  
  multiple.database <- c()
  for (i in 1:length(data.name.list)) {
    one.multiple.database <- read.delim2(data.name.list[i], header = T, sep = '\t')
    multiple.database <- rbind(multiple.database, one.multiple.database)
  }
  
  cat('OK dim :')
  cat(dim(unique(multiple.database)))
  return(unique(multiple.database))
  
}

############################################################################################

indices_min_jaccard <- function(mat.JaccardDistance) {
  
  ind <- c(-1, -1)
  for (i in 1:(dim(mat.JaccardDistance)[1] - 1)) {
    for (j in (i + 1):dim(mat.JaccardDistance)[2]) {
      if (mat.JaccardDistance[i,j] == min(mat.JaccardDistance, na.rm = TRUE)) {
        ind[1] <- i
        ind[2] <- j
      }
    }
  }
  
  return(ind)
  
}

############################################################################################

remove_redundants <- function(inputtable) {
  
  # Recuperation du reseau d'interactions
  tab_PPI <- as.matrix(inputtable)
  
  # On suit l'evolution du parcours
  prep_file <- tab_PPI[duplicated(tab_PPI[,c(4,5)]) == F,]
  
  # Dedoublement du reseau pour avoir les interactions dans les deux sens
  prep_file_inverse  <- as.matrix(data.frame(prep_file[,1:3], prep_file[,5], prep_file[,4], prep_file[,6:10], prep_file[,12], prep_file[,11], stringAsFactors=F))
  
  # Parcours du reseau
  i <- 1
  while (i <= dim(prep_file)[1]) {
    prsbar <- txtProgressBar(min = 1, max = dim(prep_file)[1], style = 3)
    j <- 1
    nom_pubmed_j <- c()
    while (j <= dim(prep_file_inverse)[1]) {
      # On enleve les lignes avec des interactions entre les meme proteines, sans enlever les interactions d'une proteine avec elle meme
      if (prep_file[i,4] == prep_file_inverse[j,4] && prep_file[i,5] == prep_file_inverse[j,5] && as.character(prep_file[i,4]) != as.character(prep_file[i,5])) {
        prep_file <- prep_file[-j,]
        prep_file_inverse  <-  prep_file_inverse[-j,]
      }
      else {
        j <- j + 1
      }
    }
    i <- i + 1
    setTxtProgressBar(prsbar, i)
  }
  
  close(prsbar)
  cat('>Redundants removed')
  
  return(prep_file)
  
}

############################################################################################

remove_unique_links <- function(inputtable) {
  
  cat('\n>Removing unique links ... ')
  
  tab_PPI <- as.matrix(inputtable)
  i <- 1
  while (i <= dim(tab_PPI)[1]) {
    # Suppression des interactions d'une proteine avec elle meme
    if (length(which(tab_PPI[i,4] == tab_PPI[,4:5])) == 1 || length(which(tab_PPI[i,5] == tab_PPI[,4:5])) == 1 ) {
      tab_PPI <- tab_PPI[-i,]
    }
    else {
      i <- i + 1
    }
  }
  
  cat('OK')
  
  return(tab_PPI)
  
}

############################################################################################

saving <- function(network2, not_founds, network.path, Os, r.i, r.s.i, r.u.l, ex.type, inter.type, os, IPL, Th, DB) {
  
  cat('\n\n>> Saving network files in : ')
  cat(network.path)
  cat(' directory :')
  
  # Recuperation des parametres de construction du reseau
  Cor <- "Network"
  OS <- gsub(" ", "-", Os)
  REMOVE <- ''
  if (r.i == TRUE) {
    REMOVE <- paste(REMOVE, "Second degree interactions for proteins with unique link in first degree network\n", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-\n", sep = '')
  }
  if (r.s.i == TRUE) {
    REMOVE <- paste(REMOVE, "Self interactions\n", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-\n", sep = '')
  }
  if (r.u.l == TRUE) {
    REMOVE <- paste(REMOVE, "All unique links", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-", sep = '')
  }
  if (is.null(dim(not_founds)[1]) == T) {
    not_founds <- "All IDs are found in thesaurus"
  }
  
  cat('\n>Saving network ... ')
  
  # Sauvegarde du reseau
  out.put.name <- paste(ex.type, Cor, OS, 'degree', inter.type, 'interactions.txt', sep = '_')
  write.table(network2, file = paste(network.path, out.put.name, sep = '/'), row.names = F, col.names = T, quote = F, sep = "\t")
  
  setwd(network.path)
  
  cat('OK\nNetwork : ')
  cat(paste(network.path, out.put.name, sep = '/'))
  cat('\n>Saving summary ... ')
  
  # Sauvegarde du recapitulatif des parametres utilises et des modifications apportees au reseau (et aux bases)
  out.put.name2 <- paste(OS, 'Summary_network.txt', sep = '_')
  write.table(paste("\n\n\n\n\nDATE :", Sys.time(), "NETWORK FILE NAME :", out.put.name, "ORGANISM :", os, "EXPERIENCE TYPE :", ex.type, "INTERACTIONS TYPE :", inter.type, "REMOVE :", REMOVE, "INPUTLIST :", IPL, "THESAURUS :", Th, "DATABASE(S) :", DB, "FUNCTION USED :\nbuild_network()", "ID(S) not found in thesaurus :", sep = '\n'), file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table(DB, file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table("\nFUNCTION USED :\nbuild_network()\nID(S) not found in thesaurus :\n", file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table(not_founds, file = paste(network.path, out.put.name2, sep = '/'), append = T)
  
  cat('OK\nSummary : ')
  cat(paste(network.path, out.put.name2, sep = '/'))
  cat('\n\n>Construction of network is done.')
  
}


############################################################################################
saving_window <- function(network2, network.path, Os, r.i, r.s.i, r.u.l, not_founds, ex.type, inter.type, os, IPL, Th, DB) {
  
  cat('\n\n>> Saving network files in : ')
  cat(network.path)
  cat(' directory :')
  
  # Recuperation des parametres de construction du reseau
  Cor <- "Network"
  OS <- gsub(" ", "-", Os)
  REMOVE <- ''
  if (r.i == "yes") {
    REMOVE <- paste(REMOVE, "Second degree interactions for proteins with unique link in first degree network\n", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-\n", sep = '')
  }
  if (r.s.i == "yes") {
    REMOVE <- paste(REMOVE, "Self interactions\n", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-\n", sep = '')
  }
  if (r.u.l == "yes") {
    REMOVE <- paste(REMOVE, "All unique links", sep = '')
  }
  else {
    REMOVE <- paste(REMOVE, "-", sep = '')
  }
  if (is.null(dim(not_founds)[1]) == T) {
    not_founds <- "All IDs are found in thesaurus"
  }
  
  
  cat('\n>Saving network ... ')
  
  # Sauvegarde du reseau
  out.put.name <- gfile('Save the network', type = "save", initial.dir = network.path, initial.filename = paste(ex.type, Cor, OS, inter.type, 'interactions.txt', sep = '_'))
  network2<-remove_redundants(network2)
  write.table(network2, file = out.put.name, row.names = F, col.names = T, quote = F, sep = "\t")
  
  setwd(network.path)
  
  cat('OK\nNetwork : ')
  cat(paste(out.put.name))
  cat('\n>Saving summary ... ')
  
  # Sauvegarde du recapitulatif des parametres utilises et des modifications apportees au reseau (et aux bases)
  out.put.name2 <- paste(OS, 'Summary_network.txt', sep = '_')
  write.table(paste("\n\n\n\n\nDATE :", Sys.time(), "NETWORK FILE NAME :", out.put.name, "ORGANISM :", os, "EXPERIENCE TYPE :", ex.type, "INTERACTIONS TYPE :", inter.type, "REMOVE :", REMOVE, "INPUTLIST :", IPL, "THESAURUS :", Th, "DATABASE(S) :", sep = '\n'), file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table(DB, file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table("\nFUNCTION USED :\nbuild_network_window()\nID(S) not found in thesaurus :\n", file = paste(network.path, out.put.name2, sep = '/'), append = T)
  write.table(not_founds, file = paste(network.path, out.put.name2, sep = '/'), append = T)
  
  cat('OK\nSummary : ')
  cat(paste(network.path, out.put.name2, sep = '/'))
  cat('\n\n>Construction of network is done.\n\n\n\n')
  
}

############################################################################################

remove_ppi <- function(protRang1, newPPI, inputlist) {
  
  # Supression des interactions non pertinentes
  # Les interactions avec des proteines de rang 2 sont supprimees si la proteine de rang 1 correspondante n'interagit qu'avec une seule proteine du reseau direct
  
  cat('\n>Removing second degree interactions for unique links in first degree interaction...')
  
  # Separation des interactions entre proteines differentes des autointeractions (boucles)
  newPPI <- as.matrix(newPPI)
  newPPI2 <- newPPI[newPPI[,1] == newPPI[,2],]
  if (dim(as.matrix(newPPI2))[2] == 1) {
    newPPI2 <- t(as.data.frame(newPPI2))
  }
  newPPI <- newPPI[newPPI[,1] != newPPI[,2],]
  # On recupere les proteines en interaction avec deux proteines de rang1.
  
  for(i in 1:length(protRang1)){	
    if(length(grep(protRang1[i],inputlist))>0){
      protRang1<-protRang1[-i]
    }
  }	
  

  # Recherche des proteines de rang1 n'ayant qu'une seule interaction dans le reseau direct
  SUP <- rep(FALSE, dim(newPPI)[1])
  supp <- FALSE
  for (m in 1:dim(newPPI)[1]) {
    n <- 0
    
    if((length(grep(newPPI[m,1],protRang1))>0)&&(length(grep(newPPI[m,2],protRang1))>0)) # si dans la colonne 1 et dans la colonne 2, il y a une prot de rang 1
    {	
      n<-2
    }
    

    if((length(grep(newPPI[m,1],inputlist))>0)&&(length(grep(newPPI[m,2],inputlist))>0)) # si dans la colonne 1 et dans la colonne 2, il y a une prot de rang 1
    {	
      n<-2
    }
    
    if((length(grep(newPPI[m,1],protRang1))>0)&&(length(grep(newPPI[m,2],inputlist))>0)) # si dans la colonne 1, il y a une prot de rang 1 et dans la colonne 2 une prot de l'inputlist
    {	
      listcol<-unique(c(as.vector(newPPI[newPPI[,1]==as.vector(newPPI[m,1]),2]),as.vector(newPPI[newPPI[,2]==as.vector(newPPI[m,1]),1])))
      protlist<-unique(c(protRang1,inputlist))
      for(l in 1:length(protlist))
      {
        n<-n+length(grep(protlist[l],listcol))
      }
    }
    if((length(grep(newPPI[m,2],protRang1))>0) && (length(grep(newPPI[m,1],inputlist))>0)) # si dans la colonne 2, il y a une prot de rang 1 et dans la colonne 1, une prot de l'inputlist 
    {
      listcol<-unique(c(as.vector(newPPI[newPPI[,1]==as.vector(newPPI[m,1]),2]),as.vector(newPPI[newPPI[,2]==as.vector(newPPI[m,1]),1])))
      protlist<-unique(c(protRang1,inputlist))
      for(l in 1:length(protlist))
      {
        n<-n+length(grep(protlist[l],listcol))
      }
      
    }				
    if(n<2)
    {
      SUP[m]<-TRUE
    }
    
  }

  newPPI <- newPPI[SUP == FALSE,]
  listProtNetwork <- unique(c(as.vector(newPPI[,1]), as.vector(newPPI[,2])))
  
  # Ajout des autointeractions, s'il y en a
  if (dim(newPPI2)[1] != 0) {
    for (i in 1:dim(newPPI2)[1]) {
      if (length(grep(newPPI2[i], listProtNetwork)) == 0) {
        newPPI2 <- newPPI2[-i,]
      }
    }
    newPPI2 <- rbind(newPPI2, newPPI)
  }
  else{
    newPPI2 <- newPPI
  }
  
  cat('OK')
  
  return(newPPI2)
  
}

############################################################################################

finish <- function(network, network2, r.s.i, r.u.l, UpDate, selected.database4, mainpath, network.path, Os, r.i, not_founds, ex.type, inter.type, os, IPL, Th, DB, organism.path) {
  # Memorisation des modifications faites sur le reseau pour mettre a jour les bases de donnees
  resume <<- as.matrix(data.frame(network[,4:5], network2[,4:5], stringAsFactors = F))
  resume <- resume[,-c(5)]
  colnames(resume) <- c("Old uidA", "Old uidB", "New uidA", "New uidB")
  # Supression des lignes pour lesquelles on ne modifie aucun identifiant
  duplicats <- c()
  for (i in 1:dim(resume)[2]) {
    if (resume[1,i] == resume[3,i] && resume[2,i] == resume[4,i]) {
      duplicats <- rbind(duplicats, i)
    }
  }
  if (length(duplicats) > 0) {
    resume <- resume[,-duplicats]
  }
  
  ### Toutes les corrections manuelles sont faites sur le reseau et memorisees ###
  
  
  # Second trie des interactions du reseau en fonction des parametre de construction selectionnes
  # Remove redundants
  cat('\n>Removing redundants\n')
  selected.database3 <- remove_redundants(network2)
  network2 <<- as.matrix(selected.database3)
  # Remove self-interactant
  if (r.s.i == 'yes') {
    cat('\n>Removing self-interactants ... ')
    network2 <<- network2[network2[,4] != network2[,5],]
    cat('OK')
  }
  else {
    cat('\n>Proteins which interact with itselves are kept')
  }
  # Remove unique links
  if (r.u.l == 'yes') {
    network2 <<- remove_unique_links(network2)
  }
  else {
    cat('\n>Proteins with only one interaction are kept')
  }
  
  ### Le dernier trie des interactions du reseau est effectue ###
  
  
  # Mise a jour des bases de donnees a partir des corrections effectuees sur le reseau (automatiques et manuelles)
  if (UpDate == 'yes') {
    update_db(selected.database4, resume, os, Os, organism.path)
  }
  
  ### Les bases de donnees sont corrigees en fonction des modifications apportees sur le reseau ###
  
  
  # Sauvegarde des resultats
  saving_window(network2, network.path, Os, r.i, r.s.i, r.u.l, not_founds, ex.type, inter.type, os, IPL, Th, DB)
  
  ### Le reseau, la matrice d'adjacence et le recapitulatif de construction sont sauvegardes dans des fichiers distincts ###
  
  setwd(mainpath)
  
}

update_db <- function(db, resume, os, Os, organism.path) {
  
  # Recuperation des bases et de la liste des modifications a approter
  database <- as.matrix(db)
  resume2 <- as.matrix(resume)
  
  # S'il y a eu des modifications dans le reseau on fait les memes dans les bases
  if (dim(resume2)[2] > 0) {
    cat('\n>Update')
    cat('\n')
    pb1 <- txtProgressBar(min = 0, max = dim(database)[1], style = 3)
    
    sources <- database[,11]
    sources <- unique(sources)
    
    # L'ensemble des bases est parcouru en entier pour la mise a jour et pour l'affectation a la base d'origine
    for (j in 1:dim(database)[1]) {
      # On corrige les lignes des bases qui sont dans la liste des corrections effectuees dans le reseau
      for (i in 1:dim(resume2)[2]) {
        if (resume2[1,i] == database[j,1] & resume2[2,i] == database[j,2]) {
          database[j,1] <- resume2[3,i]
          database[j,2] <- resume2[4,i]
        }
        if (resume2[1,i] == database[j,2] & resume2[2,i] == database[j,1]) {
          database[j,2] <- resume2[3,i]
          database[j,1] <- resume2[4,i]
        }
      }
      setTxtProgressBar(pb1, j)
    }
    
    # Creation d'un repertoire pour les bases mises a jour
    tax <- os
    setwd(organism.path)
    if (file.exists('Updated_databases') == F) dir.create('Updated_databases')
    updated.path <- paste(organism.path, '/Updated_databases', sep = '')
    setwd(updated.path)
    
    cat('\n\n>> Saving updated databases files in : ')
    cat(updated.path)
    cat(' directory')
    cat('\n>Saving ')
    cat(paste(length(sources)))
    cat(' updated database(s) :')
    
    # Sauvegarde des bases mises a jour
    OS <- gsub(" ",  "-", Os)
    if (length(sources) == 1) {
      This.base <- database
      out.put.name <- paste(OS, sources, "updated-database.txt", sep = "_")
      write.table(This.base[,-16], file = paste(organism.path, "Updated_databases", out.put.name, sep = "/"), append = F, row.names = F, quote = F, col.names = T, sep = "\t")
      
      cat('\n- Updated database : ')
      cat(paste(updated.path, out.put.name, sep = '/'))
    }
    else {
       for (s in 1:length(sources)) {
        This.base <- database[sources[s] == database[,11],]

        if((is.null(dim(This.base))==FALSE)){
        cat(paste("\nTaille de la base : "))
        cat(paste(dim(This.base)))
        out.put.name <- paste(OS, sources[s], "updated-database.txt", sep = "_")
        write.table(This.base[,-16], file = paste(organism.path, "Updated_databases", out.put.name, sep = "/"), append = F, row.names = F, quote = F, col.names = T, sep = "\t")
        
        cat('\n- Updated database ')
        cat(s)
        cat(' : ')
        cat(paste(organism.path, out.put.name, sep = '/'))
        }
      }
    }
  }
  
  else {
    cat('\n>No updates required')
  }
  
}

# fonction de recherche de l'uniprot ID de reference et association du bon nom de proteine et de gene
cherche_uniprotID<-function(data_i,data,thesaurus){
  # Colonne A
  resultatA <- search_id(data_i[4],thesaurus)
  if (length(resultatA) == 3) {
    ID <- resultatA[1]
    Proteine <- resultatA[2]
    Gene <- resultatA[3]
    data_i[4] <- ID
    data_i[1] <- Proteine
    data_i[11] <- Gene
  }
  v1<-c(grep(data_i[11],data[,11]))
  v2<-c(grep(data_i[11],data[,12]))
  if((length(unique(c(data[v1,4],data[v2,5])))>1)||(length(resultatA) != 3))
  {
    # On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
    not_found <- data_i[c(4,1,11)]
    nf<<-rbind(nf, not_found)
  }
  
  # Colonne B
  resultatB <- search_id(data_i[5], thesaurus)
  if (length(resultatB) == 3) {
    ID <- resultatB[1]
    Proteine <- resultatB[2]
    Gene <- resultatB[3]
    
    data_i[5] <- ID
    data_i[3] <- Proteine
    data_i[12] <- Gene
  }
  v3<-c(grep(data_i[12],data[,11]))
  v4<-c(grep(data_i[12],data[,12]))
  if((length(unique(c(data[v3,4],data[v4,5])))>1)||(length(resultatB) != 3))
  {
    # On recupere les identifiants que le thesaurus ne sait pas remplacer tout seul
    not_found <- data_i[c(5,3,12)]
    nf <<- rbind(nf, not_found)
  }
  nbpassage<<-nbpassage+1
  setTxtProgressBar(pb1, nbpassage)
}

