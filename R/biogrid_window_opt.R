rechercheDansThesaurus<-function(i,thesaurus_j,biogridA,biogridB,putative){
  uniprotthesaurus <- unlist(strsplit(thesaurus_j[2], ";"))
  for (u in 1:length(uniprotthesaurus)) {
  if(length(grep(biogridA,uniprotthesaurus[u]))>0 & length(grep(biogridB,uniprotthesaurus[u]))>0){
        
    # Recherche pour la proteine A
    if (biogridA == irefindex2[i,1] && biogridA == uniprotthesaurus[u]) {
      irefindex2[i,1] <- thesaurus_j[1]
      # Insertion des noms de proteine et de gene
      irefindex2[i,3] <- thesaurus_j[3]
      irefindex2[i,14] <- thesaurus_j[4]
      # Recherche des transposons
      if (thesaurus_j[6] == "TRUE") {
        identifiant <- paste("Transposon", irefindex2[i,1], sep = ":")
        irefindex2[i,1] <- gsub(irefindex2[i,1], identifiant, irefindex2[i,1])
        nomprot <- paste("Transposon", irefindex2[i,3], sep = ":")
        irefindex2[i,3] <- gsub(irefindex2[i,3], nomprot, irefindex2[i,3])
      }
    # Recherche des proteines putative
      if (putative == "FALSE") {
        if (thesaurus_j[5] == "TRUE") {
          identifiant <- paste("Putative", irefindex2[i,1], sep = ":")
          irefindex2[i,1] <- gsub(irefindex2[i,1], identifiant, irefindex2[i,1])
          nomprot <- paste("Putative", irefindex2[i,3], sep = ":")
          irefindex2[i,3] <- gsub(irefindex2[i,3], nomprot, irefindex2[i,3])
        }
      }
    }
    # Recherche pour la proteine B
    if (biogridB == irefindex2[i,2] && biogridB == uniprotthesaurus[u]) {
      irefindex2[i,2] <- thesaurus_j[1]
      # Insertion des noms de proteine et de gene
      irefindex2[i,4] <- thesaurus_j[3]
      irefindex2[i,15] <- thesaurus_j[4]
      # Recherche des transposons
      if (thesaurus_j[6] == "TRUE") {
        identifiant <- paste("Transposon", irefindex2[i,2], sep = ":")
        irefindex2[i,2] <- gsub(irefindex2[i,2], identifiant, irefindex2[i,2])
        nomprot <- paste("Transposon", irefindex2[i,4], sep = ":")
        irefindex2[i,4] <- gsub(irefindex2[i,4], nomprot, irefindex2[i,4])
      }
      # Recherche des proteines putative
      if (putative == "FALSE") {
        if (thesaurus_j[5] == "TRUE") {
          identifiant <- paste("Putative", irefindex2[i,2], sep = ":")
          irefindex2[i,2] <- gsub(irefindex2[i,2], identifiant, irefindex2[i,2])
          nomprot <- paste("Putative", irefindex2[i,4], sep = ":")
          irefindex2[i,4] <- gsub(irefindex2[i,4], nomprot, irefindex2[i,4])
        }
      }
    }
  } 
  }  
}


biogrid_othero <- function(organism, organismID, putative, file, uniprot, pathracine) {
  
  cat("\n\n>FORMATING BIOGRID DATABASE")
  
  # Creation du dossier correspondant a l'organisme
  organism <- gsub(" ", "-", organism)
  organisme <- gsub("-", "+", organism)
  organisme <- tolower(organisme)
  
  organism.path <- paste(pathracine, organism, sep = '/')
  dir.create(organism.path, showWarnings = FALSE)
  
  uniprotfile <- uniprot
  nomsortie <- paste("ThesaurusBiogrid", "_", organism, ".txt", sep = '')
  sortie <- paste(organism.path, nomsortie, sep = '/')
  
  # Construction du thesaurus
  cat("\n\n>Thesaurus construction ... ")
  path <- paste(system.file(package = "appinetwork"), "biogridPy.py", sep = "/")
  file.copy(path, pathracine)
  path.copy <- paste(pathracine, "biogridPy.py", sep = "/")
  
  command <- paste("python", "biogridPy.py", uniprotfile, sortie, putative, sep = " ")
  system(command)
  
  file.remove(path.copy)
  
  # Lecture du thesaurus
  thesaurus <- read.delim(file = sortie, header = T, sep = "\t")
  thesaurus <- as.matrix(thesaurus)
  file.remove(sortie)
  
  cat("OK\n>Loading database ... ")
  
  # Lecture de la base
  biogridfile <- file
  Biogrid <- read.delim(file = biogridfile, header = T, sep = "\t")
  
  if (dim(Biogrid)[2] != 24) {
    cat(paste("\nIncorrect data dimension, read userguide for more informations"))
    stop()
  }
  
  numParticipants <- rep(2, dim(Biogrid)[1])
  
  irefIndex <- cbind(Biogrid[,2], Biogrid[,3], Biogrid[,8], Biogrid[,9], Biogrid[,12], Biogrid[,14], Biogrid[,15], Biogrid[,16], Biogrid[,17], Biogrid[,13], Biogrid[,24], Biogrid[,19], numParticipants, Biogrid[,8], Biogrid[,9])
  colnames(irefIndex) <- c("uidA", "uidB", "aliasA", "aliasB", "method", "author", "pmids", "taxa", "taxb", "interactionType", "sourceDB", "confidence", "numParticipants", "GeneNameA", "GeneNameB")
  irefindex <- as.matrix(irefIndex)
  
  cat("OK")
  cat("\n>Formatting database : ")
  cat(paste(dim(Biogrid)[1]))
  cat(" lines in database\n")
  
  # Remplacement des Id biogrid par les uniprot et suppression des interactions non souhaitees
  prsbar <- txtProgressBar(min = 1, max = dim(Biogrid)[1], style = 3)
  N <- dim(irefindex)[1]
  for (i in 1:N) {
    # Verification de l'organisme
    organism <- gsub("-", " ", organism)
    if (organismID == irefindex[i,8] && organismID == irefindex[i,9]) {
      irefindex[i,8] <- paste("taxid:", organismID, "(", organism, ")", sep = '')
      irefindex[i,9] <- paste("taxid:", organismID, "(", organism, ")", sep = '')
      irefindex[i,7] <- paste("pubmed:", irefindex[i,7], sep = '')
    }
    
    irefindex[i,5] <- paste("", as.character(Biogrid[i,12]), sep = '')
    irefindex[i,6] <- paste("", as.character(Biogrid[i,14]), sep = '')
    irefindex[i,10] <- paste("", as.character(Biogrid[i,13]), sep = '')
    irefindex[i,11] <- paste("", as.character(Biogrid[i,24]), sep = '')
    irefindex[i,12] <- paste("", as.character(Biogrid[i,19]), sep = '')
    
    # Remplacement de l'identifiant
    biogridA <- irefindex[i,1]
    biogridB <- irefindex[i,2]
    irefindex2<<-irefindex
    apply(thesaurus,1,rechercheDansThesaurus,i=i,biogridA=biogridA,biogridB=biogridB,putative=putative)
    irefindex<-irefindex2
    rm(irefindex2)
    if (irefindex[i,1] == biogridA) {
      irefindex[i,1] <- paste("Protein", biogridA, sep = ":")
      irefindex[i,3] <- paste("Protein", irefindex[i,3], sep = ":")
    }
    if (irefindex[i,2] == biogridB) {
      irefindex[i,2] <- paste("Protein", biogridB, sep = ":")
      irefindex[i,4] <- paste("Protein", irefindex[i,4], sep = ":")
    }
    setTxtProgressBar(prsbar, i)
  }
  
  cat("\n>Saving database ... ")
  
  setwd(organism.path)
  database.path <- paste(organism.path, "Databases", sep = '/')
  dir.create(database.path, showWarnings = FALSE)
  setwd(database.path)
  
  organism <- gsub(" ", "-", organism)
  outfilename <- paste(organism, "_biogrid.txt", sep = "")
  write.table(irefindex, file = outfilename, row.names = F, col.names = T, quote = F, sep = "\t")
  
  cat(paste("OK\n\n>Formating biogrid database is done.\n\nDatabase file is saved in", database.path, sep = " : "))
  cat(paste("\n\n"))
  setwd(pathracine)
  
  visible(mainpanel) <- T
  
}

biogrid_window2 <- function(f_pos, mainpanel, pana, mainpath) {
    file <- c()
    uniprot <- c()
    return.parameter <- c()

    panel_para <- gwindow("Biogrid database file formatting :", parent = f_pos, visible = T)
    pp <- gvbox(container = panel_para)

    flyt <- gformlayout(container = pp, expand = TRUE)

    # Selection des parametres
    gcombobox(c('Caenorhabditis elegans', 'Drosophila melanogaster', 'Escherichia coli', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Other'), label = "Organism", selected = 7, container = flyt)
    gradio(c("Yes", "No"), selected = 1, horizontal = TRUE, label = "Conserve interactions with putative proteins", container = flyt)

    # Selection du fichier
    chdb <- ggroup(container = pp, horizontale = T)
    addSpring(chdb)
    bouton1 <- gbutton("Biogrid file", container = chdb, handler = function(...) {
        file <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
        if (is.null(file) == T) {
            gmessage('Selected Biogrid file is null', icon = 'error')
        }
        if (is.null(file) == F) {
            bouton1$set_value(paste(length(file), 'biogrid file selected'))
            cat(paste('\n>Biogrid file selected : ', file, sep = ''))
        }
    })
    bouton2 <- gbutton("Uniprot file", container = chdb, handler = function(...) {
        uniprot <<- gfile(text = "Select a file", type = "open", multi = F, container = chdb)
        if (is.null(uniprot) == T) {
            gmessage('Selected Uniprot file is null', icon = 'error')
        }
        if (is.null(uniprot) == F) {
            bouton2$set_value(paste(length(uniprot), 'uniprot file selected'))
            cat(paste('\n>Uniprot file selected : ', uniprot, sep = ''))
        }
    })

    ppb <- ggroup(container = pp)
    addSpring(ppb)
    gbutton("Formate file", handler = function(h,...) {

        return.parameter <<- svalue(flyt)
        visible(panel_para) <- F
        # Memorisation des parametres de recherche selectionnes
        organism <- as.character(return.parameter[1])
        putative <- as.character(return.parameter[2])


        # Execution de la fonction lorsque tous les parametres sont fournis
        if (is.null(file) == F) {
            ######################################################################################
            ########################### Execution de la fonction #################################

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

                # 1 : Autre organism
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

                    biogrid_othero(organism, organismID, putative, file, uniprot, pathracine)

                    dispose(bpc)
                }, container = bpc)

            }

            else {

                cat("\n\n>FORMATING BIOGRID DATABASE")

                # Creation du dossier correspondant a l'organisme
                organism <- gsub(" ", "-", organism)
                organisme <- gsub("-", "+", organism)
                organisme <- tolower(organisme)

                organism.path <- paste(pathracine, organism, sep = '/')
                dir.create(organism.path, showWarnings = FALSE)

                uniprotfile <- uniprot
                nomsortie <- paste("ThesaurusBiogrid", "_", organism, ".txt", sep = '')
                sortie <- paste(organism.path, nomsortie, sep = '/')

                # Construction du thesaurus
                cat("\n\n>Thesaurus construction ... ")
                path <- paste(system.file(package = "appinetwork"), "biogridPy.py", sep = "/")
                file.copy(path, pathracine)
                path.copy <- paste(pathracine, "biogridPy.py", sep = "/")
                command <- paste("python", "biogridPy.py", uniprotfile, sortie, putative, sep = " ")
                system(command)

                file.remove(path.copy)

                # Lecture du thesaurus
                thesaurus <- read.delim(file = sortie, header = T, sep = "\t")
                thesaurus <- as.matrix(thesaurus)
                #file.remove(sortie)

                cat("OK\n>Loading database ... ")

                # Lecture de la base
                biogridfile <- file
                Biogrid <- read.delim(file = biogridfile, header = T, sep = "\t")


                if (dim(Biogrid)[2] != 24) {
                    cat(paste("\nIncorrect data dimension, read userguide for more informations"))
                    stop()
                }
                Biogrid<-as.matrix(Biogrid)
                numParticipants <- rep(2, dim(Biogrid)[1])
                irefIndex <- cbind(Biogrid[,2], Biogrid[,3], Biogrid[,8], Biogrid[,9], Biogrid[,12], Biogrid[,14], Biogrid[,15], Biogrid[,16], Biogrid[,17], Biogrid[,13], Biogrid[,24], Biogrid[,19], numParticipants, Biogrid[,8], Biogrid[,9])
                colnames(irefIndex) <- c("uidA", "uidB", "aliasA", "aliasB", "method", "author", "pmids", "taxa", "taxb", "interactionType", "sourceDB", "confidence", "numParticipants", "GeneNameA", "GeneNameB")
                irefindex <- as.matrix(irefIndex)

                organism <- gsub("-", " ", organism)
                list.organism<-cbind(c("Saccharomyces cerevisiae","Homo sapiens","Caenorhabditis elegans","Drosophila melanogaster","Escherichia coli","Mus musculus","Rattus norvegicus"),c("559292","9606","6239","7227","562","10090","10116"))
                organismID<-list.organism[list.organism[,1]==organism,2]

                cat("\n>Formatting database : ")
                cat(paste(dim(Biogrid)[1]))
                cat(" lines in database\n")

                prsbar <- txtProgressBar(min = 1, max = dim(Biogrid)[1], style = 3)


                irefindex[,5] <- paste("", as.character(Biogrid[,12]), sep = '')
                irefindex[,6] <- paste("", as.character(Biogrid[,14]), sep = '')
                irefindex[,10] <- paste("", as.character(Biogrid[,13]), sep = '')
                irefindex[,11] <- paste("", as.character(Biogrid[,24]), sep = '')
                irefindex[,12] <- paste("", as.character(Biogrid[,19]), sep = '')
                irefindex<-irefindex[organismID == irefindex[,8] && organismID == irefindex[,9],]
                irefindex[,8] <- paste("taxid:", organismID, "(", organism, ")", sep = '')
                irefindex[,9] <- paste("taxid:", organismID, "(", organism, ")", sep = '')
                irefindex[,7] <- paste("pubmed:", irefindex[,7], sep = '')

                # eliminate transposons
                thesaurusNotranNoPut<-thesaurus[thesaurus[,6]=="FALSE",]
                  
                if(length(thesaurus[thesaurus[,6]==" TRUE",6])>0){
                  Tab.transposons<-thesaurus[thesaurus[,6]==" TRUE",]
                  Tab.transposons[,2]<-unlist(apply(Tab.transposons,c(1,2), strsplit,split=";")[,2])
                  for (i in 1:dim(Tab.transposons)[1]){
                    if(length(grep(Tab.transposons[i,2],irefindex[,1]))>0)
                      irefindex<-irefindex[-grep(Tab.transposons[i,2],irefindex[,1]),]
                    if(length(grep(Tab.transposons[i,2],irefindex[,2]))>0)
                      irefindex<-irefindex[-grep(Tab.transposons[i,2],irefindex[,2]),]
                  }
                }
                # eliminate putative proteins
                if (putative == " TRUE") {
                    thesaurusNotranNoPut<-thesaurusNotranNoPut[thesaurusNotranNoPut[,5]=="FALSE",]
                    Tab.putatives<-thesaurus[thesaurus[,5]==" TRUE",]
                    Tab.putatives[,2]<-unlist(apply(Tab.putatives,c(1,2), strsplit,split=";")[,2])
                    for (i in 1:dim(Tab.putatives)[1]) {
                        if(length(grep(Tab.putatives[i,2],irefindex[,1]))>0)
                        irefindex<-irefindex[-grep(Tab.putatives[i,2],irefindex[,1]),]
                        if(length(grep(Tab.putatives[i,2],irefindex[,2]))>0)
                        irefindex<-irefindex[-grep(Tab.putatives[i,2],irefindex[,2]),]
                    }
                }
                uniprotthesaurusNotranNoPut<-unlist(apply(thesaurusNotranNoPut,c(1,2), strsplit,split=";")[,2]) # new
                biogridA<-irefindex[,1]
                biogridB<-irefindex[,2]

                for(i in 1:dim(thesaurusNotranNoPut)[1]){
                    for (u in 1:length(uniprotthesaurusNotranNoPut[i])){
                        irefindex[irefindex[,1]==uniprotthesaurusNotranNoPut[i][u],3]<-thesaurusNotranNoPut[i,3]
                        irefindex[irefindex[,1]==uniprotthesaurusNotranNoPut[i][u],14]<-thesaurusNotranNoPut[i,4]
                        irefindex[irefindex[,1]==uniprotthesaurusNotranNoPut[i][u],1]<-thesaurusNotranNoPut[i,1]

                        irefindex[irefindex[,2]==uniprotthesaurusNotranNoPut[i][u],4]<-thesaurusNotranNoPut[i,3]
                        irefindex[irefindex[,2]==uniprotthesaurusNotranNoPut[i][u],15]<-thesaurusNotranNoPut[i,4]
                        irefindex[irefindex[,2]==uniprotthesaurusNotranNoPut[i][u],2]<-thesaurusNotranNoPut[i,1]
                    }
                }
                #irefindex<-irefindex[(irefindex[,1] != biogridA) &(irefindex[,2] != biogridB),]
                cat("\n>Saving database ... ")

                setwd(organism.path)
                database.path <- paste(organism.path, "Databases", sep = '/')
                dir.create(database.path, showWarnings = FALSE)
                setwd(database.path)

                organism <- gsub(" ", "-", organism)
                outfilename <- paste(organism, "_biogrid.txt",sep="")
                write.table(irefindex, file = outfilename, row.names = F, col.names = T, quote = F, sep = "\t")

                cat(paste("OK\n\n>Formating biogrid database is done.\n\nDatabase file is saved in", database.path, sep = " : "))
                cat(paste("\n\n"))
                setwd(pathracine)

                visible(mainpanel) <<- T

            }


            dispose(panel_para)
            dispose(ppb)

        }

        # Tests de la presence de tous les elements necessaires a l'execution de la fonction
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

	# Retour
	gbutton("Return", handler = function(h,...) {
		dispose(panel_para)
		visible(pana) <<- T
	}, container = ppb)

	visible(panel_para) <- T

}
