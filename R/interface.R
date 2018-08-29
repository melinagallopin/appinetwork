interface <- function() {

	mainpath <- getwd()
	f_pos <- c(1000, 500)

	# Main navigation window allowing the access to the different package's options
	mainpanel <<- gwindow("PPI interface" , parent = f_pos, visible = F)
	g <- ggroup(container = mainpanel, horizontal = F)

    # Construction of an ID correspondence file for a selected organism
	bouton0 <- gbutton("Construct an ID corespondence file", container = g, handler = function(...) {
		visible(mainpanel) <- F
		thesaurus_window(f_pos, mainpanel, mainpath)

	})

	# Formatting of database files used for the PPI network building
	bouton1 <- gbutton("Format data files", container = g, handler = function(...) {
		visible(mainpanel) <- F
		# Display possible options
		pana <<- gwindow("Databases formatting", parent = f_pos, width = '300', height = '300', visible = F)
		pag <- ggroup(container = pana, horizontal = F)

		# Irefindex database
		bouton1 <- gbutton("IrefIndex", container = pag, handler = function(...) {
			visible(pana) <- F
			irefindex_window(f_pos, mainpanel, pana, mainpath)
		})

		# Intact database
		bouton2 <- gbutton("Intact", container = pag, handler = function(...) {
			visible(pana) <- F
			intact_window(f_pos, mainpanel, pana, mainpath)
		})

		# Biogrid database
		bouton3 <- gbutton("Biogrid", container = pag, handler = function(...) {
		  visible(pana) <- F
		  biogrid_window2(f_pos, mainpanel, pana, mainpath)
		})

		# Return to the main navigation window (principal menu)
		gbutton("Return", handler = function(h,...) {
			dispose(pana)
			visible(mainpanel) <- T
			}, container = pag)

		visible(pana) <- T

	})

	# Building the network
	bouton3 <- gbutton("Building the network", container = g, handler = function(...) {
		visible(mainpanel) <- F
		build_network_window(f_pos, mainpanel, mainpath)
	})

	# Modeling of the assembly intermediaries of a protein complex
	bouton4 <- gbutton("Modeling the assembly intermediaries of a protein complex", container = g, handler = function(...) {
		visible(mainpanel) <- F
		assembly_intermediary_window(f_pos, mainpanel)
	})

	# PPI network clustering to identify proteins of a biological process
	bouton5 <- gbutton("Identify proteins of a biological process", container = g, handler = function(...) {
		visible(mainpanel) <- F
		# Diplay Clustering
		panb <<- gwindow("Clustering", parent = f_pos, width = '300', height = '300', visible = F)
		pagb <- ggroup(container = panb, horizontal = F)

		# Clustering with TFit
		bouton1 <- gbutton("TFit", container = pagb, handler = function(...) {
		  visible(panb) <- F
		  clustering_tfit_window(f_pos, mainpanel, panb, mainpath)
		})

		# Return to the main navigation window (principal menu)
		gbutton("Return", handler = function(h,...) {
			dispose(panb)
			visible(mainpanel) <- T
			}, container = pagb)

		visible(panb) <- T

	})

	visible(mainpanel) <- T

}

