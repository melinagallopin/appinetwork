#!/usr/bin/python

import sys
import re

try:
	uniprot = sys.argv[sys.argv.index("thesaurusPy.py")+1]
except:    
	print ("ERROR: please, enter uniprot file")
	sys.exit()

try:
	thesaurus = sys.argv[sys.argv.index("thesaurusPy.py")+2]
except:    
	sys.exit()

try:
	organism = sys.argv[sys.argv.index("thesaurusPy.py")+3]
except:    
	sys.exit()


### Recuperation des donnees

fichier = open(uniprot,'r')

donnees = []
test = 0
ligne = 'NA'
liste = []
ID = 'NA'
Biogrid = 'NA'
BIOGRID = []
NameProt = 'NA'
Name = 'NA'
Ref = 'NA'
Syn = 'NA'
ORF = 'NA'
Iso = 'NA'
ISO = 'NA'
NIso = 'NA'

for i in fichier :

	i = i.strip("\n")

	# Recuperation du nom principal de la proteine (29, 30-37)
	if i[0:2] == "ID" :
		
		review_state = i[29:39]

		NameProt = re.search(r"(ID)(\s)*(?P<id>\w+)(_)", i)
		if NameProt is not None:
			NameProt = NameProt.group('id')
			NameProt = str(NameProt)
		else :
			NameProt = 'NA'

	if review_state == "Unreviewed" :
		continue

	# Recuperation de l'identifiant uniprot de la proteine et de ses anciens ID
	if i[0:2] == "AC" :
		if ID == "NA" :
			ID = re.search(r"(AC)(\s)*(?P<ac>\w+)(;)", i)
			if ID is not None:
				ID = ID.group('ac')
				ID = str(ID)
				line = str(i)
				line = line.strip("AC   ")
				liste = line.split()
				if len(liste) > 1 :
					Iso = ''
					for j in range(1, len(liste)) :
						Iso = Iso + liste[j]
			else :
				ID = 'NA'
		else :
			line = str(i)
			line = line.strip("AC   ")
			liste = line.split()
			if len(liste) > 0 :
				for j in range(0, len(liste)) :
					Iso = Iso + liste[j]



	# Recuperation de la reference du locus (NP)
	elif i[0:12] == "DR   RefSeq;" :
		if Ref == "NA" :
			Ref = re.search(r"(DR   RefSeq; NP_)(?P<rsq>\d+)(.)", i)
			if Ref is not None:
				Ref = Ref.group('rsq')
				Ref = "NP_" + str(Ref)
			else :
				Ref = re.search(r"(DR   RefSeq; NP_)(?P<rsq>\d+)(;)", i)
				if Ref is not None:
					Ref = Ref.group('rsq')
					Ref = "NP_" + str(Ref)
				else :
					Ref = re.search(r"(DR   RefSeq; XP_)(?P<rsq>\d+)(.)", i)
					if Ref is not None:
						Ref = Ref.group('rsq')
						Ref = "XP_" + str(Ref)
					else :
						Ref = re.search(r"(DR   RefSeq; XP_)(?P<rsq>\d+)(;)", i)
						if Ref is not None:
							Ref = Ref.group('rsq')
							Ref = "XP_" + str(Ref)
						else :
							Ref = 'NA'

	# Recuperation de l'identifiant biogrid
	elif i[0:13] == "DR   BioGrid;" :
		Biogrid = re.search(r"(DR   BioGrid;)(\s)*(?P<bg>\w+)(;)", i)
		if Biogrid is not None:
			Biogrid = Biogrid.group('bg')
			Biogrid = str(Biogrid)
			BIOGRID.append(Biogrid)
		else :
			Biogrid = 'NA'

	# Recuperation du nombre d'isoformes
	elif i[0:2] == "CC" :
		if NIso == 'NA' :
			NIso = re.search(r"(Named isoforms=)(?P<is>\d+)(;)", i)
			if NIso is not None :
				NIso = NIso.group('is')
				NIso = str(NIso)
			else :
				NIso = 'NA'

	# Recuperation de l'identifiant du gene : Hsapiens, Celegans, Dmelanogaster, Mmusculus, Rnorvegicus, Athaliana
	#elif i[0:13] == "DR   UniGene;" :
	#	if organism == "homo+sapiens" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(Hs.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'Hs.' + str(ORF)
	#			else :
	#				ORF = 'NA'	
	#	if organism == "caenorhabditis+elegans" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(Cel.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'Cel.' + str(ORF)
	#			else :
	#				ORF = 'NA'
	#	if organism == "drosophila+melanogaster" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(Dm.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'Dm.' + str(ORF)
	#			else :
	#				ORF = 'NA'	
	#	if organism == "mus+musculus" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(Mm.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'Mm.' + str(ORF)
	#			else :
	#				ORF = 'NA'
	#	if organism == "rattus+norvegicus" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(Rn.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'Rn.' + str(ORF)
	#			else :
	#				ORF = 'NA'	
	#	if organism == "arabidopsis+thaliana" :
	#		if ORF == 'NA' :
	#			ORF = re.search(r"(DR   UniGene;)(\s)*(At.)(?P<ug>\w+)(;)", i)
	#			if ORF is not None:
	#				ORF = ORF.group('ug')
	#				ORF = 'At.' + str(ORF)
	#			else :
	#				ORF = 'NA'		

	# Recuperation de l'identifiant du gene : Hsapiens, Celegans, Dmelanogaster, Mmusculus, Rnorvegicus, Athaliana + Ecoli, Scerevisiae
	elif i[0:10] == "DR   KEGG;" :
		if organism == "homo+sapiens" :
			if ORF == 'NA' :
				ORF = re.search(r"(hsa:)(?P<ug>\S+)(;)", i)
				if ORF is not None:
					ORF = ORF.group('ug')
					ORF = 'hsa:' + str(ORF)
				else :
					ORF = 'NA'
		if organism == "caenorhabditis+elegans" :
			if ORF == 'NA' :
				ORF = re.search(r"(cel:)(?P<ug>\S+)(;)", i)
				if ORF is not None:
					ORF = ORF.group('ug')
					ORF = 'cel:' + str(ORF)
				else :
					ORF = 'NA'
		if organism == "drosophila+melanogaster" :
			if ORF == 'NA' :
				ORF = re.search(r"(dme:)(?P<ug>\S+)(;)", i)
				if ORF is not None:
					ORF = ORF.group('ug')
					ORF = 'dme:' + str(ORF)
				else :
					ORF = 'NA'
		if organism == "mus+musculus" :
			if ORF == 'NA' :
				ORF = re.search(r"(mmu:)(?P<ug>\S+)(;)", i)
				if ORF is not None:
					ORF = ORF.group('ug')
					ORF = 'mmu:' + str(ORF)
				else :
					ORF = 'NA'
		if organism == "rattus+norvegicus" :
			if ORF == 'NA' :
				ORF = re.search(r"(rno:)(?P<ug>\S+)(;)", i)
				if ORF is not None:
					ORF = ORF.group('ug')
					ORF = 'rno:' + str(ORF)
				else :
					ORF = 'NA'
		#if organism == "escherichia+coli" :
		#	if ORF == 'NA' :
		#		ORF = re.search(r"(eco:)(?P<ug>\S+)(;)", i)
		#		if ORF is not None:
		#			ORF = ORF.group('ug')
		#			ORF = 'eco:' + str(ORF)
		#		else :
		#			ORF = 'NA'
		#if organism == "saccharomyces+cerevisiae" :
		#	if ORF == 'NA' :
		#		ORF = re.search(r"(sce:)(?P<ug>\S+)(;)", i)
		#		if ORF is not None:
		#			ORF = ORF.group('ug')
		#			ORF = 'sce:' + str(ORF)
		#		else :
		#			ORF = 'NA'	

	# Recuperationdu nom du gene, des synonymes et de l'identifiant du gene
	if i[0:2] == "GN" :

		# Recuperation des synonymes
		if test == 1 :
			ligneSyn = i.split(";")
			n = len(i)
			listeLigne = ligneSyn[0].strip("GN   ")
			laLigne = listeLigne.replace(" ", "")
			listeSyn = laLigne.split(",")
			Supprimer = []
			for k in range(0, len(listeSyn)) : # enlever {...}
				sup = re.search(r"({)(?P<supp>\S+)(})", listeSyn[k])
				if sup is not None :
					sup = sup.group('supp')
					sup = "{" + str(sup) + "}"
					listeSyn[k] = listeSyn[k].replace(sup, "")
				sup = re.search(r"({)(?P<supp>\S+)$", listeSyn[k])
				if sup is not None :
					sup = sup.group('supp')
					sup = "{" + str(sup)
					listeSyn[k] = listeSyn[k].replace(sup, "")
				sup = re.search(r"^(?P<supp>\S+)(})", listeSyn[k])
				if sup is not None :
					Supprimer.append(k)
			for s in Supprimer :
				del listeSyn[s]
			m = len(listeSyn)
			if i[n-1] == ";" and m > 0 : # tester si c'est la derniere ligne
				test = 0
				listeSyn[m - 1] = listeSyn[m - 1].replace(";", "")
			elif m == 0 :
				test = 0
				listeSyn = []
			else :
				del listeSyn[m - 1]
			if len(listeSyn) > 0 :
				ajout = ";".join(listeSyn)
				Syn = Syn + ajout + ";"
		SYN = re.search(r"(Synonyms=)", i)
		if SYN is not None :
			listeLigne = i.split(";")
			if len(listeLigne) == 2 : # la ligne ne contient que les synonymes, et la liste est complete sur une seule ligne
				listeLigne[0] = listeLigne[0].strip("GN   Synonyms=")
				laLigne = listeLigne[0].replace(" ", "")
				listeSyn = laLigne.split(",")
				Syn = ''
				Supprimer = []
				for k in range(0, len(listeSyn)) : # enlever {...}
					sup = re.search(r"({)(?P<supp>\S+)(})", listeSyn[k])
					if sup is not None :
						sup = sup.group('supp')
						sup = "{" + str(sup) + "}"
						listeSyn[k] = listeSyn[k].replace(sup, "")
					sup = re.search(r"({)(?P<supp>\S+)$", listeSyn[k])
					if sup is not None :
						sup = sup.group('supp')
						sup = "{" + str(sup)
						listeSyn[k] = listeSyn[k].replace(sup, "")
					sup = re.search(r"^(?P<supp>\S+)(})", listeSyn[k])
					if sup is not None :
						Supprimer.append(k)
				for s in Supprimer :
					del listeSyn[s]
				Syn = ";".join(listeSyn)
				Syn = Syn + ";"
			elif len(listeLigne) == 1 : # la ligne ne contient que les synonymes, et la liste ne se termine pas a la premiere ligne, elle est donc sur au moins deux lignes
				test = 1
				listeLigne[0] = listeLigne[0].strip("GN   Synonyms=")
				laLigne = listeLigne[0].replace(" ", "")
				listeSyn = laLigne.split(",")
				Syn = ''
				Supprimer = []
				for k in range(0, len(listeSyn)) : # enlever {...}
					sup = re.search(r"({)(?P<supp>\S+)(})", listeSyn[k])
					if sup is not None :
						sup = sup.group('supp')
						sup = "{" + str(sup) + "}"
						listeSyn[k] = listeSyn[k].replace(sup, "")
					sup = re.search(r"({)(?P<supp>\S+)$", listeSyn[k])
					if sup is not None :
						sup = sup.group('supp')
						sup = "{" + str(sup)
						listeSyn[k] = listeSyn[k].replace(sup, "")
					sup = re.search(r"^(?P<supp>\S+)(})", listeSyn[k])
					if sup is not None :
						Supprimer.append(k)
				for s in Supprimer :
					del listeSyn[s]
				del listeSyn[-1]
				Syn = ";".join(listeSyn)
				Syn = Syn + ";"
			else : # la ligne contient diverses informations en plus des synonymes, il faut donc retrouver le fragment concernant les synonymes
				for j in range(0, len(listeLigne)) :
					SYN = re.search(r"(Synonyms=)", listeLigne[j])
					if SYN is not None :
						if listeLigne[-1] == "" : # la liste se termine sur cette ligne
							if j == 0 :
								listeLigne[j] = listeLigne[j].strip("GN   Synonyms=")
							else :
								listeLigne[j] = listeLigne[j].strip(" Synonyms=")
							laLigne = listeLigne[j].replace(" ", "")
							listeSyn = laLigne.split(",")
							Syn = ''
							Supprimer = []
							for k in range(0, len(listeSyn)) : # enlever {...}
								sup = re.search(r"({)(?P<supp>\S+)(})", listeSyn[k])
								if sup is not None :
									sup = sup.group('supp')
									sup = "{" + str(sup) + "}"
									listeSyn[k] = listeSyn[k].replace(sup, "")
								sup = re.search(r"({)(?P<supp>\S+)$", listeSyn[k])
								if sup is not None :
									sup = sup.group('supp')
									sup = "{" + str(sup)
									listeSyn[k] = listeSyn[k].replace(sup, "")
								sup = re.search(r"^(?P<supp>\S+)(})", listeSyn[k])
								if sup is not None :
									Supprimer.append(k)
							for s in Supprimer :
								del listeSyn[s]
							Syn = ";".join(listeSyn)
							Syn = Syn + ";"
						else : # la liste est sur plusieurs lignes
							test = 1
							listeLigne[j] = listeLigne[j].strip(" Synonyms=")
							laLigne = listeLigne[j].replace(" ", "")
							listeSyn = laLigne.split(",")
							Syn = ''
							Supprimer = []
							for k in range(0, len(listeSyn)) : # enlever {...}
								sup = re.search(r"({)(?P<supp>\S+)(})", listeSyn[k])
								if sup is not None :
									sup = sup.group('supp')
									sup = "{" + str(sup) + "}"
									listeSyn[k] = listeSyn[k].replace(sup, "")
								sup = re.search(r"({)(?P<supp>\S+)$", listeSyn[k])
								if sup is not None :
									sup = sup.group('supp')
									sup = "{" + str(sup)
									listeSyn[k] = listeSyn[k].replace(sup, "")
								sup = re.search(r"^(?P<supp>\S+)(})", listeSyn[k])
								if sup is not None :
									Supprimer.append(k)
							for s in Supprimer :
								del listeSyn[s]
							del listeSyn[-1]
							Syn = ";".join(listeSyn)
							Syn = Syn + ";"
					else:
						SYN = 'NA'
		else :
			SYN = 'NA'

		# Recuperation de l'identifiant du gene : Scerevisiae + Ecoli
		if organism == "saccharomyces+cerevisiae" or organism == "escherichia+coli" :
			if ORF == "NA" :
				ligneORF = i.split(";")
				for k in ligneORF :
					k = k.replace(" ","")
					orfs = re.search(r"(.)*(OrderedLocusNames=)(?P<orf>\S+)$", k)
					if orfs is not None:
						ORF = orfs.group('orf')
						ORF = str(ORF)
						sup = re.search(r"({)(?P<SUP>\S+)(})", ORF)
						if sup is not None :
							sup = sup.group('SUP')
							sup = '{' + sup + '}'
							ORF = ORF.replace(sup, "")
					else :
						orfs = 'NA'
		if organism == "arabidopsis+thaliana" :
			if ORF == "NA" :
				ligneORF = i.split(";")
				for k in ligneORF :
					k = k.replace(" ","")
					k = k.split("{")[0]
					orfs = re.search(r"(.)*(OrderedLocusNames=)(?P<orf>\S+)$", k)
					if orfs is not None:
						ORF = orfs.group('orf')
						ORF = str(ORF)
						ORF = ORF.upper()
						sup = re.search(r"({)(?P<SUP>\S+)(})", ORF)
						if sup is not None :
							sup = sup.group('SUP')
							sup = '{' + sup + '}'
							ORF = ORF.replace(sup, "")
					else :
						orfs = 'NA'
		# Recuperation du nom du gene
		if Name == "NA" :
			ligneName = i.split(";")
			for k in ligneName :
				k = k.replace(" ", "")
				NAME = re.search(r"(Name=)(?P<id>\S+)$", k)
				if NAME is not None:
					Name = NAME.group('id')
					Name = str(Name)
					sup = re.search(r"({)(?P<SUP>\S+)$", Name)
					if sup is not None :
						sup = sup.group('SUP')
						sup = '{' + sup
						Name = Name.replace(sup, "")
				else :
					NAME = 'NA'

	# A la fin de chaque bloc (chaque proteine) les donnees sont rassemblees sur une meme ligne
	elif i[0:2] == "//" :
		if Name == 'NA' :
			Name = ORF
		if NameProt == 'NA' :
			NameProt = ID
		if Name == 'NA' :
			Name = NameProt
		if NIso == 'NA' :
			NIso = '0'
		b = len(BIOGRID)
		if b == 0 :
			Biogrid = 'NA'
		else:
			Biogrid = ''
			for B in range(0,b) :
				Biogrid = str(BIOGRID[B]) + ';'

		ligne = ID + "\t" + Biogrid + "\t" + Name + "\t" + Ref + "\t" + NameProt + "\t" + ORF + "\t" + Iso + "\t" + NIso
		donnees.append(ligne)

		ligne = 'NA'
		ID = 'NA'
		Biogrid = 'NA'
		BIOGRID = []
		Name = 'NA'
		NameProt = 'NA'
		Ref = 'NA'
		Syn = 'NA'
		ORF = 'NA'
		Iso = 'NA'
		ISO = 'NA'
		NIso = 'NA'
		test = 0

fichier.close()

# Recuperation des resultats dans un fichier de sortie

sortie = open(thesaurus,'w')

sortie.write("Uniprot-ID\tBiogrid-ID\tGene-Name\tRef-Seq\tProtein-Name\tGene-ID\tOld-Uniprot-IDs\tIsoformes")

N = len(donnees)
for i in range(0,N) :
	sortie.write("\n")
	sortie.write(donnees[i])
	
sortie.close()

