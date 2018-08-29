#!/usr/bin/python

import sys
import re

### python biogridPy.py uniprotfile sortie putative

try:
	uniprot = sys.argv[sys.argv.index("biogridPy.py")+1]
except:    
	print ("ERROR: please, enter uniprot file")
	sys.exit()

try:
	thesaurus = sys.argv[sys.argv.index("biogridPy.py")+2]
except:    
	sys.exit()


### Recuperation des donnees

fichier = open(uniprot,'r')

donnees = []
ligne = 'NA'
ID = 'NA'
NameProt = 'NA'
Name = 'NA'
Biogrid = 'NA'
BIOGRID = []
Putative = "FALSE"
Transposon = "FALSE"

for i in fichier :

	i = i.strip("\n")

	# Recuperation de l'identifiant uniprot de la proteine
	if i[0:2] == "AC" :
		if ID == "NA" :
			ID = re.search(r"(AC)(\s)*(?P<ac>\w+)(;)", i)
			if ID is not None:
				ID = ID.group('ac')
				ID = str(ID)
			else :
				ID = 'NA'

	# Recuperation du nom principal de la proteine
	if i[0:2] == "ID" :
		NameProt = re.search(r"(ID)(\s)*(?P<id>\w+)(_)", i)
		if NameProt is not None:
			NameProt = NameProt.group('id')
			NameProt = str(NameProt)
		else :
			NameProt = 'NA'
			
	# Recuperation du nom du gene
	if i[0:2] == "GN" :
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

	# Recuperation de l'identifiant biogrid
	elif i[0:12] == "DR   GeneID;" :
		Biogrid = re.search(r"(DR   GeneID;)(\s)*(?P<bg>\w+)(;)", i)
		if Biogrid is not None:
			Biogrid = Biogrid.group('bg')
			Biogrid = str(Biogrid)
			BIOGRID.append(Biogrid)
		else :
			Biogrid = 'NA'

	# Recuperation de l'identifiant biogrid : id alternatif
	elif i[0:12] == "DR   BioGrid;" :
		Biogrid = re.search(r"(DR   GeneID;)(\s)*(?P<bg>\w+)(;)", i)
		if Biogrid is not None:
			Biogrid = Biogrid.group('bg')
			Biogrid = str(Biogrid)
			BIOGRID.append(Biogrid)
		else :
			Biogrid = 'NA'

	# Recuperation du RecName : Putative, Transposon
	elif i[0:13] == "DE   RecName:" :
		if Putative == "FALSE" :
			Putative = re.search(r"(DE   RecName:)(\s)*(Full=Putative)(.)*", i)
			if Putative is not None:
				Putative = "TRUE"
			else :
				Putative = "FALSE"
		if Transposon == "FALSE" :
			Transposon = re.search(r"(DE   RecName:)(\s)*(Full=Transposon)(.)*", i)
			if Transposon is not None:
				Transposon = "TRUE"
			else :
				Transposon = "FALSE"

	# A la fin de chaque bloc (chaque proteine) les donnees sont rassemblees sur une meme ligne
	elif i[0:2] == "//" :
		b = len(BIOGRID)
		if b > 0 : 
			Biogrid = ''
			for B in range(0,b) :
				Biogrid = str(BIOGRID[B]) + ';'
			ligne = ID + "\t" + Biogrid + "\t" + NameProt + "\t" + Name + "\t" + Putative + "\t" + Transposon
			donnees.append(ligne)

		ligne = 'NA'
		ID = 'NA'
		NameProt = 'NA'
		Name = 'NA'
		Biogrid = 'NA'
		BIOGRID = []
		Putative = "FALSE"
		Transposon = "FALSE"

fichier.close()

# Recuperation des resultats dans un fichier de sortie

sortie = open(thesaurus,'w')

sortie.write("Uniprot-ID\tBiogrid-ID\tProtein-Name\tGene-Name\tPutative\tTransposon")

N = len(donnees)
for i in range(0,N) :
	sortie.write("\n")
	sortie.write(donnees[i])
	
sortie.close()

