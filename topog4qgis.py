# -*- coding: utf-8 -*-

"""
/***************************************************************************
topog4qgis		A QGIS plugin Tools for managing Topographic tool on vector
				layers

                             -----------------------------
        begin                : 2013-10-30
        copyright            : (C) 2013 by Giuliano Curti (orinal author)
        email                : giulianc51@gmail.com
        
        updated on           : 2020-04-29
        maintainer           : Marco Lombardi
        email                : marco.lombardi.rm@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

# Import standard libraries
import os,math,copy,functools,operator
# Import codecs library
import codecs
# Import the PyQt and QGIS libraries
from PyQt5 import QtCore, QtWidgets, QtGui
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *
from qgis.core import *
from qgis.gui import *

# ----- help functions -----------

def about(mw,parent):
	"""
		Visualizza info sulla procedura
		presume che parent abbia le variabili: vers,build_date,author,mail,copyright,license
	"""
	QMessageBox.about(
		mw,
		'About',
		"TOPOGRAPHIC tools for qGis"
+ "\n----------------------------"
+ "\nversion:      %s" % (parent.vers)
+ "\nbuild_date:   %s" % (parent.build_date)
+ "\nauthor:       %s" % (parent.author)
+ "\ncontributor:  %s" % (parent.contributor)
+ "\nmaintainer:   %s" % (parent.maintainer)
+ "\ncopyright:    %s" % (parent.copyright)
+ "\nlicense:      %s" % (parent.license)
	)

def info(mw,vers,mail):
	"""
		First advice to users
	"""
	msg = """

	"""

	QMessageBox.about(
		mw,
		'Info',
		"TOPOGRAPHIC tools for QGIS ("
+ vers
+ ")"
+ "\n-----------------------------------------------"
+ "\nThe procedure allows the survey management;"
+ "\nnow this is limited to reading of the site"
+ "\nrecords (Libretto di campagna), in the forms."
+ "\nof optical, gps and alignment records."
+ "\n"
+ "\nReadings of EdM, TAF and DIS files of the"
+ "\nItalian Cadastral Service are allowed also."
+ "\n"
+ "\nThe user provided data are matched to compose"
+ "\nthe topographical shape of the area of interest"
+ "\nand displayed as a series of vector layers."
+ "\n"
+ "\nThe availability of the terrain control points"
+ "\n(in italian: PF punti fiduciali) allows the"
+ "\ntrasformation of the land survey to correct"
+ "\nCRS."
+ "\n"
+ "\nThe user can investigate their survey by"
+ "\ninquiring all information on stations, vertices,"
+ "\nPF, etc. by plugin functions in addition to"
+ "\nthe QGIS native tools, as in example to"
+ "\nsave them in different formats too."
+ "\n"
+ "\nOpening the python console you can enjoy some"
+ "\nuseful info from the system."
+ "\n" 
+ "\nThis procedure is EXPERIMENTAL; it might"
+ "\ncontains many bugs, few duplications and some"
+ "\nmistakes, only in part known to the author;"
+ "\n"
+ "\nplease let us know about any encountered"
+ "\nproblems"
	)

# ----- graphic functions -----------

def searchFeat(parent,point):
	"""
		Seleziona le features individuate con il mouse;
		presume che il parent abbia la variabile: cLayer.eps (raggio di ricerca) 
	"""
	# setup the provider select to filter results based on a rectangle
	pntGeom = QgsGeometry.fromPointXY(point)  
	# scale-dependent buffer of 5 pixels-worth of map units
	pntBuff = pntGeom.buffer(parent.eps,0) 
	rect = pntBuff.boundingBox()
	# create the select statement
	parent.cLayer.select(rect,True)	# prende quelli che intersecano

def pointGetCds(feat):
	"""
		Restituisce id e coordinate di una feature PUNTO;
		in caso di fallimento restituisce un id = -1.
	"""
	fid = feat.id()
	geom = feat.geometry()
	if geom.type() == QGis.Point:
		cod = feat.attributes()[0]
		pnt = geom.asPoint()
		x,y = pnt.x(),pnt.y()
		z = feat.attributes()[1]
		note = feat.attributes()[2]
		stazione = feat.attributes()[3]
		nRiga = feat.attributes()[4]
	else:
		fid,cod,x,y,z,note,nRiga = -1,0,0.0,0.0,0.0,'',-1
	return fid,cod,x,y,z,note,stazione,nRiga

def singleSymbol(geoTyp,props):
	"""
		prepara la simbologia singola
		- geoTyp	type of geometry
		- props properties of symbol
	"""
	namCol = props
	symb = QgsSymbol.defaultSymbol(geoTyp)
	col = symb.color()
	col.setNamedColor(namCol)
	symb.setColor(col)
	# crea il renderer
	render = QgsSingleSymbolRenderer(symb)
	return render

def catSymbol(geoTyp,attName,cats):
	"""
		prepara la simbologia categorizzata
		- geoTyp	type of geometry
		- attName	nome dell'attributo
		- cats lista delle categorie [catName,namCol,label]
	"""
	# crea la lista delle categorie
	catList = []
	# parsing delle categorie
	for i in cats:
		name,namCol,lab = i
		if name == 'NT':
			symb = QgsLineSymbol.createSimple({'line_style': 'dash', 'color': 'black'})
		elif name == 'RT':
			symb = QgsLineSymbol.createSimple({'line_style': 'dash', 'color': 'red'})
		elif name == 'RP':
			symb = QgsLineSymbol.createSimple({'line_style': 'dot', 'color': 'red'}) 
		elif name == 'NP':
			symb = QgsLineSymbol.createSimple({'line_style': 'dot', 'color': 'black'})             
		else:
			symb = QgsSymbol.defaultSymbol(geoTyp)
		col = symb.color()
		col.setNamedColor(namCol)
		symb.setColor(col)
		cat = QgsRendererCategory(name,symb,lab)
		catList.append(cat)
	# crea il renderer
	render = QgsCategorizedSymbolRenderer(attName,catList)
	return render

def annotationText(canvas,text,pos):
	"""
		disegna una didascalia a video (purtroppo c'è uno spostamento incomprensibile)
		- canvas	area di disegno
		- text	testoda visualizzare
		- pos		posizione (QgsPoint)
	"""
	myDoc = QTextDocument(text)
	myTxt = QgsTextAnnotation(canvas)
	myTxt.setDocument(myDoc)
	myTxt.setMapPosition(pos) 
	canvas.refresh()
	return myTxt

def drawStar(parent,center,vertices,archivio):
	"""
		disegna una stella con il rubberBand dal center ai vertices
		questi indicano la posizione dei records che contengomo le coordinate
	"""
	l,x,y,z = pointArchivioCds(archivio,center)
	base = QgsPoint(x,y)
	for p in vertices:
		l,x,y,z = pointArchivioCds(archivio,p)
		pnt = QgsPoint(x,y)
		# attiva il RB
		parent.rubBnd.addPoint(base)
		parent.rubBnd.addPoint(pnt)
	QMessageBox.information(
		parent.iface.mainWindow(),
		"star",
		"Controlla"
	)
	# rimuove la rubber band
	parent.rubBnd.reset()

# ----- trigonometric functions -----------

def toRad(vIn,a_giro):
	"""
		converte gradi sessadecimali dalla notazione corrente (a_giro) in radianti
	"""
	return (vIn/a_giro)*2*math.pi

def fromRad(angRad,a_giro):
	"""
		Converte i radianti in gradi sessadecimali nella notazione corrente (a_giro)
	"""
	a = angRad * a_giro / (2*math.pi)
	if a < 0:
		a = a + a_giro
	return a

# ----- matrix functions -----------

def printMatrix(mat):
	nr = len(mat)
	nc = len(mat[0])
	print(' ')
	for c in range(nc):
		print('%5s ' % (c+1))
	print('------' * (nc+2))
	for r in range(nr):
		print(r+1, '|')
		for c in range(nc):
			print('%5.2f ' % mat[r][c])
	print('------' * (nc+2))

def nullMatrix(nr,nc):
	"""
		resituisce la matrice nulla di nr x nc
	"""
	mat = []
	for i in range(nr):
		tmp = []
		for j in range(nc):
			tmp.append(0.0)
		mat.append(tmp)
	return mat

def identity(n):
	"""
		Restituisce la matrice identità
		di dimensione n
	"""
	mat = nullMatrix(n,n)
	for i in range (0,n):
		mat[i][i] = 1.0
	return mat

def matrixTranspose(mat):
	"""
		esegue la trasposta di una matrice
	"""
	nr = len(mat)
	nc = len(mat[0])
	matT = nullMatrix(nc,nr)
	for r in range(nr):
		for c in range(nc):
			matT[c][r] = mat[r][c]
	return matT

def adjoint(mat1,nr,mat2):
	"""
		Costruisce la matrice aggiunta composta dalle due matrici in input
		(presume le matrice di pri numero di righe)
	"""
	nc2 = len(mat2[1])
	adj = []
	for i in range(nr):
		tmp = mat1[i][:]
		for j in mat2[i]:
			tmp.append(j)
		adj.append(tmp)
	return adj

def matTranslation3D(dx,dy,dz):
	"""
		Calcola la matrice di traslazione;
	"""
	return [
		[1.,0.,0.,dx],
		[0.,1.,0.,dy],
		[0.,0.,1.,dz],
		[0.,0.,0.,1.]
	]

def matRotation3DinZ(a):
	"""
		rotazione dell'angolo a [radians]
	"""
	c= math.cos(a)
	s = math.sin(a)
	return [
		[c,-s,0.,0.],
		[s, c,0.,0.],
		[0.,0.,1.,0.],
		[0.,0.,0.,1.]
	]

def matrixMultiplication(mat1,mat2):
	"""
		Calcola il prodotto delle matrici mat1*mat2
	"""
#	print "matrice 1:"
#	self.printMatrix(mat1)
#	print "matrice 2:"
#	self.printMatrix(mat2)
	nr1 = len(mat1)
	nc1 = len(mat1[0])
	nr2 = len(mat2)
	nc2 = len(mat2[0])
	mat = []
	if nc1 == nr2:
		for i in range(nr1):
			tmp = []
			for j in range (nc2):
				val = 0
				for k in range (nc1):
			  		val += mat1[i][k] * mat2[k][j]
				tmp.append(val)
			mat.append(tmp)
		return mat
	else:
		print ('matrici non congruenti per la moltiplicazione')
		return -1

def matRotoTrasla(s1,s2,d1,d2):
	"""
		matrice di collimazione a 2 punti
		s1	vecchio centro
		s2	vecchio punto di allineamento
		d1	nuovo centro
		d2	nuovo punto di allineamento
		rototrasla in XY e trasla in Z
	"""
	s1x,s1y,s1z = s1
	s2x,s2y,s2z = s2
	d1x,d1y,d1z = d1
	d2x,d2y,d2z = d2
	# calcola matrice di trasformazione
	mat = matTranslation3D(-s1x,-s1y,-s1z)
#	print "dopo 1.a traslazione"
#	printMatrix(mat)
	# ruota la linea s1-s2 sull'orizzontale
	a1 = math.atan2(s2y-s1y,s2x-s1x)
	# ruota per allinearla a d1-d2
	a2 = math.atan2(d2y-d1y,d2x-d1x)
#	print "ruota di",a2-a1
	mat1 = matRotation3DinZ(a2-a1)
	mat = matrixMultiplication(mat1,mat)
	#print("dopo rotazione")
	#printMatrix(mat)
	#print("trasla nel punto d1",d1x,d1y,d1z)
	mat1 = matTranslation3D(d1x,d1y,d1z)
	mat = matrixMultiplication(mat1,mat)
	#print("dopo 2.a traslazione")
	#printMatrix(mat)
	return mat
    
def matRotoTraslaTS(s1,s2,d1,d2,ang):
	"""
		matrice di collimazione a 2 punti
		s1	vecchio centro
		s2	vecchio punto di allineamento
		d1	nuovo centro
		d2	nuovo punto di allineamento
		rototrasla in XY e trasla in Z
	"""
	s1x,s1y,s1z = s1
	s2x,s2y,s2z = s2
	d1x,d1y,d1z = d1
	d2x,d2y,d2z = d2
	# calcola matrice di trasformazione
	mat = matTranslation3D(-s1x,-s1y,-s1z)
	#print("dopo 1.a traslazione")
	#printMatrix(mat)
	# ruota la linea s1-s2 sull'orizzontale
	#a1 = math.atan2(s2y-s1y,s2x-s1x)
	# ruota per allinearla a d1-d2
	#a2 = math.atan2(d2y-d1y,d2x-d1x)
	#print("ruota di a2-a1=",0)
	#print("a2=",a2)
	#print("a1=",a1)
	mat1 = matRotation3DinZ(ang)
	mat = matrixMultiplication(mat1,mat)
	#print("dopo rotazione")
	#printMatrix(mat)
	#print("trasla nel punto d1",d1x,d1y,d1z)
	mat1 = matTranslation3D(d1x,d1y,d1z)
	mat = matrixMultiplication(mat1,mat)
	#print("dopo 2.a traslazione")
	#printMatrix(mat)
	return mat

def trasformaPunti2D(points,mat):
	"""
		trasforma i points secondo la matrice mat (NB: usa coordinata omogenea)
		conservando gli parametri della lista
	"""
	nDim = len(mat)
	for i in points:
#		print "il nodo %s (%7.3f %7.3f)" % (i[0],i[1],i[2]),
		cds = []
		for [a,b,c] in mat:
			cds.append(a*i[1]+b*i[2]+c)
		i[1],i[2] = cds[0],cds[1]
#		print 'diventa (%7.3f %7.3f)' % (i[1],i[2])
	return points

def trasformaPunti3D(points,mat):
	"""
		trasforma i points secondo la matrice mat
	"""
	for i in points:
#		print "il nodo %s (%7.3f %7.3f%7.3f)" % (i[0],i[1],i[2],i[3]),
		cds = []
		for [a,b,c,d] in mat:
			cds.append(a*i[1]+b*i[2]+c*i[3]+d)
		i[1],i[2],i[3] = cds[0],cds[1],cds[2]
#		print 'diventa (%7.3f %7.3f%7.3f)' % (i[1],i[2],i[3])
	return points

def EchelonNF(mat):
	"""
		2013-02-12
		calcola la Echelon Normal Form
		NB: gestisce matrici rettangolari
	"""
	nr = len(mat)
	nc = len(mat[0])
	# print "discesa"
	for i in range(nr):
		# print "riga",i
		# determina il pivot
		# dovrebbe essere da i in avanti perchè dietro dovrebbero essere nulli
		# però, poichè non facciamo lo swap delle righe, ripartiamo dall'inizio
		for j in range(nc):
			if abs(mat[i][j]) > 1E-10:
				# print "controllo",i,j
				cp = j
				piv = mat[i][cp]
				# print "trovato pivot riga",i,"in colonna",cp,"e vale",piv,"normalizzo la riga"
				for l in range(nc):	# forse si può partire da cp che è la prima colonna nonnulla
					mat[i][l] = mat[i][l]/piv
				# elimina dalla colonna pivot in avanti
				for j in range(i+1,nr):
					# print "elimino all'ingiu riga",j
					k = mat[j][cp]
					for l in range(cp,nc):
						mat[j][l] = mat[j][l]-k*mat[i][l]
				break
		# printMatrix(mat)
	# print "salita"
	for i in reversed(range(nr)):
		#print "riga",i
		# determina il pivot (vedi sopra a proposito della ricerca del pivot)
		for j in range(nc):
			if abs(mat[i][j]) > 1E-10:
				cp = j
				piv = mat[i][j]
				#print "trovato pivot riga",i,"in colonna",cp,"e vale",piv
				for j in range(i):
					#print "elimino riga",j
					k = mat[j][cp]/piv
					for l in range(cp,nc):
						mat[j][l] = mat[j][l]-k*mat[i][l]
				break
		#printMatrix(mat,nr,nc)
	return mat

def rank(mat):
	"""
		calcola i vettori colonna linearmente indipendenti
		e registra la riga e colonna del pivot;
		a rigore non sarebbe necessario registrare la riga, però voglio sperimentare
		la possibilità di operare senza swappare le righe pertanto le registro
	"""
	nr = len(mat)
	nc = len(mat[0])
	rng = []
	for i in range(nr):
		# print "riga",i
		for j in range(nc):
			# print "colonna",j
			if mat[i][j] > 1E-10:
				# print "      ",mat[i][j],"è un pivot"
				rng.append([i,j])
				break
	return rng

def prelievoColonne(mat,list):
	"""
		preleva le colonne indicate in list dalla matrice mat
		e le restituisce in newMat;
		presume che list sia corretta, cioè non contenga riferimenti a colonne inesistenti;

		è tagliata sul columnSpace per il quale arriva la lista rng i cui items sono duple
		riga,colonna del pivot;
	"""
	nr = len(mat)
	newMat = []
	for r in range(nr):
		tmp = []
		for k,c in list:
			tmp.append(mat[r][c])
		newMat.append(tmp)
	return newMat

def matrixInvert(mat):
	"""
		calcola l'inversa della matrice

		codici di errore:
			-1	matrice non quadrata
			-2	matrice singolare
	"""
	nr = len(mat)
	nc = len(mat[0])
	if nr == nc:
		if nr == 1:
			mat[0][0] = 1/mat[0][0]
			return mat
		else:
			# crea la matrice aggiunta
			matId = identity(nr)
			mat = adjoint(mat,nr,matId)
			# printMatrix(mat)
			# genera la echelon form
			mat = EchelonNF(mat)
			# printMatrix(mat)
			# calcola il rank (nelle prime nc colonne!)
			rng = rank(mat)
			rk = len(rng)
			# print "Rank =",rk,"colonne del range",rng
			if rk == nr:
				# preleva la parte destra della matrice ridotta
				list = []
				for i in range(nr):
					list.append([i,nr+i])
				# print "lista è",list
				# deve prelevare la inversa
				mat = prelievoColonne(mat,list)
				return mat
			else:
				# print "la matrice non è invertibile"
				return -2
	else:
		# print "matrice non quadrata: non esiste l'inversa!"
		return -1

def lsm(old,new):
	"""
		l(est) s(quare) m(ethod)

		questa routine non è ottimizzata ma è ripresa da uno studio di algebra lineare
		e si porta dietro tutto l'armamentario lì definito; potrebbe essere pesantemente
		riorganizzata
	"""
#	print "matrice di origine"
#	printMatrix(old)
#	print "matrice di destinazione"
#	printMatrix(new)
	# calcolo la trasposta
	oldT = matrixTranspose(old)
	mat = matrixMultiplication(oldT,old)
#	print "matrice normale"
#	printMatrix(mat)
#	print "Termine noto"
	tn = matrixMultiplication(oldT,new)
#	print tn
#	print "matrice inversa"
	iMat = matrixInvert(mat)
#	printMatrix(iMat)
	sol = []
	for i in range(len(tn)):	# per ogni colonna del termine noto
#		print "---------- calcolo ----------",i
		myList = [[0,i]]	# questa strana forma è dovuta a precedenti (vedi rank())
		tmp = prelievoColonne(tn,myList)
#		print "termine noto"
#		printMatrix(tmp)
		solTmp = matrixMultiplication(iMat,tmp)
#		print "soluzione",solTmp
		tmp = matrixTranspose(solTmp)
#		print tmp
		sol.append(tmp[0])	# questa operazione è empirica, serve ad eliminare un livello di parentesi, controllare
#	print "soluzione"
#	printMatrix(sol)
	return sol

def erroreLsm(res,att):
	"""
		calcola l'errore medio dei RESultati rispetto agli ATTesi
		il nome riferisce a LSM ma in realtà può essere usata per qualsiasi situazione
	"""
	nP = len(res)	# numero di punti
	nD = len(res[0])		# dimensione
	errMin,errMax,errTot = 9.99e99,-9.99e99,0.0	# errore minimo, massimo, totale
	for i,r in enumerate(res):
		err = 0
		for k,rv in enumerate(r):
			d = rv-att[i][k]
#			print i,k,rv,att[i][k]
			err += d**2
		err = math.sqrt(err)
#		print "errore",i,err
		errTot += err
		if err < errMin:
			errMin = err
		if err > errMax:
			errMax = err
	return errMin,errMax,errTot/nP

# ----- projection functions -----------
"""
	contiene alcune routines per la riproiezione di cds
	geografiche <-> geocentriche <-> topocentriche;
	al momento è usata solo la trasformazione geocentriche->topocentriche
	per il trattamento delle misure gps
"""

def wgs84():
	"""
		fornisce i dati dell'ellissoide wgs84
	"""
	a = 6378137.00	# semiasse maggiore
	f = 298.2572236
	r = 1-1/f
	b = a*r	# semiasse minore
	return a,b

def geocentriche2wgs84(x,y,z):
	"""
		trasforma coordinate geocentriche in coordinate WGS84
		(controllata con EPSG Guide 7-2 pag.94)
	"""
	# ellissoide WGS84
	a,b = wgs84()
#	print "semiassi:",a,b
	# dati preliminari
	e2 = (a**2-b**2)/a**2
#	print "e2:",e2
	eps = e2 / (1-e2)
#	print "eps:",eps
	p = math.sqrt(x**2 + y**2)
#	print"p:",p
	q = math.atan2(z*a,p*b)
#	print "q:",q
	# cordinate geografiche
	lat = math.atan2(z+eps*b*math.sin(q)**3,p-e2*a*math.cos(q)**3)
#	print "lat:",lat
	lon = math.atan2(y,x)	# anche formula Vanicek/Krakiwsky o Vermeille in Ligas/Banasik
	v = a/math.sqrt(1-e2*math.sin(lat)**2)
#	print "v:",v
	h = p/math.cos(lat)-v 
	return lon,lat,h

def geocentriche2topocentriche(X,Y,Z,X0,Y0,Z0):
	"""
		geocentric to topocentric
		(controllata con EPSG Guide 7-2 pag.97)
	"""
	lon0,lat0,h0 = geocentriche2wgs84(X0,Y0,Z0)
#	print lon0,lat0
	U = -(X-X0)*math.sin(lon0) + (Y-Y0)*math.cos(lon0)
	V = -(X-X0)*math.sin(lat0)*math.cos(lon0) - (Y-Y0)*math.sin(lat0)*math.sin(lon0) + (Z-Z0)*math.cos(lat0) 
	W = (X-X0)*math.cos(lat0)*math.cos(lon0) + (Y-Y0)*math.cos(lat0)*math.sin(lon0) + (Z-Z0)*math.sin(lat0)
	return U,V,W
    
def topocentriche2geocentriche(U,V,W,X0,Y0,Z0):
	"""
		topocentric to geocentric
		(controllata con EPSG Guide 7-2 pag.98)
	"""
	lon0,lat0,h0 = geocentriche2wgs84(X0,Y0,Z0)
#	print lon0,lat0
	X = (X0) - (U)*math.sin(lon0) - (V)*math.sin(lat0)*math.cos(lon0) + (W)*math.cos(lat0)*math.cos(lon0)
	Y = (Y0) + (U)*math.cos(lon0) - (V)*math.sin(lat0)*math.sin(lon0) + (W)*math.cos(lat0)*math.sin(lon0) 
	Z = (Z0) + (V)*math.cos(lat0) + (W)*math.sin(lat0)
	return X,Y,Z    

# ----- I/O functions -------------------

def loadFile(fname):
	"""
	Importa il file fname
	"""
	lib = []
	f = open(fname, 'r')
	try:
		f.readline()
	except:
		f = codecs.open(fname, 'r', 'cp1252')
	for data in f:
		if '***** Relazione  Tecnica *****' in data:
			break	# finisce la lettura
		# pulisce la riga
		data = data.rstrip('\n')
		data = data.rstrip('\r')
		lib.append(data)    
	return lib

def openLibretto_vertici(libretto):
	"""
		legge un libretto di campagna con le codifiche pregeo
		al momento legge solo linee:
		- 1 stazione celerimetrica
		- 2 punto osservato celerimetrico;
		le letture vengono lasciate immutate (saranno interpretate successivamente)
		viene solo aggiunto il campo della stazione di osservazione di ogni vertice;
		ora legge da lista (2013-10-22)
		viene registrata anche la riga del libretto (2013-10-22)
	"""
	celerimensura = []
	stazione = []
	nRiga = -1	# così il primo valore è zero
	isFirst =True
	codStaz = 0	# si potrebbe usare anche questo come flag della prima stazione
	for line in libretto:
		nRiga += 1	
		tmp = line.split('|')
		cod = tmp.pop(0)		# questo elimina il primo campo	
		if cod == '1'or cod == '2':	# così scarta tutto quello che non interessa
			tmp.pop()				# elimina l'ultimo campo sempre vuoto
#			print "Riga:",nRiga,tmp[1]
			# è una lettura gps?
			if not ',' in tmp[1]:	# altrimenti sarebbe una lettura gps
				# E' una stazione?
				if cod == '1':
					codStaz = tmp[0]
#					print "stazione",codStaz
					# E' successiva alla prima?
					if not isFirst:
						celerimensura.append(stazione)
					isFirst = False
					# pulisce la lista stazione
					stazione = []
				# aggiunge il campo STAZIONE
				tmp.append(codStaz)
				# aggiunge la riga del libretto
				tmp.append(nRiga)
				# salva la lettura
				stazione.append(tmp)
	# salva l'ultima stazione se c'è
	if stazione:
		celerimensura.append(stazione)
	return celerimensura

def openLibretto_gps(libretto):
	"""
		legge un libretto di campagna con le codifiche pregeo
		- 1 baseline GPS
		- 2 lettura GPS;
		le letture vengono lasciate immutate (saranno interpretate successivamente)
		viene solo aggiunto il campo della stazione di osservazione di ogni vertice;
		ora legge da lista (2013-10-22)
		viene registrata anche la riga del libretto (2013-10-22)
		forse non ha senso aspettarsi più baseline gps (2013-10-22)
	"""
	gps = []
	lista = []
	nRiga = -1	# così il primo valore è zero
	note = 'gps'
	isFirst =True
	for line in libretto:
		nRiga += 1	
		tmp = line.split('|')
		cod = tmp.pop(0)		# questo elimina il primo campo	
		if cod == '1'or cod == '2':	# così scarta tutto quello che non interessa
			tmp.pop()				# elimina l'ultimo campo sempre vuoto
			# è una lettura gps?
			if ',' in tmp[1]:
				pnt = tmp[0]
				tmp2 = tmp[1].split(',')
				dx,dy,dz = float(tmp2[0]),float(tmp2[1]),float(tmp2[2])
				if cod == '1':
#					print "baseline",basLne
					if tmp[2]:
						hp = float(tmp[2])
					else:
						hp = 0.0
					# Ã¨  successiva alla prima?
					if not isFirst:
						gps.append(lista)
					isFirst = False
					# pulisce la lista
					lista =[]
				else:
					#print('non è baseline')
					if tmp[4]:
						hp = float(tmp[4])
					else:
						hp = 0.0
#				print "Riga:%s Gps:%s (%f %f %f) h=%f note:%s" % (nRiga,pnt,dx,dy,dz,hp,note)
				# salva la lettura
				lista.append([pnt,dx,dy,dz,hp,note,nRiga])
	# salva l'ultima lista se c'è
	if lista:
		gps.append(lista)
#	print 'fine lettura',fname
	return gps

def openLibretto_allSquad(libretto):
	"""
		legge un libretto di campagna con le codifiche pregeo
		- 4 punti di allineamento
		- 5 punto generato, distanza da p0 e squadro
		ora legge da lista (2013-10-22)
		viene registrata anche la riga del libretto (2013-10-22)
	"""
	RilAllin = []
	nRiga = -1	# così il primo valore è zero
	isFirst =True
	p0,p1 = -1,-1
	codStaz = 0			# gli RilAllin hanno stazione nulla
	for line in libretto:
		nRiga += 1	
		tmp = line.split('|')
		cod = tmp.pop(0)
		# allineamento?
		if cod == '4':
			p0,p1 = tmp[0],tmp[1]
		# squadro?
		if cod == '5':
			p,a,s,note = tmp[0],tmp[1],tmp[2],tmp[3]
			RilAllin.append([p,p0,p1,a,s,note,codStaz,nRiga])
	return RilAllin

def openLibretto_contorni(file):
	"""
		restituisce le liste:
		- ctrn[]	contorni particelle (lista [v0,v1,,....])
		- style[]	stile dei contorni (lista [NC/NT/RC/RT])

		mancano da leggere:
		- identificativi particelle
	"""
	ctrn = []
	style = []
	nv = 0
	tmp2 = []
	for data in file:
		# (righe con codice iniziale 7)
		if data[0:1] == '7':
			# decodifica i dati
			tmp = data.split('|')
			# elimina il campo codice
			tmp.pop(0)
			# registra, eliminando, il numero di vertici
			num = int(tmp.pop(0))
			# setta il contatore dei vertici se è una riga nuova
			if num:
				nv = num  
			# definisce il numero max di vertici della riga
			k = min(10,nv)
			# li parcheggia in tmp2
			for v in tmp[0:k]:   
				if (v != 'RC' and v != 'RT' and v != 'NC' and v != 'NT' and v != 'RP' and v != 'NP'): 
					tmp2.append(v)
				else:
					k = k -1                        		
			#print("nv=",nv,"num=",num,"tmp2=",tmp2)	
			tmpSty = tmp[k]
			# aggiorna il contatore dei vertici
			nv -= k
			# se è finito il contorno trasferisce nelle liste
			if nv == 0:
				ctrn.append(tmp2)
				style.append(tmpSty)
				#print(tmpSty)
				tmp2 = []		            
	# elimina il segno percentuale ove necessario
	for k in ctrn:
		for i,v in enumerate(k):
			if v[-1] in ['%','L','T']:
				k[i] = v	#NUOVO
	return ctrn,style

def openEDM_pf(file):
	"""
		restituisce la lista pf[] dei punti fiduciali di un EdM (lista di liste: [pfId,x,y,z=0])
	"""
	pf = []
	for data in file:
		# (righe con codice iniziale 8)
		if data[0:1] == '8':
			# (e secondo campo "PFxx/xxxx/xxxx")
			if data[2:4] == 'PF':
				tmp = data.split("|")
				note = str(tmp[1])
				pf.append([note,float(tmp[3]),float(tmp[2]),0])	# NB: mette prima Nord e poi Est
	return pf

def openEDM_vertici(file):
	"""
		restituisce la lista lut[]	vertici della particella (lista di liste: [pfId,x,y,z=0])
	"""
	lut = []
	for data in file:
		# (con codice iniziale 8)
		if data[0:1] == '8':
			tmp = data.split("|")
			cod = tmp[1][-1]
			pnt = tmp[1]
			# (e secondo campo con "%" o 'T' o 'L'))
			if cod in ('%','T','L'):
				note = 'EdM-'+str(pnt)
				# controllare numerazione di pregeo
				lut.append([pnt,float(tmp[3]),float(tmp[2]),0])	# NB: mette prima Nord e poi Est
	#print(lut)
	return lut

# ------- inquiry/report functions ---------

def pointArchivioCds(archivio,p):
	"""
		fornisce il progressivo e le coordinate del punto p
		forza il parametro "p" a stringa
		in caso di fallimento restituisce un id = -1
	"""
	for k,i in enumerate(archivio):
#		print "controllo se",str(p),"è uguale a",i[0]
		if i[0] == str(p):
#			print "trovato",p
			return k,i[1],i[2],i[3]
#	print "non trovato",p
	return [-1,0,0,0]

def printList(lib):
	"""
		generica routine di stampa delle liste
		stampa anche la numerazione degli item
	"""
	for i,l in enumerate(lib):
		print(i,l)

def stazioniLista(libretto):
	"""
		elenco delle stazioni
	"""
	list = []
	for line in libretto:
		tmp = line.split('|')
		if tmp[0] == '1':
			if ',' not in tmp[2]:	# controlla che non sia una lettura gps
				list.append(tmp[1])                
			else:
				list.append('')
				#print("trattasi di libretto gps")
	#print("trattasi di libretto celerimetrico")
	return list

def lettureFraStazioni(libretto):
	"""
		fornisce tutte le letture fra stazioni
		potrei usare la modularità rispetto a 100, però preferisco usare il codice di inizio riga
	"""
	print('cerco le stazioni')
	stazList = stazioniLista(libretto)
	# cerca le lettura fra stazioni
	lst = []
	for s0 in stazList:
#		print "cerco la stazione",s0
		for k,l0 in enumerate(libretto):
			tmp0 = l0.split("|")
#			print tmp0
			if tmp0[0] == '1' and tmp0[1] == s0:
				break
#		print "cerco la stazione mirata",s1
		for j in range(k+1,len(libretto)):	# dalla stazione precedente alla fine del libretto
			l1 = libretto[j]
			tmp1 = l1.split("|")
#			print tmp1
			if tmp1[0] == '1':
				break	# ha trovato un'altra stazione
			if tmp1[0] == '2' and tmp1[1] in stazList:
#				print s0,tmp1[1]
				lst.append([s0,tmp1[1]])
	return lst

def osservazioniCelerimetriche(libretto):
	#print('cerco le stazioni')
	stazList = stazioniLista(libretto)
	# cerca le lettura fra stazioni
	lst = []
	for s0 in stazList:
		for k,l0 in enumerate(libretto):
			tmp0 = l0.split("|")
			if tmp0[0] == '1' and tmp0[1] == s0:
				for j in range(k+1,len(libretto)):
					l1 = libretto[j]
					tmp1 = l1.split("|")
					if tmp1[0] == '2':
						print('Stazione',s0,'->',l1)

def distanzeRidotte(archivio,a_giro):
	"""
		fornisce per ogni stazione le distanza ridotte di ogni vertice misurato
		parte di codice uguale a polar2rect()
	"""
	list = []
	for line in archivio:
		tmp = line.split('|')
		cod = tmp.pop(0)		# questo elimina il primo campo	
		if cod == '1'or cod == '2':	# così scarta tutto quello che non interessa
			tmp.pop()				# elimina l'ultimo campo sempre vuoto
#			print "Riga:",nRiga,tmp[1]
			if ',' not in tmp[1]:	# altrimenti sarebbe una lettura gps
				# E' una stazione?
				if cod == '1':
					staz = tmp[0]
#					print "stazione",codStaz
				else:
					if len(tmp) == 6:	# punto dotato di altimetria
						pnt = tmp[0]
						# converte zenith in radianti
						av = toRad(float(tmp[2]),a_giro)
						av -= 3*math.pi/2		# NB: ricordare angolo nullo allo zenith
						# calcolo la distanza ridotta
						d = float(tmp[3])
						dr = d*abs(math.cos(av))
						list.append("distanza %5s - %20s: %12.3f" % (staz,pnt,dr))
	return list

def pfLista(libretto):
	"""
		elenco dei PF rilevati
	"""
	list = []
	for line in libretto:
		tmp = line.split('|')
		if tmp[0] == '2':
			if '/' in tmp[1]:
				list.append(tmp[1])
	return list

# ---------- operation functions ---------

def elaborazioneGps(myGps):
	"""
		elaborazione dei gps (procedura sperimentale)
		trasformazione le cds geocentriche delle misure gps in cds topocentriche;
		una formalizzazione corretta presupporrebbe la conoscenza del punto di
		emanazione del sistema di riferimento locale; poichè il dato non è noto a
		priori surroghiamo con una approssimazione: adottiamo il centro di emanazione
		nella baseline del rilievo, lasciando alla collimazione con i fiduciali la
		referenziazione nel corretto CRS;
		rimane il problema dell'altitudine da controllare bene;
	"""
	tmp = []
	note = 'gps'
	for g in myGps:
		# preleva la baseline (cds geocentriche)
		reg = g[0]
		basLn = reg[0]
		x0,y0,z0 = reg[1:4]
		nRiga = reg[-1]
		# trasforma in cds topocentriche
		u,v,w = geocentriche2topocentriche(x0,y0,z0,x0,y0,z0)
		tmp.append([basLn,u,v,w,note,basLn,nRiga])
		# calcola tutti i vertici
		for i in range(1,len(g)):
			# preleva i punti (cds geocentriche)
			reg = g[i]
			pId = reg[0]
			dx,dy,dz = reg[1:4]
			# trasforma in cds topocentriche
			u,v,w = geocentriche2topocentriche(x0+dx,y0+dy,z0+dz,x0,y0,z0)
			nRiga = reg[-1]
			tmp.append([pId,u,v,w,note,basLn,nRiga])	
	return tmp

def polar2rect(rilievo,a_giro):
	"""
		riceve il rilievo di una stazione del libretto in coordinate polari
		e la converte in coordinate rettangolari
		NB: il senso di rotazione del teodolite è inverso, quindi l'angolo va cambiato di segno
	"""
	coordRet = []
	# inserisce la stazione
	stazione = rilievo[0]
	stazDiMira = stazione[0]
	hs = float(stazione[1])
	nRiga = stazione[-1]
	coordRet.append([stazDiMira,0.0,0.0,0.0,stazione[2],stazDiMira,nRiga])
	for i in range(1,len(rilievo)):
		punto = rilievo[i]
		cod = punto[0]
		x,y,z = 0.0,0.0,0.0
		nRiga = punto[-1]
		if len(punto) == 8 or len(punto) == 7:	# punto dotato di altimetria
			# converte zenith in radianti
			av = toRad(float(punto[2]),a_giro)
			av -= 3*math.pi/2		# NB: ricordare angolo nullo allo zenith
			# calcolo la distanza ridotta
			d = float(punto[3])
			dr = d*abs(math.cos(av))
			# altitudine
			hp = float(punto[4])
			z = (hs + d*math.sin(av) - hp)
			if len(punto) == 8:			
				nota = punto[5]
			elif len(punto) == 7:
				nota = ""                 
		if len(punto) == 6:	# punto senza altimetria (NB: c'è sempre un campo vuoto alla fine)
			dr = float(punto[2])
			nota = punto[3]
		#else:
		#	dr = 0.0	# dovrebbe mettere tutto a zero
		#	nota = 'polar2rect: numero di parametri non compatibili'
		#	stazDiMira = 0
		# converte azimuth in radianti
		ah = toRad(-float(punto[1]),a_giro)
		# coordinate rettangolari nel piano orizzontale
		x = dr*math.cos(ah)
		y = dr*math.sin(ah)
		# salva
		coordRet.append([cod,x,y,z,nota,stazDiMira,nRiga])
#		print 'punto %s %f %f %f' % (str(punto),x,y,z)
	return coordRet

def rect2polar(x,y,z,x0,y0,z0,a_giro):
	"""
		trasforma coordinate rettangolari in polari relative al punto(x0,y0,z0)
		nella notazione angolare selezionata (a_giro)
	"""
	dx,dy,dz = x-x0,y-y0,z-z0
	ah = fromRad(math.atan2(dy,dx),a_giro)
	dr = math.sqrt(dx**2+dy**2)
	# NB1: ricordare angolo nullo allo zenith
	# NB2: questo va bene per il mio teodolite che ha il goniometro verticale ruotato
	av = fromRad(math.atan2(dz,dr)+3*math.pi/2,a_giro)
	d = math.sqrt(dr**2+dz**2)
	return ah,av,d

def collimazioneStazione(rilievo,archivio,tipologia):
	"""
		Esegue la trasformazione di punti in modo che i punti sorgenti [srcList]
		vadano a collimare con i destinatari [destList].

		questa si occupa del trattamento di ogni stazione (rilievo) per adattarla al
		riferimento  della prima stazione (archivio);

		la collimazione si limita al solo piano XY; la Z subisce una pura traslazione in media;

		la collimazione viene eseguita con il metodo dei 2 punti, la nuova stazione prende
		la posizione mirata dalla precedente stazione ed il rilievo viene ruotato per portare
		la posizione mirata della precedente stazione sulla sua posizione già in archivio;
		questa modalità previene l'uso di scalature e consente quindi di apprezzare la qualità 
		del rilievo nella sua interezza;

		questo metodo impone che nella poligonale due stazioni contigue si devono sempre vedere;
		questo metodo però non funziona quando c'è una baseline gps, pertanto quando non trova
		la stazione mirante fra i punti della stazione, prende il primo punto ribattuto;
	"""
	az_xy = 0
	az_en = 0
	xa = 0
	ya = 0
	xb = 0
	yb = 0
	ea = 0
	na = 0
	eb = 0
	nb = 0
	ang = 0
	cod = ""
	#print("collimazione: ricevo archivio",archivio,"rilievo",rilievo)
	collimati = []
	# controlla se la stazione è già stata osservata
	newStaz = rilievo[0][0]
	s1 = [0,0,0]
#	print "nuova stazione",newStaz,"dalla posizione",0,0,0
	pos,x,y,z = pointArchivioCds(archivio,newStaz)
	if pos >= 0:
		d1 = [x,y,z]
#		print "alla posizione",x,y,z
		# stazione mirante
		prevStaz = archivio[pos][5]
		pos,x,y,z = pointArchivioCds(rilievo,prevStaz)
		if pos >= 0:
			s2 = [x,y,z]
#			print "la stazione mirante",prevStaz,"dalla posizione",x,y,z
			pos,x,y,z = pointArchivioCds(archivio,prevStaz)
			if pos >= 0:
				d2 = [x,y,z]
#				print "alla posizione",x,y,z
				# matrice di trasformazione
				if tipologia == 2:
					for k in archivio:
						if newStaz == k[0]:
							#print("trovato il punto stazione in archivio",newStaz)
							#print('trovata',k)
							xa = k[1]
							ya = k[2]
							#break
					for k in rilievo:
						if newStaz == k[0]:
							#print("trovato il punto stazione in rilievo",newStaz)
							#print('trovata',k)
							ea = k[1]
							na = k[2]
							#break
					for k in archivio:		
						if cod == k[0]:
							#print("trovato il punto orientamento in archivio",cod)
							#print('trovata',k)
							xb = k[1]
							yb = k[2]
							#break
					for k in rilievo:		
						if cod == k[0]:
							#print("trovato il punto orientamento in rilievo",cod)
							#print('trovata',k)
							eb = k[1]
							nb = k[2]
							#break
					az_xy = math.atan2((xb-xa),(yb-ya))
					az_en = math.atan2((eb-ea),(nb-na))
					#print('az_xy',az_xy)
					#print('az_en',az_en)
					ang = az_en - az_xy
					#print('az_en - az_xy',ang)
					#print(ang)
				mat = matRotoTraslaTS(s1,s2,d1,d2,ang)
#				print "matrice di trasformazione",mat
				collimati = trasformaPunti3D(rilievo,mat)
			else:
				print('la stazione mirante non è in archivio (questo errore non si dovrebbe MAI verificare)')
		else:
			#print('la stazione mirante',prevStaz,'non è stata osservata dalla stazione',newStaz) 
			#print('cerchiamo il primo punto ribattuto')
			for i in range(1,len(rilievo)):
				cod = rilievo[i][0]
				for k in archivio:
					if cod == k[0]:
						#print("trovato il punto ribattuto",cod)
						pos,x,y,z = pointArchivioCds(rilievo,i)
						s2 = [x,y,z]
						pos,x,y,z = pointArchivioCds(archivio,k)
						d2 = [x,y,z]
						if tipologia == 2:
							for k in archivio:
								if newStaz == k[0]:
									#print("trovato il punto stazione in archivio",newStaz)
									#print('trovata',k)
									xa = k[1]
									ya = k[2]
									#break
							for k in rilievo:
								if newStaz == k[0]:
									#print("trovato il punto stazione in rilievo",newStaz)
									#print('trovata',k)
									ea = k[1]
									na = k[2]
									#break
							for k in archivio:		
								if cod == k[0]:
									#print("trovato il punto orientamento in archivio",cod)
									#print('trovata',k)
									xb = k[1]
									yb = k[2]
									#break
							for k in rilievo:		
								if cod == k[0]:
									#print("trovato il punto orientamento in rilievo",cod)
									#print('trovata',k)
									eb = k[1]
									nb = k[2]
									#break
							az_xy = math.atan2((xb-xa),(yb-ya))
							az_en = math.atan2((eb-ea),(nb-na))
							#print('az_xy',az_xy)
							#print('az_en',az_en)
							ang = az_en - az_xy
							#print('az_en - az_xy',ang)
							#print(ang)
						# matrice di trasformazione
						mat = matRotoTraslaTS(s1,s2,d1,d2,ang)
#						print "matrice di trasformazione",mat
						collimati = trasformaPunti3D(rilievo,mat)	# si può usare matrixMultiplication()
						break
	else:
		print('la stazione',newStaz,'non è mirata da alcuna stazione precedente in archivio')
	return collimati

def allEsquad(archivio,allin):
	"""
		gestione di un allineamento e squadro
	"""
	cod,n1,n2,a,s,note,staz,nRiga = allin
	a = float(a)
	s = -float(s)	# la routine usa una convenzione diversa per il segno dello squadro
	# legge le cds dei punti
	pos,x1,y1,z1 = pointArchivioCds(archivio,n1)
	pos,x2,y2,z2 = pointArchivioCds(archivio,n2)
	# costruisce i vettori allineati
	ux,uy = x2-x1,y2-y1
	mu = math.sqrt(ux**2+uy**2)
	if mu == 0.0:
		#verificare print('ERRORE: allineamento',%s,'su',%s(%f-%f)-%s(%f-%f),'da un vettore nullo',%(cod,n1,x1,y1,n2,x2,y2)')
		return cod,-1,-1,-1,'errore',0,0
	else:
		# vettore ortogonale (nel piano orizzontale)
		vx,vy = -uy,ux
		# calcolo il modulo di u
		# punto allineato (piede H)
		hx,hy = x1+(a/mu)*ux,y1+(a/mu)*uy
		# punto P
		x,y = hx+(s/mu)*vx,hy+(s/mu)*vy
		# NB: altimetria non gestita
		z = 0.0
		return [cod,x,y,z,note,staz,nRiga]

def collimazione2PF(oldCds,newCds,archivio):
	"""
		esegue la georeferenzazione a due PF, il primo di centratura ed
		il secondo di allineamento;
		- oldC	coordinate origine
		- newC	coodinate destinazione
	"""
	# coordinate del centro 
	oldC = oldCds[0]
	newC = newCds[0]
#	print "centra da",oldC,"a",newC
	# coordinate dell'allineamento
	oldA = oldCds[1]
	newA = newCds[1] 
#	print "allinea da",oldA,"a",newA
	# esegue la rototraslazione (collimazione a 2 punti)
	mat = matRotoTrasla(oldC,oldA,newC,newA)
#	print "matrice di rototraslazione"
#	printMatrix(mat)
	archivio = trasformaPunti3D(archivio,mat)	# non si può usare matrixMultiplication() perchè ci sono altri elementi da conservare
	return archivio

def collimazione3PF(oldCdsList,newCdsList,rilievo):
	"""
		Esegue la trasformazione di punti nello spazio dei PF
		(lavora nel solo piano XY, valutare se adottare lo spazio XYZ)

		NB: occorre per forza LSM perchè sono 3 PF per 2 dim
	"""
	# mette a 0.0 la z
	for v in oldCdsList:
		v[2] = 0.0
	for v in newCdsList:
		v[2] = 0.0
	mat = lsm(oldCdsList,newCdsList)
#	printMatrix(mat)
	collimati = trasformaPunti2D(rilievo,mat)	# non si può usare matrixMultiplication() perchè ci sono altri elementi da 
	return collimati

# ======================== classe celerimensura ========================

class celerimensuraDlg(QDialog):
	"""
		Dialogo per misure celerimetriche

		traccia il layout della misura celerimetrica
		l'origine è nel punto a terra sotto la stazione
	"""

	#	style for texts
	pen = QPen()
	pen.setWidthF(2.0)		# controlla la dimensione del punto
	pen.setColor(QColor(0,0,0))
	font = QFont()
	penText = [pen,font]
	#	style for points
	pen = QPen()
	pen.setWidthF(5.0)		# controlla la dimensione del punto
	pen.setColor(QColor(0,0,255))
	font = QFont()
	font.setPointSize(10)
	penPoint = [pen,font]
	#	per for edges
	pen = QPen()
	pen.setWidthF(1.5)
	pen.setColor(QColor(0,0,0))
	font = QFont()
	font.setPointSize(10)
	penEdge = [pen,font]

	maxVprt = 500	# massima dimensione della viewport

	def __init__(self,a_giro,s,hs,pnt,av,d,h):
		"""
			Inizializza la maschera
			a_giro	notazione angolare (360/400)
			s	stazione id
			hs	altezza stazione
			pnt	point id
			av	vertical angle
			d	distance
			h	prisma height
		"""
		QDialog.__init__(self)
		# impostazione interfaccia utente
		self.setWindowTitle('Misura celerimetrica')
		self.resize(self.maxVprt,self.maxVprt+120)

		# definisce la scala
		self.scalH = 0.65*self.maxVprt/d
		self.scalV = 5*self.scalH

		self.a_giro = a_giro
		self.s = s
		self.hs = hs
		self.h = h
		self.pnt = pnt
		self.av = av
		self.d = d
		self.dr = 0.0
		self.dz = 0.0
		self.disl = 0.0

		# --------layout verticale ------------
		vBox = QVBoxLayout()
		self.setLayout(vBox)

		# ------------ canvas ------------------
		self.canvas = QGraphicsScene(self)
		view = QGraphicsView(self.canvas)
		view.show()
		vBox.addWidget(view)

		# ----- parametri ---------
		self.compute()

		lbl = QLabel(self)
		lbl.setText('stazione')
		lbl.setGeometry(QRect(15,480,100,20))
		lbl = QLabel(self)
		lbl.setText(self.s)
		lbl.setAlignment(Qt.AlignRight)
		lbl.setGeometry(QRect(120,480,120,20))

		lbl = QLabel(self)
		lbl.setText('punto mirato')
		lbl.setGeometry(QRect(255,480,100,20))
		lbl = QLabel(self)
		lbl.setText(self.pnt)
		lbl.setAlignment(Qt.AlignRight)
		lbl.setGeometry(QRect(360,480,120,20))

		lbl = QLabel(self)
		lbl.setText('altezza stazione')
		lbl.setGeometry(QRect(15,500,100,20))
		lbl = QLabel(self)
		lbl.setText("%12.3f"%(self.hs))
		lbl.setAlignment(Qt.AlignRight)
		lbl.setGeometry(QRect(120,500,120,20))

		lbl = QLabel(self)
		lbl.setText('altezza di mira')
		lbl.setGeometry(QRect(255,500,100,20))
		lbl = QLabel(self)
		lbl.setText("%12.3f"%(self.h))
		lbl.setAlignment(Qt.AlignRight)
		lbl.setGeometry(QRect(360,500,120,20))

		lbl = QLabel(self)
		lbl.setText('zenith')
		lbl.setGeometry(QRect(15,520,100,20))
		self.eAv = QLabel(self)
		self.eAv.setAlignment(Qt.AlignRight)
		self.eAv.setGeometry(QRect(120,520,120,20))

		lbl = QLabel(self)
		lbl.setText('distanza')
		lbl.setGeometry(QRect(255,520,100,20))
		self.eDi = QLabel(self)
		self.eDi.setAlignment(Qt.AlignRight)
		self.eDi.setGeometry(QRect(360,520,120,20))

		lbl = QLabel(self)
		lbl.setText('distanza ridotta')
		lbl.setGeometry(QRect(15,540,100,20))
		self.eDr = QLineEdit(self)
		self.eDr.setText("%12.3f"%(self.dr))
		self.eDr.setAlignment(Qt.AlignRight)
		self.eDr.setGeometry(QRect(120,540,120,20))

		lbl = QLabel(self)
		lbl.setText('dislivello')
		lbl.setGeometry(QRect(255,540,100,20))
		self.eDsl = QLineEdit(self)
		self.eDsl.setText('%12.3f' %(self.disl))
		self.eDsl.setAlignment(Qt.AlignRight)
		self.eDsl.setGeometry(QRect(360,540,120,20))

		lbl = QLabel(self)
		lbl.setText('NB: scala in Y 5x')
		vBox.addWidget(lbl)

		# ----- bottoni -------------------
		hBox = QHBoxLayout()
		vBox.addLayout(hBox)

		self.btn = QPushButton("Refresh")
		self.btn.clicked.connect(self.refresh)
		hBox.addWidget(self.btn)

		self.btn = QPushButton("Close")
		self.btn.clicked.connect(self.close)
		hBox.addWidget(self.btn)

		self.redraw()
		self.parametri()

	def compute(self):
		av1 = (self.av/self.a_giro)*2*math.pi
		av1 -= 3*math.pi/2		# NB: ricordare angolo nullo allo zenith
		self.dr = self.d*abs(math.cos(av1))
		self.dz = self.d*math.sin(av1)
		self.disl = self.hs + self.dz - self.h

	def parametri(self):
		self.eDi.setText("%12.3f"%(self.d))
		self.eAv.setText("%12.5f"%(self.av))

	def drawText(self,text,x,y):
		"""
			Disegna un testo sul canvas
		"""
		lab = QGraphicsTextItem(text)
		lab.setPos(QPointF(x,-y))	# qui fa la trasformazione
		col = self.penText[0].color()
		lab.setFont(self.penText[1])
		lab.setDefaultTextColor(col)
		self.canvas.addItem(lab)

	def drawPoint(self,n,x,y):
		"""
			Disegna il punto N del rilievo.
		"""
		dim = self.penPoint[0].width()
		pnt = QGraphicsRectItem(x-dim/2,-y-dim/2,dim,dim)	# qui fa la trasformazione
		pnt.setPen(self.penPoint[0])
		self.canvas.addItem(pnt)
		self.drawText(str(n),x,y)

	def drawEdge(self,lab,x1,y1,x2,y2):
		"""
			Disegna il lato (x1,y1) (x2,y2).
		"""
		edge = QGraphicsLineItem(x1,-y1,x2,-y2)	# qui fa la trasformazione
		edge.setPen(self.penEdge[0])
		self.canvas.addItem(edge)
		self.drawText(lab,(x1+x2)/2,(y1+y2)/2)

	def redraw(self):
		"""
			disegna lo schema grafico
		"""
		self.penEdge[0].setColor(QColor(0,0,0))	# la pen è l'elemento 0 della lista
		self.penEdge[0].setWidthF(1.5)
		self.penText[0].setColor(QColor(0,0,0))	# la pen è l'elemento 0 della lista
		# pulisce
		self.canvas.clear()
		# disegna l'origine
		self.drawPoint('O',0,0)
		# disegna la stazione
		self.drawPoint(self.s,0,self.scalV*self.hs)
		# disegna prisma
		self.drawPoint(self.pnt,self.scalH*self.dr,self.scalV*(self.hs+self.dz))
		# disegna punto a terra
		self.drawPoint('PT',self.scalH*self.dr,self.scalV*self.disl)
		# disegna rilievo
		self.penText[0].setColor(QColor(0,0,255))
		self.drawEdge("%12.3f"%(self.hs),0,0,0,self.scalV*self.hs)
		self.drawEdge("%12.3f"%(self.d),0,self.scalV*self.hs,self.scalH*self.dr,self.scalV*(self.hs+self.dz))
		self.drawEdge("%12.3f"%(self.h),self.scalH*self.dr,self.scalV*(self.hs+self.dz),self.scalH*self.dr,self.scalV*self.disl)
		# disegna distanza ridotta
		self.penEdge[0].setColor(QColor(255,0,0))	# la pen è l'elemento 0 della lista
		self.penEdge[0].setWidthF(0.5)
		self.penText[0].setColor(QColor(255,0,0))
		self.drawEdge("%12.3f"%(self.dr),0,0,self.scalH*self.dr,0)
		# disegna dislivello
		self.drawEdge("%12.3f"%(self.disl),self.scalH*self.dr,0,self.scalH*self.dr,self.scalV*self.disl)

	def refresh(self):
		self.disl = float(self.eDsl.text())
		self.dr = float(self.eDr.text())
		self.dz = self.h + self.disl - self.hs
		self.d = math.sqrt(self.dr**2+self.dz**2)
		av1 = math.atan2(self.dz,self.dr)
		print('av1'),av1
		self.av = av1*self.a_giro/(2*math.pi)
		print('av'),self.av
		self.av += (3./4.)*self.a_giro
		print('av'),self.av

		print('dislivello %12.3f' % (self.disl))
		print('distanza ridotta %12.3f' % (self.dr))
		print('dz %12.3f' % (self.dz))
		print('distanza %12.3f' % (self.d))
		print('zenith %12.5f' % (self.av))

		self.redraw()
		self.parametri()

# ======================== classe navigatore ========================

class navigatorDlg(QDialog):
	""" Dialogo per le combinazione """

	def __init__(self,a_giro,libretto,lista):
		"""
			Inizializza la maschera per la navigazione delle stazioni.
		"""
		QDialog.__init__(self)
		# impostazione interfaccia utente
		self.setWindowTitle('Navigazione poligonale')
		self.resize(200,100)

		# inizializza variabili
		self.a_giro = a_giro
		self.libretto = libretto
		self.lista = lista
		self.nStaz=len(lista)
		self.cntr = 0

#		----------- box principale ---------------
		vBox = QVBoxLayout()
		self.setLayout(vBox)

		# ------------ navigatori ---------------
		hBox = QHBoxLayout()
		vBox.addLayout(hBox)
		self.btn = QPushButton("<")
		self.btn.clicked.connect(self.prev)
		hBox.addWidget(self.btn)
		self.eLbl = QLabel(self)
		self.eLbl.setText('')
		hBox.addWidget(self.eLbl)
		self.btn = QPushButton(">")
		self.btn.clicked.connect(self.next)
		hBox.addWidget(self.btn)

		# ------------ bottoni -------
		hBox = QHBoxLayout()
		vBox.addLayout(hBox)
		self.btn = QPushButton("Go")
		self.btn.clicked.connect(self.update)
		hBox.addWidget(self.btn)

		# inizializza la maschera
		s1,s2 = self.lista[self.cntr]
		msg = "%s - %s" % (s1,s2)
		self.eLbl.setText(msg)

	def prev(self):
		if self.cntr > 0:
			self.cntr -= 1
			s1,s2 = self.lista[self.cntr]
			msg = "%s - %s" % (s1,s2)
			self.eLbl.setText(msg)

	def next(self):
		if self.cntr < self.nStaz-1:
			self.cntr += 1
			s1,s2 = self.lista[self.cntr]
			msg = "%s - %s" % (s1,s2)
			self.eLbl.setText(msg)

	def update(self):
		s1,s2 = self.lista[self.cntr]
#		print "cerco la stazione",s1
		for k,l1 in enumerate(self.libretto):
			tmp1 = l1.split("|")
#			print tmp0
			if tmp1[0] == '1' and tmp1[1] == s1:
#				print "cerco la stazione mirata",s2
				for j in range(k+1,len(self.libretto)):	# dalla stazione precedente alla fine del libretto
					l2 = self.libretto[j]
					tmp2 = l2.split("|")
#					print tmp1
					if tmp2[0] == '2' and tmp2[1] == s2:
						print(s1,':',k,l1)
						print(s2,':',j,l2)
						# chiama il dialogo (layout)
						hs = float(tmp1[2])
						av = float(tmp2[3])
						d = float(tmp2[4])
						h = float(tmp2[5])
#						print s1,hs,s2,av,d,h
						dlg = celerimensuraDlg(self.a_giro,s1,hs,s2,av,d,h)
						dlg.show()
						dlg.exec_()
						break

# ======================== classe principale ========================

class topog4qgis:
	vers = '0.2'
	build_date = '2020-29-04'
	author = 'giuliano curti (giulianc51@gmail.com)'
	contributor = 'giuseppe patti (gpatt@tiscali.it)'
	maintainer = 'marco lombardi (marco.lombardi.rm@gmail.com)'
	copyright = '2013-2014 giuliano curti'
	license = 'GPL  (http://www.gnu.org/licenses/gpl-2.0.html)'
	anTxtList = []		# lista delle annotazioni
	cLayer = ''
	provider = ''
	layGeom = ''
	eps = 0

	isClickToolActivated = False
	isFirst = True
	a_giro = 400.0	# angoli centesimali(400) o sessadecimali(360)

	def new(self):
		# ---- cancella i layer ------
		#s = QSettings()
		#crs = 'USER:000000'
		#oldValidation = s.value( "/Projections/defaultBehavior" )
		#s.setValue( "/Projections/defaultBehavior", "prompt" )
		#tmpLayer = QgsVectorLayer("Polygon?crs=EPSG:4326", "result", "memory")
		#QgsMapLayerRegistry.instance().addMapLayer(tmpLayer)
        
        # create Point layer
		#tmpLayer = QgsVectorLayer("Polygon?crs=EPSG:", "result", "memory")
		#QgsProject.instance().addMapLayer(tmpLayer, True)
        
        #QgsProject.instance().removeAllMapLayers()
		self.layLibMisur = '' # layer vertici rilievo (misurati)
		self.layLibRibat = '' # layer vertici rilievo (ribattuti)
		self.layLibCollim = '' # layer vertici rilievo (collimati)
		self.layLibCollimWGS84 = '' # layer vertici rilievo (collimati su WGS84)        
		self.layLibCtrn = '' # layer contorni rilievo
		self.layEdmPf = ''	# layer punti fiduciali
		self.layEdmVrts = ''	# layer vertici dell'EdM
		self.layEdmCtrn = ''	# layer contorni EdM
		# ----- azzera le variabili -------
		self.libretto = []	# copia del libretto in ram 
		self.misurati = []	# archivio dei vertici misurati [cod,x,y,z,note]
		self.ribattuti = []	# archivio dei vertici ribattuti
		self.collimati= []	# archivio dei punti collimati con i PF
		self.edmPf = []		# archivio PF
		self.edmVrts = []		# archivio vertici di particella dell'EdM
		# variabili globali necessarie per consentire il rigenerazione dei layer dopo georeferenzazione
		self.RilCtrn = []		# archivio dei contorni del libretto
		self.RilSty = []		# archivio degli dei contorni del libretto
		# NB edmCtrn,edmSty si possono lasciare locali
		# ------ disattiva i menu --------
		self.bImpEDM.setEnabled(True)
		self.bImpLib.setEnabled(True)
		self.bPfTaf.setDisabled(True)
		self.bViewLib.setEnabled(True)
		self.bPfRil.setDisabled(True)
		self.bDistPfRil.setDisabled(True)
		self.bVrtsPrtcEdm.setDisabled(True)
		self.bStazList.setDisabled(True)
		#self.bStazVrt.setDisabled(True)
		#self.bVrtsStaz.setDisabled(True)
		self.bMisurList.setDisabled(True)
		self.bRibatList.setDisabled(True)
		self.bCollimList.setDisabled(True)
		self.bNavPol.setDisabled(True)
		self.bDistRid.setDisabled(True)
		#self.bRAP.setDisabled(True)
		#self.bPAP.setDisabled(True)
		#self.bRRP.setDisabled(True)
		#self.bPRP.setDisabled(True)
		self.bGeoref.setDisabled(True)
		self.bErrPf.setDisabled(True)
		self.bDistPf.setDisabled(True)        
		self.bBaric.setDisabled(True)        
		self.bPfUff.setDisabled(True)
		self.bDistPfUff.setDisabled(True)
		self.bDistPfArch.setDisabled(True)

	def __init__(self, iface):
		# Save reference to the QGIS interface
		self.iface = iface
		# reference to map canvas
		self.canvas = self.iface.mapCanvas()
		# il rubber band è sempre attivo
		self.rubBnd = QgsRubberBand(self.canvas)
		self.rubBnd.setColor(QColor('#ff8800'))
		# create our GUI dialog
		self.dlg = QDialog()
		#self.dlg.setWindowFlags(Qt.CustomizeWindowHint)
#		con il clickTool qui non si riesce più a prendere il comando del mouse dopo aver eseguito un comando esterno
		# connect the layer changed handler to a signal that the TOC layer has changed
		self.iface.currentLayerChanged.connect(self.myHandleLayerChange)
		# out click tool will emit a QgsPoint on every click
		self.clickTool = QgsMapToolEmitPoint(self.canvas)
		# make our clickTool the tool that we'll use for now
		self.canvas.setMapTool(self.clickTool)

	def initGui(self):
		# Create action that will start plugin configuration
		self.action = QAction(QIcon(os.path.join(os.path.dirname(__file__), 'icon.png')),"topog4qgis",self.iface.mainWindow())
		#self.action = QtWidgets.QAction(QtGui.QIcon(os.path.join(os.path.dirname(__file__), 'icon.png')), self.dlg("Click to open"), self.iface.mainWindow())
		# connect the action to the run method
		self.action.triggered.connect(self.run)
		# Add toolbar button and menu item
		self.iface.addToolBarIcon(self.action)
		self.iface.addPluginToMenu("topog4qgis", self.action)
		self.dlg.setWindowTitle("topog4qgis v0.2")
        # -------- file menubar ------------
		mb = QMenuBar(self.dlg)
		mb.setGeometry(0,0,270,120)

		# ---------- file menu -----------------
		mFile = mb.addMenu('File')

		tmp = QAction(QIcon(''),'Nuovo',self.dlg)        
		tmp.triggered.connect(self.new)
		mFile.addAction(tmp)

		#tmp = QAction(QIcon(''),'Importa Libretto',self.dlg)
		self.bImpLib = QAction(QIcon(''),'Importa Libretto',self.dlg)
		self.bImpLib.triggered.connect(self.importaLibretto)
		#mFile.addAction(tmp)
		mFile.addAction(self.bImpLib)
		self.bImpLib.setDisabled(True)

		#tmp = QAction(QIcon(''),'Importa EdM',self.dlg)
		self.bImpEDM = QAction(QIcon(''),'Importa EdM',self.dlg)        
		self.bImpEDM.triggered.connect(self.importaEDM)
		#mFile.addAction(tmp)
		mFile.addAction(self.bImpEDM)
		self.bImpEDM.setDisabled(True)

		self.bPfTaf = QAction(QIcon(''),'Importa PF da Taf',self.dlg)        
		self.bPfTaf.triggered.connect(self.importaPfDaTaf)
		mFile.addAction(self.bPfTaf)
		self.bPfTaf.setDisabled(True)

		# ---------- referencng menu --------------
		mGeoref = mb.addMenu('Elaborazione')

		self.bGeoref = QAction(QIcon(''),'Rototrasla su PF',self.dlg)        
		self.bGeoref.triggered.connect(self.georeferencer)
		mGeoref.addAction(self.bGeoref)
		self.bGeoref.setDisabled(True)
        
		self.bGeorefWGS84 = QAction(QIcon(''),'Georeferenzia su WGS84(EPSG:4326)',self.dlg)        
		self.bGeorefWGS84.triggered.connect(self.georeferWGS84)
		mGeoref.addAction(self.bGeorefWGS84)
		self.bGeorefWGS84.setDisabled(True)        

		# -------- validation menu -------------
		mValid = mb.addMenu('Calcoli')

		self.bErrPf = QAction(QIcon(''),'Errore medio PF',self.dlg)        
		self.bErrPf.triggered.connect(self.errorePF)
		mValid.addAction(self.bErrPf)
		self.bErrPf.setDisabled(True)
        
		self.bDistPf = QAction(QIcon(''),'Distanze PF',self.dlg)        
		self.bDistPf.triggered.connect(self.distanzePF)
		mValid.addAction(self.bDistPf)
		self.bDistPf.setDisabled(True)        
        
		self.bBaric = QAction(QIcon(''),'Baricentri Geometrici',self.dlg)        
		self.bBaric.triggered.connect(self.baricentro)
		mValid.addAction(self.bBaric)
		self.bBaric.setDisabled(True)

		# ---------- inquiry menu --------------
		mInquiry = mb.addMenu('Funzioni')

		self.bViewLib = QAction(QIcon(''),'Vedi Libretto',self.dlg)        
		self.bViewLib.triggered.connect(self.librViewTool)
		mInquiry.addAction(self.bViewLib)
		self.bViewLib.setDisabled(True)

		self.bPfUff = QAction(QIcon(''),'Elenco PF ufficiali',self.dlg)        
		self.bPfUff.triggered.connect(self.elencoPfUfficiali)
		mInquiry.addAction(self.bPfUff)
		self.bPfUff.setDisabled(True)

		self.bDistPfUff = QAction(QIcon(''),'Distanze PF ufficiali',self.dlg)        
		self.bDistPfUff.triggered.connect(self.distPfUfficiali)
		mInquiry.addAction(self.bDistPfUff)
		self.bDistPfUff.setDisabled(True)

		self.bDistPfArch = QAction(QIcon(''),'Importa archivio distanze PF',self.dlg)        
		self.bDistPfArch.triggered.connect(self.leggiDis)
		mInquiry.addAction(self.bDistPfArch)
		self.bDistPfArch.setDisabled(True)

		self.bPfRil = QAction(QIcon(''),'Elenco PF rilevati',self.dlg)        
		self.bPfRil.triggered.connect(self.elencoPfRilevati)
		mInquiry.addAction(self.bPfRil)
		self.bPfRil.setDisabled(True)

		self.bDistPfRil = QAction(QIcon(''),'Distanze PF misurate',self.dlg)        
		self.bDistPfRil.triggered.connect(self.distPfMisurate)
		mInquiry.addAction(self.bDistPfRil)
		self.bDistPfRil.setDisabled(True)

		self.bVrtsPrtcEdm = QAction(QIcon(''),'Vertici particella EdM',self.dlg)        
		self.bVrtsPrtcEdm.triggered.connect(self.elencoEdmVertici)
		mInquiry.addAction(self.bVrtsPrtcEdm)
		self.bVrtsPrtcEdm.setDisabled(True)

		self.bMisurList = QAction(QIcon(''),'Elenco misurati',self.dlg)        
		self.bMisurList.triggered.connect(self.elencoMisurati)
		mInquiry.addAction(self.bMisurList)
		self.bMisurList.setDisabled(True)

		self.bCollimList = QAction(QIcon(''),'Elenco collimati',self.dlg)        
		self.bCollimList.triggered.connect(self.elencoCollimati)
		mInquiry.addAction(self.bCollimList)
		self.bCollimList.setDisabled(True)
		
		self.bStazList = QAction(QIcon(''),'Elenco stazioni',self.dlg)        
		self.bStazList.triggered.connect(self.elencoStazioni)
		mInquiry.addAction(self.bStazList)
		self.bStazList.setDisabled(True)

		#self.bStazVrt = QAction(QIcon(''),'Stazioni su Vertice',self.dlg)        
		#self.bStazVrt.triggered.connect(self.stazioniSuVerticeTool)
		#mInquiry.addAction(self.bStazVrt)
		#self.bStazVrt.setDisabled(True)

		#self.bVrtsStaz = QAction(QIcon(''),'Vertici da stazione',self.dlg)        
		#self.bVrtsStaz.triggered.connect(self.verticiDaStazioneTool)
		#mInquiry.addAction(self.bVrtsStaz)
		#self.bVrtsStaz.setDisabled(True)

		self.bRibatList = QAction(QIcon(''),'Elenco ribattuti',self.dlg)        
		self.bRibatList.triggered.connect(self.elencoRibattuti)
		mInquiry.addAction(self.bRibatList)
		self.bRibatList.setDisabled(True)

		self.bNavPol = QAction(QIcon(''),'Navigazione poligonale',self.dlg)        
		self.bNavPol.triggered.connect(self.navigazioneStazioni)
		mInquiry.addAction(self.bNavPol)
		self.bNavPol.setDisabled(True)

		self.bDistRid = QAction(QIcon(''),'Elenco distanze ridotte',self.dlg)        
		self.bDistRid.triggered.connect(self.elencoDistRidotte)
		mInquiry.addAction(self.bDistRid)
		self.bDistRid.setDisabled(True)

		self.bOssCeler = QAction(QIcon(''),'Elenco osservazioni celerimetriche',self.dlg)        
		self.bOssCeler.triggered.connect(self.elencoOssCeler)
		mInquiry.addAction(self.bOssCeler)
		self.bOssCeler.setDisabled(True)

		#self.bRAP = QAction(QIcon(''),'Rectangular absolute position',self.dlg)        
		#self.bRAP.triggered.connect(self.rectAbsPosTool)
		#mInquiry.addAction(self.bRAP)
		#self.bRAP.setDisabled(True)

		#self.bPAP = QAction(QIcon(''),'Polar absolute position',self.dlg)        
		#self.bPAP.triggered.connect(self.polAbsPosTool)
		#mInquiry.addAction(self.bPAP)
		#self.bPAP.setDisabled(True)

		#self.bRRP = QAction(QIcon(''),'Rectangular relative position',self.dlg)        
		#self.bRRP.triggered.connect(self.rectRelPosTool)
		#mInquiry.addAction(self.bRRP)
		#self.bRRP.setDisabled(True)

		#self.bPRP = QAction(QIcon(''),'Polar relative position',self.dlg)        
		#self.bPRP.triggered.connect(self.polRelPosTool)
		#mInquiry.addAction(self.bPRP)
		#self.bPRP.setDisabled(True)

		# ----------- help menu ----------------
		mHelp = mb.addMenu('Help')

		tmp = QAction(QIcon(''),'About',self.dlg)        
		tmp.triggered.connect(self.about)
		mHelp.addAction(tmp)

		tmp = QAction(QIcon(''),'Info',self.dlg)        
		tmp.triggered.connect(self.info)
		mHelp.addAction(tmp)

		# ----- control panel-------
		# ------ angle notation ------
		tmp = QLabel(self.dlg)
		tmp.setGeometry(QRect(10,30,120,20))
		tmp.setText('angle:')

		self.bAngle = QPushButton('%d' % (self.a_giro),self.dlg)
		self.bAngle.clicked.connect(self.switchAngle)
		self.bAngle.setGeometry(QRect(110,30,120,20))

		# ------ nome layer ------
		tmp = QLabel(self.dlg)
		tmp.setGeometry(QRect(10,60,120,20))
		tmp.setText('layer:')

		self.eLay = QLabel(self.dlg)
		self.eLay.setGeometry(QRect(110,60,120,20))
		self.eLay.setText('')

		# ------ nome layer ------
		#tmp = QLabel(self.dlg)
		#tmp.setGeometry(QRect(10,70,120,20))
		#tmp.setText('type:')

		#self.eType = QLabel(self.dlg)
		#self.eType.setGeometry(QRect(110,70,120,20))
		#self.eType.setText('')

		# -------- editabilità  ----------
		#tmp = QLabel(self.dlg)
		#tmp.setGeometry(QRect(10,90,120,20))
		#tmp.setText('editable:')

		#self.eLayEdit = QLabel(self.dlg)
		#self.eLayEdit.setGeometry(QRect(110,90,120,20))
		#self.eLayEdit.setText('')

		# ---- provider -----
		tmp = QLabel(self.dlg)
		tmp.setGeometry(QRect(10,90,120,20))
		tmp.setText('provider:')

		self.eProv = QLabel(self.dlg)
		self.eProv.setGeometry(QRect(110,90,120,20))
		self.eProv.setText('')

		# ---- tolleranza ricerca ----
		#tmp = QLabel(self.dlg)
		#tmp.setGeometry(QRect(10,130,120,20))
		#tmp.setText('tolerance:')

		#self.eEps = QLineEdit(self.dlg)
		#self.eEps.setGeometry(QRect(110,130,120,20))
		#self.eEps.setText('')

		# inizializza le variabili
		#self.new()

	# run method that performs all the real work
	def run(self):
#		con la gestione del clickTool quì, da errore quando si attiva un comando esterno al plugin
		# connect the layer changed handler to a signal that the TOC layer has changed
#		QObject.connect(self.iface,SIGNAL("currentLayerChanged(QgsMapLayer *)"), self.myHandleLayerChange)
		# out click tool will emit a QgsPoint on every click
#		self.clickTool = QgsMapToolEmitPoint(self.canvas)
		# make our clickTool the tool that we'll use for now
#		self.canvas.setMapTool(self.clickTool)
		# show the dialog
		self.dlg.show()
		result = self.dlg.exec_()

	def unload(self):
		# rimuove i segnali attivati
		#if self.cLayer:	# potrebbero non esserci layer attivi
			#self.cLayer.removeSelection()
			# disconnect the layer changed handler
			#self.iface.currentLayerChanged.disconnect(self.currentLayerChanged)
		# try to disconnect all signals
		if self.isClickToolActivated:
			self.clickTool.canvasClicked.disconnect()
			self.isClickToolActivated = False
		# Remove the plugin menu item and icon
		#QgsApplication.processingRegistry().removeProvider(self.provider)
		self.iface.removePluginMenu("topog4qgis",self.action)
		self.iface.removeToolBarIcon(self.action)

#	------------ utility functions -----------

	def about(self):
		about(self.iface.mainWindow(),self)

	def info(self):
		info(self.iface.mainWindow(),self.vers,self.author)

	def switchAngle(self):
		"""
			switcha la notazione angolare ed aggiorna la descrizione del menu
			NB: quest'ultima però non funziona !
		"""
		if self.a_giro == 400.0:
			self.a_giro = 360.0
		else:
			self.a_giro = 400.0
		self.bAngle.setText('%d' % (self.a_giro))

#	-------- managing layer functions -----------

	def myHandleLayerChange(self):
		"""
			versione (pesantemente) modificata

			aggiorna:
				- layer corrente
				- provider
				- tipo
				- editabilità
		"""
		# controlla il layer 
		if self.iface.activeLayer():
			# registra il layer corrente
			self.cLayer = self.iface.activeLayer()
			self.eLay.setText(self.cLayer.name())
			#self.eLayEdit = self.cLayer.isEditable()
			# registra il provider
			self.provider = self.cLayer.dataProvider()
			self.eProv.setText(self.provider.name())           
			# registra se è vettoriale
			#self.layGeom = self.cLayer.type()
			#self.eType.setText(str(self.layGeom))
			# al cambio del layer stabilisce la tolleranza di ricerca
			#self.eps = self.canvas.mapUnitsPerPixel() * 5
			#self.eEps.setText("%7.2f" % (self.eps))
		else:
			self.cLayer = ''
			self.eLay.setText('')
			#self.layGeom = ''
			#self.eLayEdit = ''
			self.provider = ''
			self.eProv.setText('')
			#self.eType.setText('')
			#self.eps = 0
			#self.eEps.setText('')

	def newLayer(self,layType,layName,attrLst):
		# create Point layer
		self.cLayer = QgsVectorLayer(layType,layName,"memory")
		self.provider = self.cLayer.dataProvider()
		# add fields
		for a in attrLst:
			self.provider.addAttributes([QgsField(a[0],a[1])])
		# aggiunge al registry
		QgsProject.instance().addMapLayer(self.cLayer, True)

	def creaPointLayer(self,title,attrs,vrts):
		"""
			crea un layer point
			- title	titolo del layer
			- attrs	attributi
			- vrts	elenco vertici
						i[0]	indice
						i[3+i]	attributo i.mo
		"""
		self.newLayer('Point',title,attrs)	# lavora sulla variabile globale self.cLayer
		nAttrs=len(attrs)
		# Enter editing mode
		self.cLayer.startEditing()
		# parsing dei vertici
		feat = QgsFeature()
		for i in vrts:            
			if len(i) < nAttrs:	# deve avere gli n attributi + la coppia x,y
				print('registrazione',i,'errata')
			else:            
				feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(i[1],i[2])))
				feat.initAttributes(nAttrs)
				feat.setAttribute(0,i[0])                
				for k in range(nAttrs):                 
					feat.setAttribute(k,i[k])
				self.cLayer.addFeatures([feat])
		# Commit changes
		self.cLayer.commitChanges()
		# attiva la label
		text_format = QgsTextFormat()
		label = QgsPalLayerSettings()
		label.fieldName = 'indice'
		label.enabled = True
		label.setFormat(text_format)
		label.placement = QgsPalLayerSettings.Line
		labeler = QgsVectorLayerSimpleLabeling(label)
		self.cLayer.setLabeling(labeler)
		self.cLayer.triggerRepaint()
		# rinfresca il video
		self.canvas.refresh()
		# update layer's extent
		self.cLayer.updateExtents()
		# cancella la selezione (che non so come sia attivata)
		self.cLayer.removeSelection()

	def creaLineLayer(self,title,lines,styles,points):
		"""
			crea un layer linestring (attributi sono fissi: indice e tratto)
			- title	titolo del layer
			- lines	contorni
			- styles	stili linea
			- points	elenco vertici
						i[0]	indice
						i[3+i]	attributo i.mo
		"""
		self.newLayer("Linestring",title,[["indice",QVariant.Int],['TRATTO',QVariant.String]])
		# Enter editing mode
		self.cLayer.startEditing()
		for i,p in enumerate(lines):
#			print p
			list = []
			for v in p:
				pos,x,y,z = pointArchivioCds(points,v)
				list.append(QgsPoint(x,y))
			feat = QgsFeature()
			feat.setGeometry(QgsGeometry.fromPolyline(list))
			feat.initAttributes(2)
			feat.setAttribute(0,str(p))
			feat.setAttribute(1,styles[i])
			self.cLayer.addFeatures([feat])
		# Commit changes
		self.cLayer.commitChanges()
		# rinfresca il video
		self.canvas.refresh()
		# cancella la selezione (che non so come sia attivata)
		self.cLayer.removeSelection()

#	--------------------- I/O functions -----------

	def importaLibretto(self):
		root = QgsProject.instance().layerTreeRoot()
		msg = ''
		self.iface.messageBar().pushMessage(
			"importaLibretto",
			"Stai usando la notazione angolare %d, assicurati che sia quella corretta" % (int(self.a_giro)),
			level=Qgis.Info
		)
		# ---------- carica il libretto delle misure -----------
		fname = QFileDialog.getOpenFileName(self.iface.mainWindow(),'Open file','~','*.dat')
		if fname[0] != "":
			isCel = False
			isGps = False
			self.libretto = loadFile(fname[0])
			print("Lette %d registrazioni" % len((self.libretto)))
			# legge registrazioni gps
			myGps = openLibretto_gps(self.libretto)
			print("Trovate %d registrazioni GPS" % (len(myGps)))
			# legge vertici del libretto
			RilVrts = openLibretto_vertici(self.libretto)
			print("Trovate %d stazioni celerimetriche" % (len(RilVrts)))
			# legge gli allineamenti e squadri
			RilAll = openLibretto_allSquad(self.libretto)
			print("Trovati %d allineamenti e squadri" % (len(RilAll)))
			# legge contorni del rilievo
			self.RilCtrn,self.RilSty = openLibretto_contorni(self.libretto)
			print("trovati %d contorni nel libretto" % (len(self.RilCtrn)))
			# ----------- tratta il libretto ------------------------
			# in questa versione trattiamo prioritariamente le letture gps;
			# se ci sono divengono lo spazio iniziale; in caso negativo allo scopo
			# viene destinata la prima stazione celerimetrica;
			rilievo = []	# inizializza l'archivio
			# ----- elaborazione dei gps ---------
			if myGps:
				tmp = elaborazioneGps(myGps)
				for g in tmp:
					rilievo.append(g)
				isFirst = False	# azzera il flag
				isGps = True
			else:
				isFirst = True
				isGps = False
			# ----- elaborazione celerimetriche ---------
			if RilVrts:
				isCel = True
			if isGps == True and isCel == False:
				tipologia = 0 #gps
			elif isGps == False and isCel == True:
				tipologia = 1 #tps
			elif isGps == True and isCel == True:
				tipologia = 2 #misto
			for s in RilVrts:
				# conversione in coordinate rettangolari
				tmp = polar2rect(s,self.a_giro)
				if not isFirst:
					# collimazione (tranne la prima stazione)
					tmp = collimazioneStazione(tmp,rilievo,tipologia)
# occorre controllare i valori di ritorno, ci potrebbero essere stazioni senza punti ribattuti
				isFirst = False
				if len(tmp):
#					print "trasferimento in archivio della stazione",tmp[0][0]
					for i in tmp:
						rilievo.append(i)
				else:
					print("Errore: polar2rect non restituisce vertici; il rilievo è:")
					for i in rilievo:
						print(i)
			# ----- elaborazione allineamenti ---------
			# (deve essere gestito iterativamente perchè un allineamento può operare
			# su allineamenti precedenti)
			# NB: gli allineamenti vanno eseguiti prima della collimazione finale con i fiduciali
			# altrimenti i punti terminali risultano sensibilmente diversi dalla loro definizione
			# celerimetrica a causa di possibili scalature nella fase di collimazione ai PF;
			for all in RilAll:
				tmp = allEsquad(rilievo,all)
				rilievo.append(tmp)
#			duplica tutto nei misurati
			self.misurati = copy.deepcopy(rilievo)	# forza la clonazione
			# ---------- separa i ribattuti dal resto dell'archivio ---------------
			i = 0
			maxPnts = len(self.misurati)
			while i < maxPnts:
				p = self.misurati[i]
				k = i+1
				while k < maxPnts:
					q = self.misurati[k]
#					print "controllo",i,p[0],k,q[0],maxPnts
					if p[0] == q[0]:
#						print "punto ribattuto"
						tmp = self.misurati.pop(k)	# toglie da misurati
						self.ribattuti.append(tmp)	# e mette in ribattuti
						maxPnts -= 1				# aggiorna il numero dei misurati
					k += 1
				i += 1
			# ---------parte grafica --------------
			# crea layer vertici misurati
			if len(self.misurati):
				self.creaPointLayer('Rilievo_vertici_misurati',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],self.misurati)
				self.layLibMisur = self.cLayer
				self.cLayer.setLabelsEnabled(True)
				print("Layer vertici misurati completato")
			else:
				print("non ci sono vertici misurati nel libretto")
			# crea layer dei ribattuti
			#if len(self.ribattuti):
			#	self.creaPointLayer('Rilievo_vertici_ribattuti',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],self.ribattuti)
			#	self.layLibRibat = self.cLayer
			#	self.cLayer.setLabelsEnabled(True)                
			#	print("Layer vertici ribattuti completato")
			#else:
			#	print("non ci sono vertici ribattuti nel libretto")
			# crea layer contorni
			if len(self.RilCtrn):
				self.creaLineLayer('Rilievo_contorni',self.RilCtrn,self.RilSty,self.misurati)
				self.layLibCtrn = self.cLayer
				# ------ attiva la simbologia categorizzata per i contorni -------
				myRen = catSymbol(
					self.cLayer.geometryType(),
					'TRATTO',
					[
						['NC','#000000','nero continuo'],
						['RC','#ff0000','rosso continuo'],
						['NT','#000000','nero tratteggiato'],
						['RT','#ff0000','rosso tratteggiato'],
						['NP','#000000','nero puntinato'],
						['RP','#ff0000','rosso puntinato']                       
					]
				)
				self.cLayer.setRenderer(myRen)
				print("Layer contorni libretto completati")
			else:
				print("non ci sono contorni nel libretto")
			# attiva le voci di menu
			self.bImpEDM.setEnabled(True)
			self.bImpLib.setEnabled(True)
			self.bPfTaf.setEnabled(True)
			self.bViewLib.setEnabled(True)
			self.bPfRil.setEnabled(True)
			self.bDistPfRil.setEnabled(True)
			self.bDistPfArch.setDisabled(True)
			self.bMisurList.setEnabled(True)
			if isGps == True:
				self.bStazList.setDisabled(True)
				self.bGeorefWGS84.setEnabled(True)
			else:
				self.bStazList.setEnabled(True)
			#self.bStazVrt.setEnabled(True)
			#self.bVrtsStaz.setEnabled(True)
			if isGps == True:
				self.bCollimList.setDisabled(True)
				self.bRibatList.setDisabled(True)
				self.bNavPol.setDisabled(True)
				self.bDistRid.setDisabled(True)
				self.bOssCeler.setDisabled(True)
			else:
				self.bCollimList.setDisabled(True)
				self.bRibatList.setEnabled(True)
				self.bNavPol.setEnabled(True)
				self.bDistRid.setEnabled(True)
				self.bOssCeler.setEnabled(True)
			if isGps == True and isCel == True:
				self.bCollimList.setDisabled(True)
				self.bRibatList.setEnabled(True)
				self.bNavPol.setEnabled(True)
				self.bDistRid.setEnabled(True)
				self.bOssCeler.setEnabled(True)
				self.bStazList.setEnabled(True)
				self.bGeorefWGS84.setEnabled(True)                
			#self.bRAP.setEnabled(True)
			#self.bPAP.setEnabled(True)
			#self.bRRP.setEnabled(True)
			#self.bPRP.setEnabled(True)
		else:
			self.iface.messageBar().pushMessage(
				"importaLibretto",
				"Errore: occorre un nome di file valido",
				level=Qgis.Warning,
				duration = 4
			)

	def importaEDM(self):
		"""
			Apre un file EdM
			gestire il colore/tratto dei segmenti
		"""
		root = QgsProject.instance().layerTreeRoot()
		#groupEDM = root.addGroup("EdM")
		if len(self.edmPf):
			result = QMessageBox.question(
				self.iface.mainWindow(),
				"importaEDM",
				"I punti fiduciali sono gia disponibili, vuoi sovrascriverli?",
				QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
			if result == QMessageBox.No:
				return
			else:
				# cancella il layer dei PF
				QgsProject.instance().removeMapLayer(self.layEdmPf)
		# dialogo selezione file
		fname = QFileDialog.getOpenFileName(self.iface.mainWindow(),'Open file','~','*.emp')
		if fname[0] != "":
			edm = loadFile(fname[0])
			print('Lette %d registrazioni' % len(edm))
			# ---- legge punti fiduciali dell'EdM ---------
			self.edmPf = openEDM_pf(edm)
			print('Lette %d punti fiduciali' % len(self.edmPf))
			# ---- legge vertici dell'EdM --------
			self.edmVrts = openEDM_vertici(edm)
			print('Letti %d vertici nell EdM' % (len(self.edmVrts)))
			# ------ legge contorni particelle/fabbricati dell'EdM -------
			EdmPrtc,EdMSty = openLibretto_contorni(edm)
			print('Letti %d contorni dell EdM' % len(EdmPrtc))
			# ----- parte grafica ----------------
			# crea layer PF
			if len(self.edmPf):
				self.creaPointLayer('EdM_pf',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double]],self.edmPf)
				self.layEdmPf = self.cLayer
				print('Layer punti fiduciali completato')
				self.cLayer.setLabelsEnabled(True)
				#child0 = root.children()[0]
				#tmpLayer = child0.clone()
				#print (tmpLayer.name())
				#groupEDM.insertChildNode(0, tmpLayer)
			else:
				print('Nessun punto fiduciale nell EdM')
			# --------- scrittura vertici della particella ----------
			# queste sono superflue, possono servire solo alla numerazione dei vertici
			if len(self.edmVrts):
				self.creaPointLayer('EdM_vertici',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double]],self.edmVrts)
				print('Layer vertici particella completato')             
				#child0 = root.children()[0]
				#tmpLayer = child0.clone()
				#print (tmpLayer.name())
				#groupEDM.insertChildNode(0, tmpLayer)
				self.layEdmVrts = self.cLayer
			else:
				print('Nessun vertice nell EdM')
			# --------- scrittura contorni della particella ----------
			if len(EdmPrtc):
				self.creaLineLayer('EdM_contorni',EdmPrtc,EdMSty,self.edmVrts)
				# setta la simbologia
				myRen = singleSymbol(self.cLayer.geometryType(),'#000000')
				#verificare self.cLayer.setRenderer(myRen)
				print('Layer contorni particella completato')
				#print(EdmPrtc)
				#child0 = root.children()[0]
				#tmpLayer = child0.clone()
				#print (tmpLayer.name())
				#groupEDM.insertChildNode(0, tmpLayer)
				self.layEdmCtrn = self.cLayer
			else:
				print('Nessun contorno nell EdM')
            # attiva i menu
			self.bGeoref.setEnabled(True)
			self.bPfUff.setEnabled(True)
			self.bDistPfUff.setEnabled(True)
			self.bDistPfArch.setEnabled(True)
			self.bVrtsPrtcEdm.setEnabled(True)
			#self.bRAP.setEnabled(True)
			#self.bPAP.setEnabled(True)
			#self.bRRP.setEnabled(True)
			#self.bPRP.setEnabled(True)
		else:
			self.iface.messageBar().pushMessage(
				"openEDM",
				"Errore: occorre un nome di file valido",
				level=Qgis.Warning,
				duration=5
			)
        #QgsProject.instance().removeMapLayer(self.layEdmPf)

	def importaPfDaTaf(self):
		"""
			i dati sono formattati in modo diverso nel libretto/EdM e Taf;
			per confrontarli componiamo i dati della taf in modo congruente a libretto
			e EdM, ovviamente si poteva anche pulire in dati del libretto/EdM per
			estrarre il numero di fiduciale (forse più veloce)
		"""
		if len(self.edmPf):
			result = QMessageBox.question(
				self.iface.mainWindow(),
				"importaPfDaTaf",
				"I punti fiduciali sono già disponibili, vuoi sovrascriverli?",
				QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
			if result == QMessageBox.No:
				return
			else:
				# cancella il layer dei PF
				QgsProject.instance().removeMapLayer(self.layEdmPf)
		# cerca i PF
		pfList = pfLista(self.libretto)
		if pfList:
			lastfgCod = ""
			n = 0
			# legge file
			fname = QFileDialog.getOpenFileName(self.iface.mainWindow(),'Open file','~','*.taf')
			for i in range(-1,len(pfList)):
				pfCod,fgCod,comCod = pfList[i].split('/')	# prende uno qualsiasi, il primo
				fgCod = fgCod[0:-1]	# elimina ultimo carattere
				try:
					tmppfCod,nextfgCod,tmpcomCod = pfList[i+1].split('/')	# prende il prossimo                
					nextfgCod = nextfgCod[0:-1]	# elimina ultimo carattere
				except IndexError:
					break                                       
				#print("fgCod",fgCod)                    
				#print("lastfgCod",lastfgCod)
				#print("nextfgCod",nextfgCod)
                #---
				#if lastfgCod == "":
				#	lastfgCod = fgCod
				#	n = 1
				#	print("sto qui 1", n)
				#	continue                    
                #---    
				if nextfgCod == fgCod and fgCod == lastfgCod and n > 1:
					#print("sto qui 2", n)
					continue
				else: n = n+1
				if nextfgCod == fgCod and n > 1:
					#print("sto qui 3", n)                
					continue
				else: n = n+1
				print("Cerco i PF del Foglio",fgCod,"del Comune",comCod)
				if fname[0] != "":
					# legge file
					import codecs
					f = codecs.open(fname[0], 'r', encoding='utf-8', errors='ignore')
					for data in f:
						if fgCod[0] == "A":
							fgCod = "10" + fgCod[1:3]
						if fgCod[0] == "B":
							fgCod = "11" + fgCod[1:3]
						# pulisce la riga
						data = data.rstrip('\n')
						data = data.rstrip('\r')
						len_comCod =(len(comCod))
						if data[0:len_comCod] == comCod:
							if int(data[6:10]) == int(fgCod):
								#print('matcha il foglio',int(fgCod))
								tmp = int(data[15:17])
								fgAll = str(data[11:12])                           
								if fgAll == " ":
									fgAll = "0"
								lastfgCod = fgCod                                    
								#print('matcha il foglio',int(fgCod),'allegato',fgAll)                            
								if fgCod[0:2] == "10":
									fgCod = "A" + fgCod[2:4]
									lastfgCod = fgCod                                    
									#print(fgCod)
								if fgCod[0:2] == "11":
									fgCod = "B" + fgCod[2:4] 
									lastfgCod = fgCod                                    
									#print(fgCod)                                
								if tmp < 10:
									txt = "PF0%1s/%3s%1s/" % (tmp,fgCod,fgAll) + comCod
								else:
									txt = "PF%2s/%3s%1s/" % (tmp,fgCod,fgAll) + comCod
								#print(txt)
								if txt in pfList:	# matcha anche il pf
									y,x = data[102:114],data[115:127]
									print("Trovato",txt,x,y)
									self.edmPf.append([txt,float(x),float(y),0])
							#lastfgCod = fgCod				
				else:
					self.iface.messageBar().pushMessage(
					"importaPfDaTaf",
					"Devi darmi un nome di file valido",
					level=Qgis.Warning,
					duration=4
				)
		# crea layer PF
			if len(self.edmPf):
				self.creaPointLayer('EdM_pf',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double]],self.edmPf)
				self.layEdmPf = self.cLayer
				self.cLayer.setLabelsEnabled(True)
				print('Layer punti fiduciali completato')
				# attiva i menu
				self.bGeoref.setEnabled(True)
				self.bPfUff.setEnabled(True)
				self.bDistPfUff.setEnabled(True)
				#self.bRAP.setEnabled(True)
				#self.bPAP.setEnabled(True)
				#self.bRRP.setEnabled(True)
				#self.bPRP.setEnabled(True)
		else:
			self.iface.messageBar().pushMessage(
				"importaPfDaTaf",
				"Attenzione: il rilievo non sembra contenere punti fiduciali",
				level=Qgis.Warning,
				duration=4
			)
            
#	------------ calcolo baricentro  --------------------

	def baricentro(self):
		listaPf = pfLista(self.libretto)
		for i0,j0 in enumerate(listaPf):
			p,x0,y0,z = pointArchivioCds(self.misurati,j0)
			p,e0,n0,z = pointArchivioCds(self.edmPf,j0)
			for i1 in range(i0,1):
				j1 = listaPf[i1]
				p,x1,y1,z = pointArchivioCds(self.misurati,j1)
				p,e1,n1,z = pointArchivioCds(self.edmPf,j1) 
			for i2 in range(i0,2):
				j2 = listaPf[i2]
				p,x2,y2,z = pointArchivioCds(self.misurati,j2)
				p,e2,n2,z = pointArchivioCds(self.edmPf,j2)                    
		# stampa     
		#print('%s %12.3f %12.3f' % (i0,e0,n0))
		#print('%s %12.3f %12.3f' % (i1,e1,n1))
		#print('%s %12.3f %12.3f' % (i2,e2,n2))  
		barPFuff = []
		barPFuff.append((e0+e1+e2)/3)
		barPFuff.append((n0+n1+n2)/3)
		barPFuff.append(0)
		print("baricentro PF (TAF)",barPFuff)
        
		#print('%s %12.3f %12.3f' % (i0,x0,y0))
		#print('%s %12.3f %12.3f' % (i1,x1,y1))
		#print('%s %12.3f %12.3f' % (i2,x2,y2))
		barPFcol = []
		barPFcol.append((x0+x1+x2)/3)
		barPFcol.append((y0+y1+y2)/3)
		barPFcol.append(0)
		print("baricentro PF (collimati)",barPFcol)        

		deltaX = ((e0+e1+e2)/3)-((x0+x1+x2)/3)
		deltaY = ((n0+n1+n2)/3)-((y0+y1+y2)/3)        
		print("deltaX",deltaX)
		print("deltaY",deltaY)       
		return barPFcol,barPFuff, deltaX, deltaY
            
#	---------- referencing functions --------------------

	def georeferWGS84(self):
		archivio = []
		if self.libretto[1][0] == '6':
			tmp0 = self.libretto[2]
			k = 4            
		else:
			tmp0 = self.libretto[1]
			k = 3                    
		id0 = tmp0.split('|')[1] 
		baseline = tmp0.split('|')[2]        
		x0,y0,z0 = float(baseline.split(',')[0]),float(baseline.split(',')[1]),float(baseline.split(',')[2])        
		lon0,lat0,h0 = geocentriche2wgs84(x0,y0,z0) 
		#print(id0,math.degrees(lon0),math.degrees(lat0),h0)        
		archivio.append([id0,math.degrees(lon0),math.degrees(lat0),h0,'gps',id0,1])        
		for i in range(k,len(self.libretto)):
			tmp1 = self.libretto[i]
			if tmp1[0] == '2':
				id1 = tmp1.split('|')[1] 
				coordinate = tmp1.split('|')[2]  
				dx,dy,dz = float(coordinate.split(',')[0]),float(coordinate.split(',')[1]),float(coordinate.split(',')[2])                
				lon1,lat1,h1 = geocentriche2wgs84(x0+dx,y0+dy,z0+dz) 
				#print(id1,math.degrees(lon1),math.degrees(lat1),h1)                                
				archivio.append([id1,math.degrees(lon1),math.degrees(lat1),h1,'gps',id0,i])                
			else:      
				break 
		RilVrts = openLibretto_vertici(self.libretto)                
		if RilVrts:
			print("Trovate %d stazioni celerimetriche" % (len(RilVrts)))
			for i in range(0,len(self.misurati)):
				tmp1 = self.misurati[i]
				if tmp1[4] != 'gps':
					id1,u,v,w = tmp1[0],tmp1[1],tmp1[2],tmp1[3]
					x,y,z = topocentriche2geocentriche(u,v,w,x0,y0,z0)                    
					lon1,lat1,h1 = geocentriche2wgs84(x,y,z) 
					#print(id1,math.degrees(lon1),math.degrees(lat1),h1)                    
					archivio.append([id1,math.degrees(lon1),math.degrees(lat1),h1,'gps',id0,i])                    
		QgsProject.instance().layerTreeRoot().findLayer(self.layLibMisur.id()).setItemVisibilityChecked(False)                
		self.creaPointLayer('Rilievo_vertici_collimati_WGS84',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],archivio)
		self.cLayer.setLabelsEnabled(True)
		self.layLibCollimWGS84 = self.cLayer
		print('Layer vertici collimati su WGS84 (EPSG:4326) completato')
		# crea layer contorni
		if len(self.RilCtrn):
			QgsProject.instance().layerTreeRoot().findLayer(self.layLibCtrn.id()).setItemVisibilityChecked(False)         
			self.creaLineLayer('Rilievo_contorni_WGS84',self.RilCtrn,self.RilSty,archivio)
			self.layLibCtrn = self.cLayer
			# ------ attiva la simbologia categorizzata per i contorni -------
			myRen = catSymbol(
				self.cLayer.geometryType(),
				'TRATTO',
				[
					['NC','#000000','nero continuo'],
					['RC','#ff0000','rosso continuo'],
					['NT','#000000','nero tratteggiato'],
					['RT','#ff0000','rosso tratteggiato'],
					['NP','#000000','nero puntinato'],
					['RP','#ff0000','rosso puntinato']                       
				]
				)
			self.cLayer.setRenderer(myRen)
			print("Layer contorni libretto completati")
		else:
			print("non ci sono contorni nel libretto")        
        
	def georeferencer(self):
		"""
			esegue la rototraslazione con >=2 PF
			e la collimazione con >=3
		"""
		# controlla i PF misurati nel rilievo
		listaPf = pfLista(self.libretto)
		print("Rilevati i PF",listaPf)
		if len(listaPf) < 2:
			self.iface.messageBar().pushMessage(
				"georeferencing",
				"Attenzione: PF insufficienti: ne occorrono almeno 2 per la rototraslazione, 3 per la collimazione",
				level=Qgis.Warning,
				duration=3
			)
		else:
			# carica le coordinate origine
			oldCds = []
			for i in listaPf:
				c,x,y,z = pointArchivioCds(self.misurati,i)
				oldCds.append([x,y,z])
#			print "coordinate origine",oldCds
			# ... e quelle di destinazione
			newCds = []
			for i in listaPf:
				c,x,y,z = pointArchivioCds(self.edmPf,i)
				newCds.append([x,y,z])
#			print "coordinate destinazione",newCds
			# ----------- collimazione a 2 PF (rototraslazione) ------------
			nM = len(self.misurati)
			tmp = collimazione2PF(oldCds,newCds,self.misurati+self.ribattuti)
			self.misurati = tmp[0:nM]
			self.ribattuti = tmp[nM:]
			# cancella i layer dei misurati e ribattuti
			#QgsProject.instance().removeMapLayer(self.layLibMisur)
			#QgsProject.instance().removeMapLayer(self.layLibRibat)
			#QgsProject.instance().removeMapLayer(self.layLibCtrn)
			QgsProject.instance().layerTreeRoot().findLayer(self.layLibMisur.id()).setItemVisibilityChecked(False)
			try:
				QgsProject.instance().layerTreeRoot().findLayer(self.layLibCtrn.id())
				layLibCtrn = True                
			except:
				layLibCtrn = False                     
			# ---------------- ricrea i layer ----------------
			# si potrebbero anche aggiornare i layer esistenti
#			self.creaPointLayer('Rilievo_vertici_misurati',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],self.misurati)
#			self.layLibMisur = self.cLayer
#			self.cLayer.setLabelsEnabled(True)
#			print('Layer vertici misurati completato')
#			self.creaPointLayer('Rilievo_vertici_ribattuti',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],self.ribattuti)
#			self.layLibRibat = self.cLayer
#			self.cLayer.setLabelsEnabled(True)
#			print('Layer vertici ribattuti completato')
#			self.creaLineLayer('Rilievo_contorni',self.RilCtrn,self.RilSty,self.misurati)
#			self.layLibCtrn = self.cLayer
#			# ------ attiva la simbologia categorizzata per i contorni -------
#			myRen = catSymbol(
#				self.cLayer.geometryType(),
#				'TRATTO',
#					[
#						['NC','#000000','nero continuo'],
#						['RC','#ff0000','rosso continuo'],
#						['NT','#000000','nero tratteggiato'],
#						['RT','#ff0000','rosso tratteggiato'],
#						['NP','#000000','nero puntinato'],
#						['RP','#ff0000','rosso puntinato']                       
#					]
#			)
#			self.cLayer.setRenderer(myRen)
#			print('Layer contorni libretto completati')
			# ----------- collimazione a 3 PF (rototraslazione) ------------
			# NB: ricordare che in caso di mancata collimazione a 3 PF
			#      i collimati coincidono con i misurati	
			if len(listaPf) >= 3:
				barPFcol,barPFuff, deltaX, deltaY = self.baricentro()
				self.collimati = copy.deepcopy(self.misurati)
				# ricarica le coordinate origine
				oldCds = []
				for i in listaPf:
					c,x,y,z = pointArchivioCds(self.collimati,i)
					oldCds.append([x+deltaX,y+deltaY,z])
#				print "coordinate origine",oldCds
				self.collimati = collimazione3PF(oldCds,newCds,self.collimati)
				# crea layer vertici collimati
				self.creaPointLayer('Rilievo_vertici_collimati',[["indice",QVariant.String],["X",QVariant.Double],["Y",QVariant.Double],["Z",QVariant.Double],["NOTE",QVariant.String],["STAZIONE",QVariant.String],["LIBRETTO",QVariant.Int]],self.collimati)
				self.cLayer.setLabelsEnabled(True)
				self.layLibCollim = self.cLayer
				print('Layer vertici collimati su PF completato')
				if layLibCtrn == True:                
					self.creaLineLayer('Rilievo_contorni_collimati',self.RilCtrn,self.RilSty,self.collimati)
					self.layLibCtrn = self.cLayer
					# ------ attiva la simbologia categorizzata per i contorni -------
					myRen = catSymbol(
						self.cLayer.geometryType(),
						'TRATTO',
						[
							['NC','#000000','nero continuo'],
							['RC','#ff0000','rosso continuo'],
							['NT','#000000','nero tratteggiato'],
							['RT','#ff0000','rosso tratteggiato'],
							['NP','#000000','nero puntinato'],
							['RP','#ff0000','rosso puntinato']                       
						]
					)
					self.cLayer.setRenderer(myRen)               
				self.bCollimList.setEnabled(True)
				self.bErrPf.setEnabled(True)
				self.bDistPf.setEnabled(True)
				self.bBaric.setEnabled(True) 

#	---------- validation functions --------------------

	def errorePF(self):
		"""
			calcola l'errore medio sia per i misurati che per i collimati
			quì è calcolato solo l'errore nel piano XY
		"""
		# controlla i PF misurati nel rilievo
		listaPf = pfLista(self.libretto)
		if len(listaPf) <= 0:				# in assenza di libretto o di PF nel libretto
			for pf in self.edmPf:		# calcola le distanze fra tutti i PF dell'EdM
				listaPf.append(pf[0])
		print('Trovati i seguenti PF:',listaPf,'nel libretto')       
		# --------- misurati --------
		oldCds = []
		for i in listaPf:
			c,x,y,z = pointArchivioCds(self.misurati,i)
			oldCds.append([x,y,0.0])
		#print("Coordinate PF misurate",oldCds)
		newCds = []
		for i in listaPf:
			c,x,y,z = pointArchivioCds(self.edmPf,i)
			newCds.append([x,y,0.0])
		#print("Coordinate PF ufficiali",newCds)
		print('Errore medio dei PF (misurati):')
		print('ErrMin %10.5f ErrMax %10.5f ErrMed %10.5f' % (erroreLsm(oldCds,newCds)))
		# --------- collimati --------
		oldCds = []
		for i in listaPf:
			c,x,y,z = pointArchivioCds(self.collimati,i)
			oldCds.append([x,y,0.0])
		#print("Coordinate PF collimate",oldCds)
		print('Errore medio dei PF (collimati):')
		print('ErrMin %10.5f ErrMax %10.5f ErrMed %10.5f' % (erroreLsm(oldCds,newCds)))
        
	def distanzePF(self): 
		# controlla i PF misurati nel rilievo
		listaPf = pfLista(self.libretto)
		if len(listaPf) <= 0:				# in assenza di libretto o di PF nel libretto
			for pf in self.edmPf:		# calcola le distanze fra tutti i PF dell'EdM
				listaPf.append(pf[0])
		print('Trovati i seguenti PF:',listaPf,'nel libretto')    
		print('Stampa delle distanze tra PF')
		for i0,p0 in enumerate(listaPf):
			j,e0,n0,z = pointArchivioCds(self.edmPf,p0)
			j,x0,y0,z = pointArchivioCds(self.collimati,p0)            
			for i1 in range(i0+1,len(listaPf)):
				p1 = listaPf[i1]
				j,e1,n1,z = pointArchivioCds(self.edmPf,p1)
				j,x1,y1,z = pointArchivioCds(self.collimati,p1)                
				# stampa
				dTAF = math.sqrt((e1-e0)**2+(n1-n0)**2)   
				d = math.sqrt((x1-x0)**2+(y1-y0)**2)  
				# stampa				
				print('distanza tra %s e %s = %12.3f (TAF)' % (p0,p1,dTAF))
				print('distanza tra %s e %s = %12.3f (collimati)' % (p0,p1,d))
				print('---------------------------------')                        

# --------- inquiry functions --------------------

	def librViewTool(self):
		printList(self.libretto)

	def librView(self,point):
		"""
			Fornisce la registrazione del punto nel libretto

			fornisce solo il primo punto
		"""
		# selezione dei punti
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			id,cod,x,y,z,note,stazione,nRiga = pointGetCds(pnts[0])
			QMessageBox.information(
				self.iface.mainWindow(),
				"librView",
				"Point %s:\n Indice:%s\n Riga:%d\n %s" % (id,cod,nRiga,self.libretto[nRiga])
			)
			# pulisce la selezione
			self.cLayer.removeSelection()

	def elencoPfUfficiali(self):
		print('Elenco dei punti fiduciali prelevati da TAF/EdM')
		printList(self.edmPf)

	def distPfUfficiali(self):
		"""
			calcola e stampa le distanze fra PF ufficiali
		"""
		pfList = pfLista(self.libretto)
		if len(pfList) <= 0:				# in assenza di libretto o di PF nel libretto
			for pf in self.edmPf:		# calcola le distanze fra tutti i PF dell'EdM
				pfList.append(pf[0])
		print('trovati i seguenti PF:',pfList,'nel libretto')
		print('Stampa delle distanze ufficiali fra PF')
		for i0,n0 in enumerate(pfList):
			j,x0,y0,z = pointArchivioCds(self.edmPf,n0)
			for i1 in range(i0+1,len(pfList)):
				n1 = pfList[i1]
				j,x1,y1,z = pointArchivioCds(self.edmPf,n1)
				# stampa
				d = math.sqrt((x1-x0)**2+(y1-y0)**2)
				print('%s %12.3f %12.3f' % (n0,x0,y0))
				print('%s %12.3f %12.3f' % (n1,x1,y1))
				print('distanza %12.3f' % (d))
				print('---------------------------------')
				# visualizzazione
				self.rubBnd.addPoint(QgsPointXY(x0,y0))
				self.rubBnd.addPoint(QgsPointXY(x1,y1))
				msg = "%12.3f" % (d)
				x,y = (x0+x1)/2,(y0+y1)/2
				tmp = annotationText(self.canvas,msg,QgsPointXY(x,y))
				self.anTxtList.append(tmp)
		# messaggio finale
		#QMessageBox.information(
		#	self.iface.mainWindow(),
		#	"distanze",
		#	"Done"
		#)
		# elimina gli annotationText
		#for i in self.anTxtList:
		#	self.canvas.scene().removeItem(i)
		#self.canvas.refresh()
		# rimuove la rubber band
		self.rubBnd.reset()

	def leggiDis(self):
		"""
			cerca le distanze fra punti dello stesso foglio dello stesso comune;
			upgrade1: fare stampa ordinata con min e max
			upgrade2: per gestire distanze su fogli/comuni diversi;
		"""
		pfList = pfLista(self.libretto)
		if pfList:
			pfCod,fgCod,comCod = pfList[0].split('/')	# prende uno qualsiasi, il primo
			fgAll = fgCod[3]	# elimina ultimo carattere
			fgCod = fgCod[0:-1]	# elimina ultimo carattere 
			print("cerco i PF del foglio",fgCod,"allegato",fgAll,"del Comune di",comCod)
			# legge file
			fname = QFileDialog.getOpenFileName(self.iface.mainWindow(),'Open file','/home/giuliano','*.dis')
			if fname[0] != "":
				# legge file
				import codecs
				f = codecs.open(fname[0], 'r', encoding='utf-8', errors='ignore')
				print('Stampa dello storico delle distanze fra PF')
				for data in f:
					# pulisce la riga
					data = data.rstrip('\n')
					data = data.rstrip('\r')
					len_comCod =(len(comCod))
					#print('matcha il comune',data[0:len_comCod],comCod)
					#print('matcha il comune',data[18:18+len_comCod])
					#print('matcha il comune',comCod)                    
					if data[0:len_comCod] == comCod:# and data[7] == fgAll:# and data[18:18+len_comCod] == comCod and data[25] == fgAll:
						#print('matchato il comune',data[0:len_comCod])
						if data[7:10] == fgCod and data[7:10] == fgCod:
							print('matcha il foglio',fgCod)
#							pf1 = int(data[15:17])
#							pf2 = int(data[33:35])
#							if pf1 < 10:
#								txt1 = "PF0%1s/%3s0/%5s" % (pf1,fgCod,comCod)
#							else:
#								txt1 = "PF%2s/%3s0/%5s" % (pf1,fgCod,comCod)
#							if pf2 < 10:
#								txt2 = "PF0%1s/%3s0/%5s" % (pf2,fgCod,comCod)
#							else:
#								txt2 = "PF%2s/%3s0/%5s" % (pf2,fgCod,comCod)
#							if txt1 in pfList and txt2 in pfList:	# matchano anche i pf
#								print(txt1,txt2,data[55:62])
			else:
					self.iface.messageBar().pushMessage(
					"importaPfDaTaf",
					"Devi darmi un nome di file valido",
					level=Qgis.Warning,
					duration=4
				)
		else:
			self.iface.messageBar().pushMessage(
				"importaPfDaTaf",
				"Attenzione: il rilievo non sembra contenere punti fiduciali",
				level=Qgis.Warning,
				duration=4
			)

	def elencoPfRilevati(self):
		print('Elenco dei PF rilevati')
		printList(pfLista(self.libretto))

	def distPfMisurate(self):
		"""
			calcola e stampa le distanze misurate
		"""
		pfList = pfLista(self.libretto)
#		print "trovati i PF:",pfList,"nel libretto"
		print('Stampa delle distanze misurate fra PF')
		for i0,n0 in enumerate(pfList):
			j,x0,y0,z = pointArchivioCds(self.misurati,n0)
			for i1 in range(i0+1,len(pfList)):
				n1 = pfList[i1]
				j,x1,y1,z = pointArchivioCds(self.misurati,n1)
				# stampa
				d = math.sqrt((x1-x0)**2+(y1-y0)**2)
				print('%s %12.3f %12.3f' % (n0,x0,y0))
				print('%s %12.3f %12.3f' % (n1,x1,y1))
				print('distanza %12.3f' % (d))
				print('---------------------------------')
				# visualizzazione
				self.rubBnd.addPoint(QgsPointXY(x0,y0))
				self.rubBnd.addPoint(QgsPointXY(x1,y1))
				msg = "%12.3f" % (d)
				x,y = (x0+x1)/2,(y0+y1)/2
				tmp = annotationText(self.canvas,msg,QgsPointXY(x,y))
				self.anTxtList.append(tmp)
		# messaggio finale
		#QMessageBox.information(
		#	self.iface.mainWindow(),
		#	"distanze",
		#	"Done"
		#)
		# elimina gli annotationText
		#for i in self.anTxtList:
		#	self.canvas.scene().removeItem(i)
		#self.canvas.refresh()
		# rimuove la rubber band
		self.rubBnd.reset()

	def elencoEdmVertici(self):
		print('Elenco dei vertici di particella dell EdM')
		printList(self.edmVrts)

	def elencoStazioni(self):
		print('Elenco delle stazioni del rilievo')
		printList(stazioniLista(self.libretto))

	def elencoOssCeler(self):
		print('Elenco delle osservazioni celerimetriche')
		osservazioniCelerimetriche(self.libretto)	

	def stazioniSuVerticeTool(self):
		if self.cLayer.name() in (self.layLibMisur.name(),self.layLibRibat.name()):
			self.archivio = self.misurati + self.ribattuti
		elif self.cLayer.name() == self.layLibCollim .name():
				self.archivio = self.collimati
		else:
			self.iface.messageBar().pushMessage(
				"stazioniSuVertice",
				"The layer is not the right type",
				level=Qgis.Critical,
				duration=3
			)
			return
		self.cLayer.removeSelection()
		# try to disconnect all signals
		if self.isClickToolActivated:
			self.clickTool.canvasClicked.disconnect()
			self.isClickToolActivated = False
		# connect to click signal
		QObject.connect(
			self.clickTool,
			SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
			self.stazioniSuVertice
		)
		self.isClickToolActivated = True
		QMessageBox.information(
			self.iface.mainWindow(),
			"rectAbsPos",
			"Dammi un vertice"
		)

	def stazioniSuVertice(self,point):
		"""
			identifica le stazioni che mirano il vertice
			la ricerca più completa va eseguita sui collimati che sono completi 
			è possibile anche passare la somma dei misurati + ribattuti
		"""
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		lst = []
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			self.cLayer.removeSelection()
			feat = pnts[0]
			cod = feat.attributes()[0]		# prende l'indice
			print('vertice',cod)
			# cerco tutti le stazioni che lo mirano
			for v in self.archivio:
# 				print "controllo",v[0],v[5]
				if v[0] == cod and v[5] != cod:	# scarto la stazione che mira sè stessa
					print('mirato dalla stazione',v[5])
					lst.append(v[5])
		if len(lst):
			drawStar(self,cod,lst,self.archivio)

	def verticiDaStazioneTool(self):
		if self.cLayer.name() in (self.layLibMisur.name(),self.layLibRibat.name()):
			self.archivio = self.misurati + self.ribattuti
#		elif self.cLayer.name() == self.layLibCollim .name():
#				self.archivio = self.collimati
		else:
			self.iface.messageBar().pushMessage(
				"verticiDaStazione",
				"The layer is not the right type",
				level=Qgis.Critical,
				duration=3
			)
			return
		self.cLayer.removeSelection()
		# try to disconnect all signals
		if self.isClickToolActivated:
			self.clickTool.canvasClicked.disconnect()
			self.isClickToolActivated = False
		# connect to click signal
		QObject.connect(
			self.clickTool,
			SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
			self.verticiDaStazione
		)
		self.isClickToolActivated = True
		QMessageBox.information(
			self.iface.mainWindow(),
			"rectAbsPos",
			"Dammi una stazione"
		)

	def verticiDaStazione(self,point):
		"""
			identifica i vertici mirati dalla stazione
			la ricerca viene limitata al layer attivo;
			la ricerca più completa va eseguita sui collimati che sono completi 
		"""
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		lst = []
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			self.cLayer.removeSelection()
			feat = pnts[0]
			cod = feat.attributes()[0]		# prende l'indice
			print('la stazione',cod)
			# cerco tutti i vertici mirati dalla stazione
			for v in self.archivio:
				if v[5] == cod and v[0] != cod:	# scarto la stazione che mira sè stessa
					print('mira',v[0])
					lst.append(v[0])
		if len(lst):
			drawStar(self,cod,lst,self.archivio)

	def elencoMisurati(self):
		print('Elenco delle coordinate misurate dei vertici del rilievo')
		printList(self.misurati)

	def elencoCollimati(self):
		print('Elenco delle coordinate collimate dei vertici del rilievo')
		printList(self.collimati)

	def elencoRibattuti(self):
		print('Elenco delle coordinate ribattute dei vertici del rilievo')
		printList(self.ribattuti)

	def navigazioneStazioni(self):
		lStaz = lettureFraStazioni(self.libretto)
		print(lStaz)
		dlg = navigatorDlg(self.a_giro,self.libretto,lStaz)
		dlg.show()
		dlg.exec_()

	def elencoDistRidotte(self):
 		print('--------- Elenco distanze ridotte ---------')
 		printList(distanzeRidotte(self.libretto,self.a_giro))

	def rectAbsPosTool(self):
		if self.cLayer.geometryType() == QGis.Point:
			self.cLayer.removeSelection()
			# try to disconnect all signals
			if self.isClickToolActivated:
				self.clickTool.canvasClicked.disconnect()
				self.isClickToolActivated = False
			# connect to click signal
			QObject.connect(
				self.clickTool,
				SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
				self.rectAbsPos
			)
			self.isClickToolActivated = True
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"Give me a point"
			)
		else:
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"The layer is not the right type"
			)

	def rectAbsPos(self,point):
		"""
			Fornisce la posizione rettangolare assoluta di un punto
			fornisce solo il primo punto

			non legge più dall'archivio interno perchè questo obbligherebbe
			a usare il plugin solo sull'archivio letto, mentre così può
			essere usato per lavorare su archivi esistenti
		"""
		# selezione dei punti
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			id,cod,x,y,z,note,stazione,nRiga = pointGetCds(pnts[0])
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"Layer: %s\n Point %s:\n Indice:%s\n X=%12.3f\n Y=%12.3f\n Z=%12.3f\n Note:%s" % (self.cLayer.name(),id,cod,x,y,z,note)
			)
		# pulisce la selezione
		self.cLayer.removeSelection()

	def polAbsPosTool(self):
		if self.cLayer.geometryType() == QGis.Point:
			self.cLayer.removeSelection()
			# try to disconnect all signals
			if self.isClickToolActivated:
				self.clickTool.canvasClicked.disconnect()
				self.isClickToolActivated = False
			# connect to click signal
			QObject.connect(
				self.clickTool,
				SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
				self.polAbsPos
			)
			self.isClickToolActivated = True
			QMessageBox.information(
				self.iface.mainWindow(),
				"polAbsPos",
				"Give me a point"
			)
		else:
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"The layer is not the right type"
			)

	def polAbsPos(self,point):
		"""
			Fornisce la posizione polare assoluta di un punto
		"""
		# selezione dei punti
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			id,cod,x,y,z,note,stazione,nRiga = pointGetCds(pnts[0])
			ah,av,d = rect2polar(x,y,z,0.0,0.0,0.0,self.a_giro)
			QMessageBox.information(
				self.iface.mainWindow(),
				"polAbsPos",
				"Layer: %s\n Point %d:\n Indice:%s\n d=%12.3f\n az=%12.4f\n ze=%12.4f\n Notazione angolare:%12.4f\n Note:%s" % (self.cLayer.name(),id,cod,d,ah,av,self.a_giro,note)
			)
		# pulisce la selezione
		self.cLayer.removeSelection()

	def rectRelPosTool(self):
		if self.cLayer.geometryType() == QGis.Point:
			self.cLayer.removeSelection()
			# try to disconnect all signals
			if self.isClickToolActivated:
				self.clickTool.canvasClicked.disconnect()
				self.isClickToolActivated = False
			# connect to click signal
			QObject.connect(
				self.clickTool,
				SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
				self.rectRelPos
			)
			self.isFirst,self.origId,self.x0,self.y0,self.z0 = True,-1,0.0,0.0,0.0
			self.isClickToolActivated = True
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectRelPos",
				"Give me the false origin"
			)
		else:
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"The layer is not the right type"
			)

	def rectRelPos(self,point):
		"""
			Fornisce la posizione rettangolare relativa di un punto

			fornisce solo il primo punto
		"""
		# selezione dei punti
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			id,cod,x,y,z,note,stazione,nRiga = pointGetCds(pnts[0])
			# è il primo punto (falsa origine) ?
			if self.isFirst:
				self.origId = id
				self.origCod = cod
				self.x0,self.y0,self.z0 = x,y,z
				self.isFirst = False
				QMessageBox.information(
					self.iface.mainWindow(),
					"rectRelPos",
					"Give me a vertex"
				)
			else:
				# attiva il RB
				self.rubBnd.addPoint(QgsPoint(self.x0,self.y0))
				self.rubBnd.addPoint(QgsPoint(x,y))
				dx,dy,dz = x-self.x0,y-self.y0,z-self.z0
				QMessageBox.information(
					self.iface.mainWindow(),
					"rectRelPos",
					"Layer: %s\n Point (relative to %s[%s]) %s:\n Indice:%s\n X=%12.3f\n Y=%12.3f\n Z=%12.3f\n Note:%s" % (self.cLayer.name(),self.origId,self.origCod,id,cod,dx,dy,dz,note)
				)
			# pulisce la selezione
			self.cLayer.removeSelection()
			# rimuove la rubber band
			self.rubBnd.reset()

	def polRelPosTool(self):
		if self.cLayer.geometryType() == QGis.Point:
			self.cLayer.removeSelection()
			# try to disconnect all signals
			if self.isClickToolActivated:
				self.clickTool.canvasClicked.disconnect()
				self.isClickToolActivated = False
			# connect to click signal
			QObject.connect(
				self.clickTool,
				SIGNAL("canvasClicked(const QgsPoint &, Qt::MouseButton)"),
				self.polRelPos
			)
			self.isFirst,self.origId,self.x0,self.y0,self.z0 = True,-1,0.0,0.0,0.0
			self.isClickToolActivated = True
			QMessageBox.information(
				self.iface.mainWindow(),
				"polRelPos",
				"Give me the false origin"
			)
		else:
			QMessageBox.information(
				self.iface.mainWindow(),
				"rectAbsPos",
				"The layer is not the right type"
			)

	def polRelPos(self,point):
		"""
			Fornisce la posizione polare relativa di un punto

			fornisce solo il primo punto
		"""
		# selezione dei punti
		searchFeat(self,point)
		nPnts = self.cLayer.selectedFeatureCount() 
		if nPnts >= 1:
			pnts = self.cLayer.selectedFeatures()
			id,cod,x,y,z,note,stazione,nRiga = pointGetCds(pnts[0])
			# è il primo punto (falsa origine) ?
			if self.isFirst:
				self.origId = id
				self.origCod = cod
				self.x0,self.y0,self.z0 = x,y,z
				self.isFirst = False
				QMessageBox.information(
					self.iface.mainWindow(),
					"polRelPos",
					"Give me a vertex"
				)
			else:
				# attiva il RB
				self.rubBnd.addPoint(QgsPoint(self.x0,self.y0))
				self.rubBnd.addPoint(QgsPoint(x,y))
				ah,av,d = rect2polar(x,y,z,self.x0,self.y0,self.z0,self.a_giro)
				QMessageBox.information(
					self.iface.mainWindow(),
					"polRelPos",
					"Layer: %s\n Point (relative to %s[%s]) %s:\n Indice:%s\n d=%12.3f\n az=%12.4f\n ze=%12.4f\n Notazione angolare:%12.4f\n Note:%s" % (self.cLayer.name(),self.origId,self.origCod,id,cod,d,ah,av,self.a_giro,note)
				)
			# pulisce la selezione
			self.cLayer.removeSelection()
			# rimuove la rubber band
			self.rubBnd.reset()

