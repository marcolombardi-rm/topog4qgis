topog4qgis3
==========
VERSIONE ITALIANA
=================
Plugin QGis per la trattazione di libretti Pregeo (atti di aggiornamento catastale)

topo4qgis3 è un plugin per QGIS v.3.xx che consente il trattamento di rilievi celerimetrici con strumentazione ottico-digitale (teodolite, teodolite + distanziometro, stazione totale, ecc.) e, sperimentalmente, con strumentazione gps.
Il plugin si occupa di
- trattare i rilievi eseguiti da diverse stazioni
- verificare i punti ribattuti e
- collegare, mediante opportune trasformazioni affini, le varie stazioni riducendo il rilievo ad un unico spazio vettoriale coerente.

La lettura, nel corso del rilievo, di capisaldi noti (punti fiduciali) e la disponibilità di un estratto di mappa digitale (*.edm) o della tabella dei punti fiduciali (*.taf) consente di georiferire il rilievo nello spazio assoluto ufficiale.
Nel rilievo possono essere trattati anche contorni di unità catastali (particelle, fabbricati, linee dividenti, ecc.).
Nella versione attuale del plugin, i dati sono inseriti mediante il listato delle letture celerimetriche e/o gps, nella forma comunemente chiamata di "libretto" di campagna"; i dati in uscita sono restituiti sotto forma di shp file di tipo puntuale e lineare.
É allo studio una estensione per consentire l'editing (creazione, modifica e cancellazione delle registrazioni) del libretto di campagna.

Il plugin è strettamente legato al formato dei dati forniti dall'Agenzia del Territorio, ora Agenzia delle Entrate, dello stato italiano, pertanto trova completa applicazione in Italia;
la gestione del rilievo è comunque sufficientemente generale da risultare utile, così almeno speriamo, in situazioni anche non nazionali.
Siamo disponibili a, anzi auspichiamo, liberi commenti, critiche e contributi; speriamo che possiate divertirvi usando questo plugin.

Per il corretto funzionamento è necessario utilizzare il giusto Sistema di Riferimento CASSINI-SOLDNER.

Giuliano Curti (original author) e Giuseppe Patti (contributor)

Marco Lombardi (contributor and actual maintainer)

ENGLISH VERSION
===============
QGis Plugin for the italian cadastre update procedure (Pregeo)

topog4qgis3 is a QGis 3.xx plugin allowing the user to manage classical and gps surveys.
Following operations are supported:
- manage surveys from different stations
- verify double points
- connect different station by affine transformation in order to obtain a single survey from many smaller ones.

Georeferencing the survey in an absolute official geographic space is allowed by reading, during surveys, the fiducials and by using the official map from the cadastre database (.edm files) or by reading the fiducials table (.taf).
Surveying of delimiting lines (e.g. between particles, or delimiting a building) is correctly interpreted by the plugin.

In the current version, datas are inserted by listing the book of mutual distances between surveyed points, named "libretto di campagna", output datas are coherently provided in point and linear shapefile.
An extension allowing to edit the "libretto di campagna" is currently under development.
Please note that this plugin is strictly related to datas provided by the italian surveying agency (former Agenzia del Territorio is now Agenzia delle Entrate), and its use is then strictly limited to the italian country, by the way the routines included in it can hopefully be useful also abroad.
We, as developers, will be glad to receive suggestions on how to improve our work. Have fun with it and Happy surveying!!

Giuliano Curti (original author) e Giuseppe Patti (contributor)

Marco Lombardi (contributor and actual maintainer)