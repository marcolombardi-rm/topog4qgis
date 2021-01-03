topog4qgis v0.3.3
==========
VERSIONE ITALIANA
=================
topog4qgis e' un plugin utile all'elaborazione di libretti PreGeo (.dat e .pdf) e alla trattazione di liste di punti (.csv) su QGIS mediante opportuna rototraslazione ai minimi quadrati.

L'elaborazione di libretti PreGeo (.dat e .pdf)  e' strettamente legata al formato dei dati forniti dall'Agenzia del Territorio, ora Agenzia delle Entrate, dello stato italiano, pertanto trova completa applicazione in Italia.

Il plugin puo' elaborare rilievi celerimetrici, GNSS, misti TPS-GNSS; verificando i punti ribattuti e collegando, mediante opportune trasformazioni affini, le varie stazioni riducendo il rilievo ad un unico spazio vettoriale coerente.  Vengono elaborati dal plugin punti definiti con allineamenti e squadri (riga 4 e 5).

Se presenti,  sono estratti ed elaborati, anche contorni di unita' catastali (particelle, fabbricati, linee dividenti, ecc.) definite con riga 7.

La trattazione di liste di punti (.csv), dalla versione 0.3.2, e' invece slegata dai formati PreGeo consentendo quindi all'utente di lavorare su un rilievo gia' elaborato da altri software.

La lettura, nel corso del rilievo, di capisaldi noti (punti fiduciali) e la disponibilita' di un estratto di mappa digitale (.edm) o della tabella dei punti fiduciali (.taf o .csv) consente di georeferire mediante rototraslazione ai minimi quadrati il rilievo nello spazio assoluto ufficiale.

Dalla versione 0.3.1 e' anche possibile l'esportazione in formato .csv del rilievo elaborato e del rilievo rototraslato ai minimi quadrati.

Siamo disponibili, anzi auspichiamo, liberi commenti, critiche e contributi; speriamo che possiate divertirvi usando questo plugin. 

All'indirizzo web https://github.com/marcolombardi-rm/topog4qgis/wiki troverete una breve guida per l'utilizzo del plugin.

Per il corretto funzionamento è necessario utilizzare il giusto Sistema di Riferimento CASSINI-SOLDNER o GAUSS-BOAGA.

Giuliano Curti (original author) e Giuseppe Patti (contributor)

Marco Lombardi (contributor and actual maintainer)

ENGLISH VERSION
===============
QGis Plugin for the italian cadastre update procedure (PreGeo)

topog4qgis is a QGis 3 plugin allowing the user to manage classical and gps surveys.
Following operations are supported:
- manage surveys from different stations
- verify double points
- connect different station by affine transformation in order to obtain a single survey from many smaller ones.

Georeferencing the survey in an absolute official geographic space is allowed by reading, during surveys, the fiducials and by using the official map from the cadastre database (.edm files) or by reading the fiducials table (.taf).
Surveying of delimiting lines (e.g. between particles, or delimiting a building) is correctly interpreted by the plugin.

In the current version, datas are inserted by listing the book of mutual distances between surveyed points, named "libretto di campagna", output datas are coherently provided in point and linear shapefile.
Please note that this plugin is strictly related to datas provided by the italian surveying agency (former Agenzia del Territorio is now Agenzia delle Entrate), and its use is then strictly limited to the italian country, by the way the routines included in it can hopefully be useful also abroad.
We, as developers, will be glad to receive suggestions on how to improve our work. Have fun with it and Happy surveying!!

Giuliano Curti (original author) e Giuseppe Patti (contributor)

Marco Lombardi (contributor and actual maintainer)
