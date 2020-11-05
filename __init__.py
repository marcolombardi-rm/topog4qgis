# -*- coding: utf-8 -*-

def name():
	return "topog4qgis"

def description():
	return "Experimental topographic tools for qGis"

def version():
	return "0.3.2"

def icon():
    return ':/plugins/topog4qgis/icon.png'

def qgisMinimumVersion():
	return "3.0"

def authorName():
	return "giuliano curti"

def classFactory(iface):
    from .topog4qgis import topog4qgis
    return topog4qgis(iface)