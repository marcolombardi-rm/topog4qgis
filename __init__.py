# -*- coding: utf-8 -*-

def name():
	return "topog4qgis3"

def description():
	return "Experimental topographic tools for qGis"

def version():
	return "0.00"

def icon():
    return ':/plugins/topog4qgis3/icon.png'

def qgisMinimumVersion():
	return "3.0"

def authorName():
	return "giuliano curti"

def classFactory(iface):
    from .topog4qgis3 import topog4qgis3
    return topog4qgis3(iface)