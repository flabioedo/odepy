# -*- coding: utf-8 -*-
"""
Prueba del funcionamiento de la libreria odepy
"""
# importacion del modulo a utilizar
from matplotlib import pyplot
from Odepy import odepy

ejem1 = odepy(0,0,0.1,2)
a , b = ejem1.RK4()
pyplot.plot(a,b)
pyplot.show()
