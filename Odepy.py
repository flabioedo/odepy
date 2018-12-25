# -*- coding: utf-8 -*-
"""
odepy v1: Empleo de métodos numéricos para resolver ecuaciones diferenciales ordinarias
Este módulo esta dirigido a estudiantes de ingeniería y ciencias como una alternativa de libre distribución ante MatLab
El módulo es de libre distribución.
Desarrollado por Flabio César Hernández Salazar
22/12/2018
"""
# Importacion de la funciona a tratar
from funcion import fun
#Inicializacion de la clase
class odepy():
    """
    Esta clase es la parte principal de la libreria de odepy, a traves de la definicion
    del constructor se pueden dar los valores necesarios para realizar los calculos pertinentes
    """
    #Definicion del constructor
    def __init__(self, x_inicial, y_inicial, paso, limitante):
        self.x_inicial = x_inicial
        self.y_inicial = y_inicial
        self.paso = paso
        self.limitante = limitante
    
    #Este es el bloque del metodo de euler
    def euler(self):
        """
        Este es el metodo de euler, deben definirse las condiciones iniciales del
        sistema y el paso. Retorno dos tuplas con los valores a graficar
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        
        while xi < lim:
            #Estos vectores son la información para graficar, y el retorno
            funx.append(xi)
            funy.append(yi)
            # El try/except es utilizado para impedir que el programa se detenga al detectar un error en la funcion
            try:
                lam = float(fun(xi,yi))
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = lam*h
            yj = funcion+yi
            yi = yj
            xi += h
            i += 1
   
        return tuple(funx),tuple(funy)

    # Este es el bloque perteneciente al método de heun, tiene mayor exactitud que el de euler

    def heun(self):
        """
        Método de heun, recibe cuatro valores de inicio: las condiciones iniciales, xi, yi, el paso del método y el valor limite,
        tambien conocido como frontera.
        El método devuelve dos tuplas con los resultados de la iteraciones
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                lam = float(fun(xi,yi))
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = lam*h
            yj = funcion+yi
            yj = yi + ((fun(xi,yi)+fun(xi+h, yj))/(2))*h
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)

    #El siguiente bloque corresponde al metodo del punto medio, que es una variante del metodo de euler con mayor exactitud

    def poligonoMejorado(self):
        """
        El método del poligono mejorado recibe de igual manera las condiciones iniciales, xi, yi, el paso de la funion, h y el valor limite lim
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                lam = float(fun(xi,yi))
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = lam*h/2
            yj = funcion+yi
            yj = yi + fun(xi+h/2, yj)*h
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)

    # El siguiente bloque corresponde al método de ralston con a2 = 2/3

    def ralston(self):
        """
        Este método corresponde a los métodos de runge-kutta de segundo orden, éste
        minimiza el error por truncamiento de entre los diferentes valores que se le pueden
        dar a a1 y a2.
        Para desarrollar el método se necesitan las condiciones iniciales, xi y yi, el paso
        de la función, h, y el limite de la frontera, lim
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                k1 = fun(xi,yi)
                k2 = fun(xi + (3.0/4)*h, yi + (3.0/4)*k1*h)
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = ((1.0/3)*k1 + (2.0/3)*k2)*h
            yj = funcion + yi
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)

    # El siguiente bloque corresponde al método de runge-kutta de tercer orden

    def RK3(self):
        """
        Para este método se necesita especificar las condiciones iniciales del sistema
        xi, yi, además del paso del método, h, y el valor de frontera, lim.
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                k1 = fun(xi,yi)
                k2 = fun(xi + (1.0/2)*h, yi + (1.0/2)*k1*h)
                k3 = fun(xi + h, yi - k1*h + 2*k2*h)
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = ((k1 + 4*k2 +k3)*(1.0*h))/6
            yj = funcion + yi
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)

    # El siguiente bloque corresponde al método de runge-kutta de cuarto orden

    def RK4(self):
        """
        Para este método se deben especificar las condiciones iniciales del sistema, xi, yi,
        además del paso del método, h, y el valor de frontera, lim
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                k1 = fun(xi,yi)
                k2 = fun(xi + (1.0/2)*h, yi + (1.0/2)*k1*h)
                k3 = fun(xi + (1.0/2)*h, yi + (1.0/2)*k2*h)
                k4 = fun(xi + h, yi + k3*h)
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = ((k1 + 2*k2 +2*k3 + k4)*(1.0*h))/6
            yj = funcion + yi
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)

    # El siguiente bloque corresponde al método de runge-kutta de quinto orden, o método de Butcher

    def butcher(self):
        """
        El método de butcher tiene gran exactitud, sin embargo esto conlleva mayor trabajo computacional
        por lo que en muchas ocasiones se prefiere usar el método de runge-kutta de cuarto orden
        Para este método se deben especificar las condiciones iniciales del sistema, xi, yi,
        además del paso del método, h, y el valor de frontera, lim
        """
        xi = self.x_inicial
        yi = self.y_inicial
        h = self.paso
        lim = self.limitante
        
        i = 0
        funx = []
        funy = []
        while xi < lim:
            funx.append(xi)
            funy.append(yi)
            try:
                k1 = fun(xi,yi)
                k2 = fun(xi + (1.0/4)*h, yi + (1.0/4)*k1*h)
                k3 = fun(xi + (1.0/4)*h, yi + (1.0/8)*k1*h + (1.0/8)*k2*h)
                k4 = fun(xi + (1.0/2)*h, yi - (1.0/2)*k2*h + k3*h)
                k5 = fun(xi + (3.0/4)*h, yi + (3.0/16)*k1*h + (9.0/16)*k4*h)
                k6 = fun(xi + h, yi - (3.0/7)*k1*h + (2.0/7)*k2*h + (12.0/7)*k3*h - (12.0/7)*k4*h + (8.0/7)*k5*h)
            except:
                print "El método no se puede aplicar a la funcion"
            funcion = ((7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)*(1.0*h))/90
            yj = funcion + yi
            yi = yj
            xi += h
            i += 1
        return tuple(funx),tuple(funy)
    