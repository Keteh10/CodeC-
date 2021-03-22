from matplotlib.pyplot import *
import numpy as np
import math


G = 6,67*10**(-11)
c_s = 1
M = 2*10**30
rho_0 = 1*10**7


def a(t):
    return(-(32/15)*math.pi*(1/4)*((G*M/c_s**3)**(5/3))*((G*M)**(-1/2))*G*rho_0*((G*M)**(5/6)))

def affichage(fonction, t):
    Separation = []
    T = [3600*j for j in range(t)]
    for i in range(t):
        Separation.append(fonction(t))
    plot(Separation,T,'r')
    xlabel('t')
    ylabel('a(t)')
    title('Evolution de la separation en fonction du temps')
    show()
    return()

if __name__=="__main__":
    t = 10000
    affichage(a, t)
    
