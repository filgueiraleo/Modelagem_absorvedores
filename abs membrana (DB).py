# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 16:05:32 2021

@author: leofilgueira
"""

# Modelagem de absorvedor de membrana por inidência normal (Modelo Dellany-Bazley)

from numpy import pi, abs
from matplotlib import pyplot as plt
import numpy as np

# Função do modelo de Delany e Bazley
def delany_bazley(f, sig, c0, rho0):
    Zc = (rho0*c0) * ((1 + (9.08*((1e3*(f/sig))**(-0.75)))) - ( 1j * 11.9*((1e3*(f/sig))**(-0.73))))   # Impedância característica do material 
    kc =(omega/c0) * ((10.3*((1e3*(f/sig))**(-0.59))) + (1j * (1 + (10.8*((1e3*(f/sig))**(-0.7))))))      # Número de onda característico
    
    return Zc, kc

# Espaço de frequências
fini = 20                          # [Hz] primeira frequência de análise
ffim = 10000                        # [Hz] última frequência de análise
f = np.linspace(fini, ffim, 10000)  # [Hz] espaço linear de frequências
omega = 2 * pi * f                  # [rad/s] espaço linear de frequências angulares

# Propriedades do Ar
rho0 = 1.21                         # [kg/m³] densidade volumétrica do ar
c0   = 343                          # [m/s] velocidade de propagação do som no ar
k0   = omega/c0

# Propriedades do absorvedor
sig  = 25000                        # [Rayl] ou [Ns/m^4] Resistividade ao fluxo do material poroso
D    = 0.10                       # [m] Comprimento da cavidade do absorvedor
d    = 0.07                     # [m] Espessura do material poroso
dm   = 1e-4                         # [m] Espessura da membrana
rhom = 1300                        # [kg/m³] Densidade do material da membrana
denssup = dm*rhom                 # [kg/m²] Densidade superficial da membrana

f60 = 60/(denssup*D)**0.5
f50 = 50/(denssup*D)**0.5
print(f60,f50)

Zm  = 1j * omega * denssup                  # Impedância da membrana
Zp, kp = delany_bazley(f, sig, c0, rho0)    # Impedância característica do material poroso
Zsp = Zp * ((np.cosh(kp*d))/(np.sinh(kp*d))) # Impedância no topo da camada de material poroso         
Zsi  = ((-1j*Zsp*rho0*c0*(1/np.tanh(k0*(D-d)))) +
        ((rho0*c0)**2)) / (Zsp - (1j*rho0*c0*(1/np.tanh(k0*(D-d)))))  # Impedância no topo da camada de ar
Zs = Zm + Zsi                               # Impedância de superfície do absorvedor

# Coeficiente de absorção alpha                   
reflex = (Zs - (rho0*c0)) / (Zs + (rho0*c0))   # Coeficiente de reflexão
alpha = 1 - (abs(reflex) ** 2)

# Gerando gráfico
plt.title("Gráfico de Coef. de Absorção")
    
plt.semilogx(f, alpha)  # faça um gráfico com X logarítmico
plt.xlabel('Frequência [Hz]')
plt.ylabel('$\u03B1$ [-]')
plt.grid(True, which="both", ls="--", lw=0.4)
plt.show()