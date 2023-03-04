# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 23:31:10 2021

@author: leofilgueira
"""

from numpy import pi, abs
from matplotlib import pyplot as plt
import numpy as np

# Função do modelo Johnson-Champoux-Allard
def allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0):
    A = (sig*poros)/(1j*tortuos*rho0*omega)                                     # Termo auxiliar A
    B = (4j*(tortuos**2)*viscos*rho0*omega)/((sig**2)*(compvc**2)*(poros**2))   # Termo auxiliar B
    rhoc = (rho0*tortuos)*(1 + (A*((1 + B)**0.5)))                              # Densidade característica
    
    C = (8*viscos)/(1j*(compvc**2)*rho0*prandtl*omega)                       # Termo auxiliar C
    D = (1+((1j*rho0*(compvc**2)*prandtl*omega)/(16*viscos)))**0.5           # Termo auxiliar D
    K = (gama*P0)/(gama-((gama-1)*(((1+(C*D)))**-1)))                        # Módulo de Compressibiliade
    
    Zc = (K*rhoc)**0.5
    kc = omega*((rhoc/K)**0.5)

    return Zc, kc

# Espaço de frequências
fini = 20                          # [Hz] primeira frequência de análise
ffim = 10000                        # [Hz] última frequência de análise
f = np.linspace(fini, ffim, 10000)  # [Hz] espaço linear de frequências
omega = 2 * pi * f                  # [rad/s] espaço linear de frequências angulares

# Propriedades do Ar
rho0 = 1.21                         # [kg/m³] densidade volumétrica do ar
c0   = 343                          # [m/s] velocidade de propagação do som no ar
P0   = 101325                       # [Pa] pressão atmosférica
prandtl = 0.77                      # [-] Número de Prandtl para o ar
gama = 1.4                          # [-] Razão de calores específicos para o ar
viscos = 1.84e-5                    # [Pa.s] Viscosidade do ar
k0   = omega/c0                     # número de onda no ar

# Propriedades do absorvedor
D = 0.12                             # [m] Comprimento da cavidade do absorvedor
w = 0.005                          # [m] largura da ranhura
a = 0 
b = 0.014                           # [m] espaçamento entre as ranhuras 
#psi = (w*a)/(b**2)                # razão de área perfurada
psi = 0.1
l = 0.001                           # [m] espessura da placa perfurada 
lc = l + 2*w*(-1/np.pi)*(np.log(np.sin(0.5*np.pi*psi)))         # [m] comprimento corrigido (impedância de radiação)

#denssup = (rho0*(b**2)*lc)/(np.pi*(a**2))  

denssup = rho0*(1/psi)*(lc+(((8*viscos/omega)*(1+(l/2*a)))**0.5)) # [kg/m²] densidade superficial do gás contido na perfuração

# Propriedades do material poroso
sig  = 20000                        # [Rayl] ou [Ns/m^4] Resistividade ao fluxo
d    = 0.06                      # [m] Espessura do material poroso
poros = 0.96                        # [-] Porosidade
tortuos = 1.1                       # [-] Tortuosidade
compvc = 1e-4                       # [um] Comprimento viscoso característico

fres =  (c0/(2*np.pi)) * ((psi/(D*lc))**0.5)
print(fres)

Zc, kc = allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0)

Zpp = 1j*omega*denssup

Zsar =  (-1j * rho0 * c0) / (np.tan(k0 *(D-d)))

Zsi = ((-1j * Zc * Zsar / np.tan(kc*d)) + Zc**2) / (Zsar - ((1j * Zc) / np.tan(kc*d)))

Zs = Zpp + Zsi

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