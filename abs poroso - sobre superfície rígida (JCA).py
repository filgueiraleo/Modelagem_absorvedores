# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 16:49:51 2021

@author: leofilgueira
"""

# Modelagem de Absorvedor poroso sob superfície rígida (incidência normal)

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
fini = 100                          # [Hz] primeira frequência de análise
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

# Propriedades do material poroso
sig  = 25000                        # [Rayl] ou [Ns/m^4] Resistividade ao fluxo
d    = 0.025                        # [m] Espessura do material poroso
poros = 0.96                        # [-] Porosidade
tortuos = 1.1                       # [-] Tortuosidade
compvc = 1e-4                       # [um] Comprimento viscoso característico


# Impedância característica e numero de onda característico do material poroso
Zc, kc = allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0)

Zs = -1j * (Zc/(np.tan(kc*d)))   # Impedância de Superfície

                   
reflex = (Zs - (rho0*c0)) / (Zs + (rho0*c0))   # Coeficiente de reflexão
alpha = 1 - abs(reflex) ** 2                   # Coeficiente de absorção

# Gerando gráfico
plt.title("Gráfico de Coef. de Absorção")
    
plt.semilogx(f, alpha)  # faça um gráfico com X logarítmico
plt.xlabel('Frequência [Hz]')
plt.ylabel('$\u03B1$ [-]')
plt.grid(True, which="both", ls="--", lw=0.4)
plt.show()