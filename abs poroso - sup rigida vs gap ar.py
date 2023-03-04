# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 18:56:41 2021

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
k0   = omega/c0

# Propriedades do material

sig  = 25000                        # [Rayl] ou [Ns/m^4] Resistividade ao fluxo
d    = 0.025                        # [m] Espessura do material poroso
poros = 0.96                        # [-] Porosidade
tortuos = 1.1                       # [-] Tortuosidade
compvc = 1e-4                       # [um] Comprimento viscoso característico

# Propriedades do Absorvedor
D = 0.035                           # [m] espessura total do absorvedor

Zc, kc = allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0)

# Impedância de superfície 

Zsar =  (-1j * rho0 * c0) / (np.tan(k0 *(D-d)))

Zs[:,0] = ((-1j * Zc1 * Zsar / np.tan(kc1*d)) + Zc1**2) / (Zsar - (1j * (Zc1 / np.tan(kc1*d))))  # Impedância de Superfície DB
                 
reflex = (Zs - (rho0*c0)) / (Zs + (rho0*c0))   # Coeficiente de reflexão
alpha = 1 - abs(reflex) ** 2                   # Coeficiente de absorção

# Gráfico (Absorção x Frequência) 

# Legendas
legendas = ['Sobre superfície rígida (d1 = 30 [mm]']

# Gerando gráfico
plt.title("\u03A6=0.96 [-], \u03C3=25000 [Rayl], \u03B1\u221E=1.1 [-], \u039B=100 [\u03BCm]")
    
# a função enumerate retorna um índice correspondente a posição de cada valor
# da lista.
for idc, lgnd in enumerate(legendas):  # para cada índice e legenda na enumeração das legendas
    plt.semilogx(f, alpha[:, idc], label=lgnd)  # faça um gráfico com X logarítmico
    plt.xlabel('Frequência [Hz]')
    plt.ylabel('$\u03B1$ [-]')
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.4)
    plt.show()


