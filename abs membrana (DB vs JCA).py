# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 16:05:32 2021

@author: leofilgueira
"""
# Modelagem de absorvedor de membrana por inidência normal (DB vs JCA)

from numpy import pi, abs
from matplotlib import pyplot as plt
import numpy as np


# Função do modelo de Delany e Bazley
def delany_bazley(f, sig, c0, rho0):
    Zc = (rho0*c0) * ((1 + (9.08*((1e3*(f/sig))**(-0.75)))) - ( 1j * 11.9*((1e3*(f/sig))**(-0.73))))   # Impedância característica do material 
    kc = -1j*(omega/c0) * ((10.3*((1e3*(f/sig))**(-0.59))) + (1j * (1 + (10.8*((1e3*(f/sig))**(-0.7))))))      # Número de onda característico
    
    return Zc, kc

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

# Propriedades do absorvedor
D    = 0.10                        # [m] Comprimento da cavidade do absorvedor
dm   = 1e-4                         # [m] Espessura da membrana
rhom = 1300                         # [kg/m³] Densidade do material da membrana
denssup = dm*rhom                # [kg/m²] Densidade superficial da membrana

# Propriedades do material poroso
sig  = 25000                        # [Rayl] ou [Ns/m^4] Resistividade ao fluxo
d    = 0.07                       # [m] Espessura do material poroso
poros = 0.96                        # [-] Porosidade
tortuos = 1.1                       # [-] Tortuosidade
compvc = 1e-4                       # [um] Comprimento viscoso característico

# Propriedades do Ar
rho0 = 1.21                         # [kg/m³] densidade volumétrica do ar
c0   = 343                          # [m/s] velocidade de propagação do som no ar
P0   = 101325                       # [Pa] pressão atmosférica
prandtl = 0.77                      # [-] Número de Prandtl para o ar
gama = 1.4                          # [-] Razão de calores específicos para o ar
viscos = 1.84e-5                    # [Pa.s] Viscosidade do ar
k0   = omega/c0


Zp1, kp1 = delany_bazley(f, sig, c0, rho0)
Zp2, kp2 = allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0)


Zsp = np.empty((Zp1.size, 2), dtype="complex64")
Zsi = np.empty((Zp1.size, 2), dtype="complex64")
Zs = np.empty((Zp1.size, 2), dtype="complex64")

f60 = 60/(denssup*D)**0.5
f50 = 50/(denssup*D)**0.5
print(f60,f50)

Zm  = 1j * omega * denssup                  # Impedância da membrana
Zsp[:,0] = -1j*(Zp1/(np.tan(kp1*d)))
Zsp[:,1] = -1j*(Zp2/(np.tan(kp2*d)))        # Impedância no topo da camada de material poroso
         
Zsi[:,0] = ((-1j*Zsp[:,0]*rho0*c0*(1/np.tan(k0*(D-d)))) + 
        ((rho0*c0)**2)) / (Zsp[:,0] - (1j*rho0*c0*(1/np.tan(k0*(D-d)))))  # Impedância no topo da camada de ar
Zsi[:,1] = ((-1j*Zsp[:,1]*rho0*c0*(1/np.tan(k0*(D-d)))) + 
        ((rho0*c0)**2)) / (Zsp[:,1] - (1j*rho0*c0*(1/np.tan(k0*(D-d)))))  # Impedância no topo da camada de ar
Zs[:,0] = Zm + Zsi[:,0]                               # Impedância de superfície do absorvedor
Zs[:,1] = Zm + Zsi[:,1]  

                  
reflex = (Zs - (rho0*c0)) / (Zs + (rho0*c0))   # Coeficiente de reflexão
alpha = 1 - (abs(reflex) ** 2)                 # Coeficiente de absorção

# Gráfico (Absorção x Frequência) 

# Legendas
legendas = ['\u03C3 por Delany-Bazley','\u03C3 por Allard-Champoux']

# Gerando gráfico
plt.title("Gráfico de Coef. de Absorção")
    
# a função enumerate retorna um índice correspondente a posição de cada valor
# da lista.
for idc, lgnd in enumerate(legendas):  # para cada índice e legenda na enumeração das legendas
    plt.semilogx(f, alpha[:, idc], label=lgnd)  # faça um gráfico com X logarítmico
    plt.xlabel('Frequência [Hz]')
    plt.ylabel('$\u03B1$ [-]')
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.4)
    plt.show()