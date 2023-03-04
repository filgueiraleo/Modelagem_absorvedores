# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 16:46:02 2021

@author: leofilgueira
"""

# Modelagem de absorvedor de membrana por inidência normal (Modelo JCA)


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
k0   = omega/c0                     # Numero de onda caracteristico do ar

Zm  = 1j * omega * denssup                  # Impedância da membrana
Zp, kp = allard_champoux(rho0, tortuos, sig, poros, compvc, viscos, omega, prandtl, P0)   # Impedância característica do material poroso

Zsp = -1j*(Zp/(np.tan(kp*d)))        # Impedância no topo da camada de material poroso         
Zsi = ((-1j*Zsp*rho0*c0*(1/np.tan(k0*(D-d)))) + 
        ((rho0*c0)**2)) / (Zsp - (1j*rho0*c0*(1/np.tan(k0*(D-d))))) # Impedância no topo da camada de ar
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