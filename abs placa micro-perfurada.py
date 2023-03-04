# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 23:44:01 2021

@author: leofilgueira
"""

from numpy import pi, abs
from matplotlib import pyplot as plt
import numpy as np

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
D = 0.05                            # [m] Comprimento da cavidade do absorvedor
d = 0.0004                          # [m] diâmetro do furo da placa 
l = 0.0004                           # [m] espessura da placa perfurada 
y = d*(((omega*rho0)/(4*viscos))**0.5)
psi = 0.01                          # razão de área perfurada


Kr = ((1+((y**2)/32))**0.5)+(((2**0.5)/32)*y*(d/l))
R = ((32*viscos*l)/(psi*(d**2)))*Kr

Km = (1+((1+((y**2)/2))**(-0.5))+(0.85*(d/l)))
m = ((l*rho0)/psi)*Km


Zpp = 1j*omega*m

Zar =  -1j * rho0 * c0 * (1/np.tan(k0*D))

Zs = R + Zpp + Zar

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