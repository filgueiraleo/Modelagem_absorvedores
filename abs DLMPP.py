# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 00:56:01 2021

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

# Propriedades do absorvedor MPP 1
D1 = 0.025                            # [m] Comprimento da cavidade do absorvedor
d1 = 0.0004                          # [m] diâmetro do furo da placa 
l1 = 0.0004                           # [m] espessura da placa perfurada 
y1 = d1*(((omega*rho0)/4*viscos)**0.5)
psi1 = 0.01                          # razão de área perfurada

Kr1 = ((1+((y1**2)/32))**0.5)+(((2**0.5)/32)*y1*(d1/l1))
R1 = ((32*viscos*l1)/(psi1*(d1**2)))*Kr1

Km1 = 1+((1+((y1**2)/2))**-0.5)+(0.85*(d1/l1))
m1 = (l1*rho0/psi1)*Km1

# Propriedades do absorvedor MPP 2
D2 = 0.025                           # [m] Comprimento da cavidade do absorvedor
d2 = 0.0004                          # [m] diâmetro do furo da placa 
l2 = 0.0004                           # [m] espessura da placa perfurada 
y2 = d2*(((omega*rho0)/4*viscos)**0.5)
psi2 = 0.01                          # razão de área perfurada

Kr2 = ((1+((y2**2)/32))**0.5)+(((2**0.5)/32)*y2*(d2/l2))
R2 = ((32*viscos*l2)/(psi2*(d2**2)))*Kr2

Km2 = 1+((1+((y2**2)/2))**-0.5)+(0.85*(d2/l2))
m2 = (l2*rho0/psi2)*Km2


Zpp1 = 1j*omega*m1
Zpp2 = 1j*omega*m2

Zar1 =  -1j * rho0 * c0 * (1/np.tan(k0*D1))
Zar2 =  -1j * rho0 * c0 * (1/np.tan(k0*D2))

Zs = R1 + Zpp1 + (((1/Zar1)+(1/(R2+Zpp2+Zar2)))**-1)

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