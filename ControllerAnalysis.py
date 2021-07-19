import numpy as np              #import numpy library
import scipy as sp              #import scipy library
from scipy import stats
from scipy.stats import norm   
from scipy.optimize import curve_fit
from scipy.signal import medfilt

#This script will be used to perform stability analysis on the 2D inverted pendulum state space model. 

#Constants
#general
g = 9.81 #gravity, m/s^2

#Pendulum Body and rotor stator
mb  = 1.0 #mass of the body, kg
lb  = 1.0 #distance from the pivot to the body CG, m
Ibo = 1.0 #mass moment of inertia of the pendulum body about the pivot point o, kg*m^2
Cb  = 1.0 #coefficient of friction of of the pendulum body pivot,  (kg*m^2)/s

#Rotor
mm = 1.0 #mass of the rotor, kg
lm = 1.0 #distance from the pivot to the rotor CG, m
Im = 1.0 #mass moment of inertia of the rotor about its fixed axis, kg*m^2
Cm = 1.0 #coefficent of friction of the rotor bearings, (kg*m^2)/s
Km = 1.0 #bushless DC motor constant, 

#State Space model
#xdot(t) = Ax(t)+Bu(t)
#A
a11 = 0.0
a12 = 1.0 
a13 = 0.0

a21 = (g*(mb*lb+mm*lm))/(Ibo+Im+mm*lm**2)
a22 = -Cb/(Ibo+Im+mm*lm**2)
a23 = Cm/(Ibo+Im+mm*lm**2)

a31 = -(g*(mb*lb+mm*lm))/(Ibo+Im+mm*lm**2)
a32 = Cb/(Ibo+Im+mm*lm**2)
a33 = (Cm*(Ibo+mm*lm**2))/(Im*(Ibo+Im+mm*lm**2))

A = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])

#B

B = np.array([0, -Km/(Ibo+Im+mm*lm**2), (Km*(Ibo+mm*lm**2))/(Im*(Ibo+Im+mm*lm**2))])

print('A = ' + str(A))
print('')
print('B = ' + str(B))
print('')

#Stability Analysis of the systems

eigenVals = np.linalg.eig(A)
print(eigenVals)
