import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw

##3a##
##write a function that calculates the hermite polynomial for a given n and x value. try to use recursion
def H(n, x): #n = order of hermite polynomial, x = function variable. returns a single value (the computed polynomial)
    if n == 1:
        return 2*x
    elif n == 0:
        return 1
    else: 
        return 2*x*H(n-1, x) - 2*(n-1)*H(n-2, x) ##use n-1 instead of n since equation 4 gives us Hn+1

##3b##
##plot the harmonic oscillator wavefunctions for n=0-3 for x in range of -4 to 4 all on one graph

def psi(x, n): ##takes in a specific x and n value, returns a single value
    fraction = 1/np.sqrt(2**n*np.math.factorial(n)*np.sqrt(np.pi))
    exponent = np.exp(-x**2/2)
    hermite = H(n, x)
    return fraction * exponent * hermite


xvalues = [x / 10.0 for x in range(-40, 40, 1)] #create list of x values in range to use to calculate wavefunction values, same for each n value calculation
zeron = [] #create lists for the wavefunction values 
onen = []
twon = []
threen = []
#n = 0
n = 0
for x in xvalues:
    value = psi(x, n)
    zeron.append(value)

#n = 1
n = 1
for x in xvalues:
    value = psi(x, n)
    onen.append(value)

#n = 2
n = 2
for x in xvalues:
    value = psi(x, n)
    twon.append(value)

#n = 3
n = 3
for x in xvalues:
    value = psi(x, n)
    threen.append(value)

##create plot
plt.plot(xvalues, zeron, "r--",  xvalues, onen, "b--", xvalues, twon, "g--", xvalues, threen, "m--")
plt.title("Harmonic Oscillator Wavefunction Values vs x for n = 0-3")
plt.xlabel("X values")
plt.ylabel("Wavefunction values")
plt.legend(["n = 0", "n = 1", "n = 2", "n = 3"])
plt.show()


##3c##
##calculate the potential energy of the oscillator using equation 5 divided by 2 using Gaussian Quadrature for N = 100
##for hermite polynomials 0-10

##use equation 5.75 from textbook (trignometric subsitution)

def potential(z, n): #calculates the potential energy core function using equation 5 and the integration subsitution
    return (np.tan(z)**2 * np.abs(psi(np.tan(z), n))**2) / np.cos(z)**2
#calls psi() since the wavefunction value is needed. 

##gaussian quadrature function, modified from question 2
def gaussian(n):
    N = 100 #number of slices/sample points 
    a = -np.pi/2 #start point 
    b = np.pi/2 #end point 

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration
    s = 0.0
    for k in range(N):
        s += wp[k]*potential(xp[k], n)

    s = s / 2 ##since it is twice the potential energy
    print(str(n) + " " + str(s))
    return s

##print energy values
for i in range(0, 11):
    gaussian(i)












