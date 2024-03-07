import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw
from scipy.special import dawsn as scipydawson

##2a## Calculate the dawson function at x = 4 using trapezoidal rule, simpson's rule and gaussian quadrature for a range of slices between 8 and 2048

def dawsonfunction(x, t):
    '''This function allows us to calculate the core dawson function value for a given input, it returns a single value'''
    dawson = np.exp(t**2)
    return dawson

def simpson(slicenumber):
    '''This calculates the simpson approximation for a number of slices slices. it calls the dawsonfunction() function and outputs a single value'''
    #set the starting points of the dawson integral
    a = 0.0 
    b = 4

    #calculate our slice width size for 1000 slices
    h = (b-a)/slicenumber

    #calculate f(a) + f(b) according to the formula
    integralvalue = dawsonfunction(4, a) + dawsonfunction(4, b)

    #start a for loop to implmement the rest of the slices
    N = list(range(1, slicenumber))
    for k in N:
        t = a + k*h #k will increase by 1 each time 
        if k % 2 == 1:  #if k is odd
            integralvalue = integralvalue + 4.0 * dawsonfunction(4, t)
        if k % 2 == 0: #if k is even
            integralvalue = integralvalue + 2.0 * dawsonfunction(4, t)
    #multiply final term by h/3 to complete expression
    integralvalue = integralvalue * (h/3.0)

    #multiply by exp(-x**2) to get final dawson approximation
    integralvalue = integralvalue * np.exp(-4**2)

    return integralvalue


def trapezoidal(slicenumber):
    '''This calculates the trapezoidal approximation for a number of slices. it calls the dawsonfunction() function and outputs a single value'''
    #set the starting points of the dawson integral
    a = 0.0 
    b = 4.0

    #calculate our slice width size for 1000 slices
    h = (b-a)/slicenumber

    #calculate 0.5f(a) + 0.5f(b) according to the formula
    integralvalue = 0.5*dawsonfunction(4, a) + 0.5*dawsonfunction(4, b)

    #start a for loop to implmement the rest of the slices
    N = list(range(1, slicenumber))
    for k in N:
        t = a + k*h #k will increase by 1 each time 
        integralvalue = integralvalue + dawsonfunction(4, t)
    
    #multiply final term by h to complete expression
    integralvalue = integralvalue * h

    #multiply by exp(-x**2) to get final dawson approximation
    integralvalue = integralvalue * np.exp(-4**2)

    return(integralvalue)

    ##gaussian quadrature method using gaussxw.py from textbook 
def gaussian(slicenumber):
    N = slicenumber #number of slices/sample points 
    a = 0.0 #start point 
    b = 4.0 #end point 

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration
    s = 0.0
    for k in range(N):
        s += wp[k]*dawsonfunction(4, xp[k])

    s = s * np.exp(-(4)**2)
    return s

##calculate integral for slice number = 2**3-2**11 for each method & plot it with the true value from scipy

truevalue = scipydawson(4)
slicepowers = [3, 4, 5, 6, 7, 8, 9, 10, 11]
xvalues = []
trapyvalues = []
simpsonyvalues = []
gaussyvalues = []

##trapezoidal##
for i in slicepowers:
    slicenumber = 2**i
    xvalues.append(slicenumber) #build our xvalues list for plotting (same for each)
    integral = trapezoidal(slicenumber)
    trapyvalues.append(integral)

plt.plot(xvalues, trapyvalues, "r^")
plt.axhline(y = truevalue, color = "b")
plt.title("Trapezoidal Approximation vs Slice Numbers for Dawson Integral")
plt.xlabel("Slice Number")
plt.ylabel("Integral Value")
plt.legend(["Approximation Integral Values", "True Value"])
plt.show()

##simpson##
for i in slicepowers:
    slicenumber = 2**i
    integral = simpson(slicenumber)
    simpsonyvalues.append(integral)

plt.plot(xvalues, simpsonyvalues, "m^")
plt.axhline(y = truevalue, color = "b")
plt.title("Simpson Approximation vs Slice Numbers for Dawson Integral")
plt.xlabel("Slice Number")
plt.ylabel("Integral Value")
plt.legend(["Approximation Integral Values", "True Value"])
plt.show()


##gaussian quadrature##
for i in slicepowers:
    slicenumber = 2**i
    integral = gaussian(slicenumber)
    gaussyvalues.append(integral)

plt.plot(xvalues, gaussyvalues, "g^")
plt.axhline(y = truevalue, color = "b")
plt.title("Gaussian Quadrature Approximation vs Slice Numbers for Dawson Integral")
plt.xlabel("Slice Number")
plt.ylabel("Integral Value")
plt.legend(["Approximation Integral Values", "True Value"])
plt.show()


##2b##
#calculate the error of the gaussian quadrature method using equation 2 in the manual. Plot the error as a function of log2(N) 
##describe the error pattern. 
errors = []
for i in slicepowers:
    slicenumber = 2**i
    integral = gaussian(slicenumber) #calculate our initial integral
    doubleslicenumber = slicenumber * 2 #calculate 2n
    integral2 = gaussian(doubleslicenumber)
    error = integral2 - integral
    errors.append(np.abs(error)) #take magnitude of each error

print(errors)
plt.plot(slicepowers, errors, "g^") #plot the log2 of each slice number
plt.title("Gaussian Quadrature Approximation errors vs. log2(slice number)")
plt.xlabel("log2(Slice Number)")
plt.ylabel("Error Value")
plt.show()

##2c##
##plot the relative error (calculated gaussian quadrature integral - truevalue) vs. N 
relativeerrors = []
for integral in gaussyvalues:
    error = np.abs(integral - truevalue) #calculate relative error based on formula
    relativeerrors.append(error)

print(relativeerrors)
plt.plot(slicepowers, relativeerrors, "c^") #plot the log2 of each slice number
plt.title("Gaussian Quadrature Approximation relative errors vs. log2(slice number)")
plt.xlabel("log2(Slice Number)")
plt.ylabel("Error Value")
plt.show()

