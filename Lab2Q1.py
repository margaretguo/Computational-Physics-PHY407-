import numpy as np
import matplotlib.pyplot as plt

##1a## using the central difference method, calculate the derivative for f(x) = e^2x at x = 0 for h values ranging from 10-16 to 1 
def function(x): #outputs the core dawson function
    return np.exp(2*x)

def centraldifferences(x):
    derivative_values = []
    for i in range(-16, 1):
        h = 10**i #set the h value for this iteration
        derivative = (function(x + h/2) - function(x - h/2)) / h #calculate the derivative using the central difference formula
        print(str(h) + ": " + str(derivative)) #print the derivative value at each h value
        derivative_values.append(derivative)
    return derivative_values

centraldifferences(0)


    
##1b## calculate the difference between each result and exactly two and plot the error as a function of stepsize/h
def absoluteerror():
    true_derivative = 2.0
    derivative_values = centraldifferences(0)
    errors = []
    
    for i in derivative_values: 
        error = np.abs(i - true_derivative) #calculate and save each absolute error value
        errors.append(error)
    return errors

##generate plot##
##generate xvalues array
xvalues = []
for i in range(-16, 1):
    h = 10**i
    xvalues.append(h)

yvalues = absoluteerror()
plt.loglog(xvalues, yvalues, "r^") #use a log log plot as per instructions 
plt.title("Step size (h) vs absolute error")
plt.xlabel("Step size")
plt.ylabel("Absolute error")
plt.show()




##1c## calculate the frist 10 orders of derivatives f'(x) - f10(x) with x = 0 using central derivative approximation
#using the h value with the lowest error. 

#use provided recursive derivative calulation function. m = order of the derivative 
def delta(f, x, m, h):
    if m > 1:
        return (delta(f, x + h/2, m - 1, h) - delta(f, x - h/2, m-1, h))/(h)
    else:
        return (f(x + h/2) - f(x - h/2))/(h)

h = 10**-5 #use the h value with the lowest error
nonoph = 10**-1 ##this h produces the closest values to the actual derivative values, the size needed to be increased so that it wouldn't be lower than the machine error for higher derivative values

#calculate and print the derivatives
for i in range(1, 11):
    optimalhderivative = delta(function, 0, i, h)
    nonoptimalhderivative = delta(function, 0, i, nonoph)
    check = 2**i #actual derivative values 
    print(str(i) + "  " + str(check) + " " + str(optimalhderivative) + " " + str(nonoptimalhderivative))
    
