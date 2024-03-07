import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sp
from gaussxw import gaussxw

#Question 2 Problem Description##
##calculate the efficiency of the lightbulb by first the total energy using the equations provided for the wavelengths provided and dividing by the total energy emitted
##Then graph the efficiency as a function of the temperature it's operating at. 
##Then calculate the temperature at which the bulb is most efficient using a search method from ch 6.3. Finally, calculate that same 
##value but for a different range of wavelengths. 

##2a##

##2a solution description##
##use gaussian quadrature since it is the best integration method. can remove A from equation because it cancels out during integration. 
#create a function for the core power function ilambda and integrate that from the wavelengths provided. 
# Integrate again from 0-infinity and take the ratio of that as the efficiency. 
#of the efficiency equation to calculate the efficiency ratio n another funtion. 

#calculates the function values using the core power function in equation 8
def ilambda(T, l): 
    h = sp.Planck #plancks constant
    c = sp.speed_of_light
    kb = sp.Boltzmann

    neumerator = 2 * np.pi * h * c**2
    denominator = l**5 * (np.exp((h*c)/(l* kb * T)) - 1)

    return neumerator/denominator

#performs the total energy integral in equation 9, takes in the temp, wavelengths and returns a single value
def total_energy_gaussian(T, a, b): 
    N = 5000 #number of slices/sample points 
    a = a #start point in nm
    b = b #end point in nm

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration
    s = 0.0
    for k in range(N):
        s += wp[k]*ilambda(T, xp[k])
    return s

#ilambda function with the infinite integral subsitution from page 179 in textbook implemented. 
def subsitution_ilambda(T, z):
    h = sp.Planck #plancks constant
    c = sp.speed_of_light
    kb = sp.Boltzmann

    neumerator = 2 * np.pi * h * c**2
    denominator = (z/(1-z))**5 * (np.exp((h*c)/((z/(1-z))* kb * T)) - 1)

    return (neumerator/denominator) * (1/(1-z)^2)

    
#calculate the efficiency using equation 10, takes in temperature, wavelengths and returns a single efficiency value
def efficiency(T, lambda1, lambda2): 
    h = sp.Planck #plancks constant
    c = sp.speed_of_light #speed of light
    kb = sp.Boltzmann #boltzmann constant
    neumerator = total_energy_gaussian(T, lambda1, lambda2)
    #denominator = (2 * np.pi**5 * kb**4 * T**4)/(15 * c**2 * h**3) #analytical solution
    denominator = total_energy_gaussian(T, 0, 1)

    return neumerator/denominator 




##2b##

#solution description#
#calculate the efficiency values for each temperature values from 300-10,000K in intervals of 100k
#using the functions from part 2a. Use a for loop and append each value to a list. Plot these values as a function of temperature

xvalues = [x for x in range(300, 11000, 100)] #create x values 
yvalues = [] #create blank y values


for temp in xvalues:
    nu = efficiency(temp, 380e-9, 780e-9)
    yvalues.append(nu)

#plot values
plt.plot(xvalues, yvalues)
plt.title("Lightbulb Efficiency from T=300K - 10000K")
plt.xlabel("Filament Temperature (K)")
plt.ylabel("Efficiency")
plt.show()





##2c##
##solution description##
##use golden ratio (since we can see which interval the maxiumum is in from the 2b graph) to calculate the max 
#with a tolerance of 1k. Since it's designed for minimums, use -f(x) to calculate the maximum instead. 
#the interval is roughly 5000-7500K

def goldenratio(t1, t4, lambda1, lambda2): #returns the maximum temperature value within a range specified by t1-t4
    accuracy = 1000        # Required accuracy in K
    z = (1+np.sqrt(5))/2       # Golden ratio values

    #calculate x2 and x3 midinterval values 
    t2 = t4 - (t4-t1)/z
    t3 = t1 + (t4-t1)/z

    # Initial values of the function at the four points
    f1 = -efficiency(t1, lambda1, lambda2)
    f2 = -efficiency(t2, lambda1, lambda2)
    f3 = -efficiency(t3, lambda1, lambda2)
    f4 = -efficiency(t4, lambda1, lambda2)

    # Main loop of the search process
    while t4-t1>accuracy:
        if f2<f3:
            t4,f4 = t3,f3
            t3,f3 = t2,f2
            t2 = t4 - (t4-t1)/z
            f2 = -efficiency(t2, lambda1, lambda2)
        else:
            t1,f1 = t2,f2
            t2,f2 = t3,f3
            t3 = t1 + (t4-t1)/z
            f3 = -efficiency(t3, lambda1, lambda2)
    
    return 0.5 * (t1 + t4)

##estimate that the maximum is betweeen 5000-7500K from the 2b graph

maximum = goldenratio(5000, 7500, 380e-9, 780e-9)
print(maximum)





##2d##
##solution description##
#first plot the new efficiency graph to determine the intervals needed
##using the same function as from 2c but with a different range of wavelengths, calculate the maximum temperature that 
##corresponds to the max efficincy. then see if that value is lower than the melting point of tungsten to see if its possible
#to use the bulb at that temperature. 

xvalues = [x for x in range(300, 11000, 100)] #create x values 
yvalues = [] #create blank y values


for temp in xvalues:
    nu = efficiency(temp, 780e-9, 2250e-9)
    yvalues.append(nu)

#plot values
plt.plot(xvalues, yvalues)
plt.title("Lightbulb Efficiency from T=300K - 10000K")
plt.xlabel("Filament Temperature (K)")
plt.ylabel("Efficiency")
plt.show()

#calculate the maximum 
newmax = goldenratio(2000, 4000, 780e-9, 2250e-9)
print(newmax)



