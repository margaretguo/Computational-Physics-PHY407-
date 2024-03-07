import numpy as np
import matplotlib.pyplot as plt

##Problem Description##

##2a##
##Solution Description##

#set up initial parameters
deltax = 0.005
deltat = 0.001
g = 9.81

#returns topography value for a given x 
def Hx(x):
    a = 0.001
    return a + 0.1*np.exp(-7*x)

#used to calculate initial eta
A = 0.0005
mu = 0
sigma = 0.1

initialu = 0 # for t = 0 and all x, u = 0 
def initialeta(x): #for t = 0 and all x, eta = the number outputted by this gaussian function
    return A*np.exp(-(x-mu)**2/sigma**2)

#compute F(uarrow) with given u and eta values
def f(u, eta, x): #u and eta are values that make up uarrow together
    Hvalue = Hx(x * deltax)
    fu = 0.5 * u**2 + g * eta
    feta = (Hvalue + eta)*u
    return fu, feta

def feta(u, eta, x): #just calculates eta for use with boundary conditions
    Hvalue = Hx(x * deltax) 
    feta = (Hvalue + eta)*u
    return feta

def fhx(u, eta, x1, x2): #calculates the f(uarrow) value for half steps that need averaging of H(x)
    Hx1 = Hx(x1 * deltax)
    Hx2 = Hx(x2 * deltax)
    Hxmean = (Hx1 + Hx2) / 2 #calculate the mean between these two points as the H(x)
    fu = 0.5 * u**2 + g * eta
    feta = (Hxmean + eta)*u
    return fu, feta


def laxwendroffhx(tfinal): #takes in a final t value and returns the eta array for all t and x values calculated
    #find our x and t index values
    x0 = 0
    t0 = 0
    xfinal = 1
    total_x_indexes = int((xfinal-x0)/deltax) #index * delta x = the actual value we want to plug into u and nu functions
    total_t_indexes = int((tfinal-t0)/deltat) #index * delta t = the actual value we want to plug into u and nu functions
    print(total_x_indexes)
    print(deltax)

    #make uarrow and fill with t = 0 eta and u values
    uarray = np.zeros((total_t_indexes+1, total_x_indexes + 1)) #one row of eta and u for each time value 
    etaarray = np.zeros((total_t_indexes+1, total_x_indexes + 1))
    for x in range(0, total_x_indexes):
        etaarray[0][x] = initialeta(x*deltax)

    #apply the lax-wendroff method to the remaining values 
    for t in range(0, total_t_indexes):
        for x in range(0, total_x_indexes):

            #rule out boundary conditions #fix bc's later 
            if x == (total_x_indexes - 1):
                etanj = etaarray[t][x]
                etanjminusone = etaarray[t][x - 1] 
                fetanj = feta(uarray[t][x], etanj, x)      
                fetanjminusone = feta(uarray[t][x - 1], etanjminusone, x-1)

                etaarray[t+1][x] = etanj - (deltat/deltax) * (fetanj - fetanjminusone)
                uarray[t+1][x] = 0.0

            elif x == 0:
                
                #regular values
                etajplusone = etaarray[t][x + 1]
                etanj = etaarray[t][x]
                fetajplusone = feta(uarray[t][x+1], etajplusone, x+1)
                fetanj = feta(uarray[t][x], etanj, x)

                etaarray[t+1][x] = etanj - (deltat/deltax) * (fetajplusone - fetanj) #forward backward difference formula 
                uarray[t+1][x] = 0.0
            
            else:
                #calculate u 
                unjplusone = uarray[t][x+1]
                unj = uarray[t][x]
                unjminusone = uarray[t][x-1] 
                #calculate eta
                etajplusone = etaarray[t][x+1]
                etanj = etaarray[t][x]
                etanjminusone = etaarray[t][x-1]

                #calculate F(arrow values)
                funjplusone, fetajplusone = f(unjplusone, etajplusone, x+1) #Fnj+1
                funjminusone, fetajminusone = f(unjminusone, etanjminusone, x-1) #Fnj-1
                funj, fetanj = f(unj, etanj, x) #Fnj

                #calculate uhalfsteps
                #uarrow n+1/2 j + 1/2
                ujplushalf = 0.5 * (unjplusone + unj) - (deltat/(2*deltax)) * (funjplusone - funj)
                etajplushalf = 0.5 * (etajplusone + etanj) - (deltat/(2*deltax)) * (fetajplusone - fetanj)

                #uarrow n+1/2 j - 1/2
                ujminushalf = 0.5 * (unjminusone + unj) - (deltat/(2*deltax)) * (funj - funjminusone)
                etajminushalf = 0.5 * (etanjminusone + etanj) - (deltat/(2*deltax)) * (fetanj - fetajminusone)

                #calculate fhalfsteps
                #f(uarrow) n+1/2 j + 1/2

                fujplushalf, fetajplushalf = fhx(ujplushalf, etajplushalf, x+1, x)

                #f(uarrow) n+1/2 j - 1/2
                fujminushalf, fetajminushalf = fhx(ujminushalf, etajminushalf, x-1, x)

                #calculate uarrow n + 1 j
                unplusonej = unj - (deltat/deltax) * (fujplushalf - fujminushalf)
                etanplusonej = etanj - (deltat/deltax) * (fetajplushalf - fetajminushalf)

                #append new values to uarrow
                uarray[t+1][x] = unplusonej
                etaarray[t+1][x] = etanplusonej 
        
    return etaarray

##2b##

##Solution Description##

#create x values for plotting
xvalues = np.arange(0, 1+deltax, deltax)

#generate y values using lax wendroff method
one_second = laxwendroffhx(1.0)
two_second = laxwendroffhx(2.0)
three_second = laxwendroffhx(3.0)
four_second = laxwendroffhx(4.0)

zero_second = one_second[0]
one_second = one_second[-1]
two_second = two_second[-1]
three_second = three_second[-1]
four_second = four_second[-1]

#plot values and label plot
plt.plot(xvalues, zero_second)
plt.plot(xvalues, one_second)
plt.plot(xvalues, two_second)
plt.plot(xvalues, three_second)
plt.plot(xvalues, four_second)
plt.title("position vs eta (height) of a 1D tsunami wave at different times")
plt.xlabel("position (m)")
plt.ylabel("eta (m)")
plt.legend(["t = 0s", "t = 1s", "t = 2s", "t = 3s", "t = 4s"])
plt.show()