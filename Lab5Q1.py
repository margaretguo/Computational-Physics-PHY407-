import numpy as np
import matplotlib.pyplot as plt

##1b##
##implement the Lax Wendroff method dscribed in the lab manual using the parameters given in question 1b. 

#set up initial parameters
deltax = 0.02
deltat = 0.01
Hx = 0.01 # constant H(x) here
g = 9.81

#used to calculate initial eta
A = 0.002
mu = 0.5
sigma = 0.05

initialu = 0 # for t = 0 and all x, u = 0 
def initialeta(x): #for t = 0 and all x, eta = the number outputted by this gaussian function
    return A*np.exp(-(x-mu)**2/sigma**2)

#compute F(uarrow) with given u and eta values
def f(u, eta): #u and eta are values that make up uarrow together
    Hx = 0.01
    fu = 0.5 * u**2 + g * eta
    feta = (Hx + eta)*u
    return fu, feta

def feta(u, eta): #just calculates eta for use with boundary conditions
    Hx = 0.01
    feta = (Hx + eta)*u
    return feta


def laxwendroff(tfinal): #takes in a final t value and returns the eta array for all t and x values calculated
    #find our x and t index values
    x0 = 0
    t0 = 0
    xfinal = 1
    total_x_indexes = int((xfinal-x0)/deltax) #index * delta x = the actual value we want to plug into u and nu functions
    total_t_indexes = int((tfinal-t0)/deltat) + 1 #index * delta t = the actual value we want to plug into u and nu functions
    print(total_x_indexes)
    print(deltax)

    #make uarrow and fill with t = 0 eta and u values
    uarray = np.zeros((total_t_indexes, total_x_indexes)) #one row of eta and u for each time value 
    etaarray = np.zeros((total_t_indexes, total_x_indexes))
    for x in range(0, total_x_indexes):
        etaarray[0][x] = initialeta(x*deltax + (deltax/2))
    
    print(etaarray)

    #apply the lax-wendroff method to the remaining values 
    for t in range(0, total_t_indexes - 1):
        for x in range(0, total_x_indexes):

            #rule out boundary conditions #fix bc's later 
            if x == (total_x_indexes - 1):
                etanj = etaarray[t][x]
                etanjminusone = etaarray[t][x - 1] 
                fetanj = feta(uarray[t][x], etanj)      
                fetanjminusone = feta(uarray[t][x - 1], etanjminusone)

                uarray[t+1][x] = 0.0
                etaarray[t+1][x] = etanj - (deltat/deltax) *(fetanj - fetanjminusone)
                

            elif x == 0:
                
                #regular values
                etajplusone = etaarray[t][x + 1]
                etanj = etaarray[t][x]
                fetajplusone = feta(uarray[t][x+1], etajplusone)
                fetanj = feta(uarray[t][x], etanj)

                uarray[t+1][x] = 0.0
                etaarray[t+1][x] = etanj - (deltat/deltax) *(fetajplusone - fetanj) #forward backward difference formula 
                
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
                funjplusone, fetajplusone = f(unjplusone, etajplusone) #Fnj+1
                funjminusone, fetajminusone = f(unjminusone, etanjminusone) #Fnj-1
                funj, fetanj = f(unj, etanj) #Fnj

                #calculate uhalfsteps
                #uarrow n+1/2 j + 1/2
                ujplushalf = 0.5 * (unjplusone + unj) - 0.5*(deltat/deltax) * (funjplusone - funj)
                etajplushalf = 0.5 * (etajplusone + etanj) - 0.5*(deltat/deltax) * (fetajplusone - fetanj)

                #uarrow n+1/2 j - 1/2
                ujminushalf = 0.5 * (unjminusone + unj) - 0.5*(deltat/deltax) * (funj - funjminusone)
                etajminushalf = 0.5 * (etanjminusone + etanj) - 0.5*(deltat/deltax) * (fetanj - fetajminusone)

                #calculate fhalfsteps
                #f(uarrow) n+1/2 j + 1/2
                fujplushalf, fetajplushalf = f(ujplushalf, etajplushalf)

                #f(uarrow) n+1/2 j - 1/2
                fujminushalf, fetajminushalf = f(ujminushalf, etajminushalf)

                #calculate uarrow n + 1 j
                unplusonej = unj - (deltat/deltax) * (fujplushalf - fujminushalf)
                etanplusonej = etanj - (deltat/deltax) * (fetajplushalf - fetajminushalf)

                #append new values to uarrow
                uarray[t+1][x] = unplusonej
                etaarray[t+1][x] = etanplusonej 
        
    return etaarray


##1c##


#create x values for plotting
xvalues = np.arange(0, 1, deltax)

#generate y values using lax wendroff method
one_second = laxwendroff(1.0)
four_second = laxwendroff(4.0)

zero_second = one_second[0]
one_second = one_second[-1]
four_second = four_second[-1]

#plot values and label plot
plt.plot(xvalues, zero_second)
plt.plot(xvalues, one_second)
plt.plot(xvalues, four_second)
plt.title("position vs eta (height) at different times")
plt.xlabel("position (m)")
plt.ylabel("eta (m)")
plt.legend(["t = 0s", "t = 1s", "t = 4s"])
plt.show()


    

