import numpy as np
import matplotlib.pyplot as plt

# 3b Parameters
epsilon = 1
deltax = 0.02
deltat = 0.005
beta = (deltat/deltax) * epsilon
xfinal = 2*np.pi
boundarycondition = 0.0

def initialcondition(x): #function that defines and computes the initial condition 
    return np.sin(x)

#takes in the final t value that the algorithm will solve to and runs the leapfrog algorithm. Returns an array of u values at every t and x.
def burger(tfinal, boundarycondition):
    xfinal = 2*np.pi

    total_x_indexes = int((xfinal)/deltax)  #index * delta x = the actual value we want to plug into u functions
    total_t_indexes = int((tfinal)/deltat) + 1 #index * delta t = the actual value we want to plug into u functions
    print(total_x_indexes)
    print(total_t_indexes)

    uarray = np.zeros((total_t_indexes, total_x_indexes)) #each row = u values at a given time t

    #load uarray with initial condition values (t boundary conditions)
    #for x in range(0, total_x_indexes):
    x = np.linspace(0, xfinal, int(xfinal/deltax))
    uarray[0][:] = initialcondition(x)

    # apply x boundary conditions
    for t in range(0, total_t_indexes):
        uarray[t][0] = boundarycondition
        uarray[t][-1] = boundarycondition

    # perform the forward euler step to start the leapfrogging
    for x in range(1, total_x_indexes - 1):
        uarray[1][x] = uarray[0][x] - (beta/2) * (uarray[0][x+1]**2 - uarray[0][x-1]**2)

    # leapfrog algorithm 
    for t in range(1, total_t_indexes - 2):
        for x in range(1, total_x_indexes - 1):
            uarray[t + 1][x] = uarray[t - 1][x] - (beta / 2) * (uarray[t][x + 1] ** 2 - uarray[t][x - 1] ** 2)
    
    return uarray


##3b##

#call the function and have it run until the tfinal value indicated 
uarray=burger(2, boundarycondition)

#pull out the different arrays that correspond to each time t needed and plot them           
x = np.linspace(0, xfinal, int(xfinal/deltax))
for i in [0, 0.5, 1, 1.5]:
    plt.plot(x, uarray[int(i/deltat), :])

plt.title("Burger Equation Solution")
plt.xlabel("position (m)")
plt.ylabel("u")
plt.legend(["t = 0s", "t = 0.5s", "t = 1s", "t = 1.5s"])
plt.xlim(0, 2*np.pi)
plt.show()


