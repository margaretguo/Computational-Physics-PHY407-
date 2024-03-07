import numpy as np
import matplotlib.pyplot as plt

##Problem Description##
##Fit the data given using the linear algebra based method described in the lab manual, with the goal being to find Sinverse * Y and 
##calculate the eigenvalues and vectors of that, which would include the solution. Ensure they fit the constraint by checking that S = positive definite
# Then plot that fitted elipse on top of the datapoints themselves. 

##1a##
##generate X matrix, calculate its transpose, multiple those together to get S. Take inverse of S. Generate Y matrix, multiply with S
#calculate eigenvectors and values, pick the largest eigenvalue, find it's index and then pick the column eigenvector with the same index
#to be the value with the ellipse coeffcients. 
#  
# original data points from table 1
x = np.array([-38.04, -35.28, -25.58, -28.80, -30.06])
y = np.array([27.71, 17.27, 30.68, 31.50, 10.53])

X = np.zeros((5, 6), float)
width = X.shape[0]

# generate X matrix
for i in range(width):
    X[i] = np.array([x[i] * x[i], x[i] * y[i], y[i] * y[i], x[i], y[i], 1])

Xtranspose = np.transpose(X) #calculate X transpose 

S = np.matmul(Xtranspose, X) #calculate S matrix 

S_inv = np.linalg.inv(S)

Y = np.zeros((6, 6), float)

Y[0][2] = 2
Y[1][1] = -1/4
Y[2][0] = 2

E, V = np.linalg.eig(np.dot(S_inv, Y)) #calculate eigenvectors 

print(E)  #eigenvalue at index 3 has highest value. 

a = V[:,3] #isolate the corresponding eigenvector

#equation 2 coefficients 
A = a[0]
B = a[1]
C = a[2]
D = a[3]
F = a[4]
G = a[5]



##1b##
#create a grid for plotting the elipse, calculate ellipse values using the function parameters
# found in 1a and plot them, plot the data points on top as normal

#generate blank arrays for plotting 
xfit = np.linspace(-40, 5, 100)
yfit = np.linspace(-5, 35, 100)

X, Y = np.meshgrid(xfit, yfit) #create a 2d matrix for plotting the elipse 

F = a[0] * X * X + a[1] * X * Y + a[2] * Y * Y + a[3] * X + a[4] * Y + a[5]

plt.contour(X, Y, F, [0]) #calculate and plot the elipse on the matrix created
plt.plot(x, y, "ro") #plot data points
plt.title("Best Fit Eclipse and Data Points")
plt.xlabel("X axis (AU)")
plt.ylabel("Y axis (AU)")
plt.legend(["Data Points"])
plt.show()





