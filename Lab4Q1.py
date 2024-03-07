import numpy as np
import matplotlib.pyplot as plt

##Problem Description##
#In question 1, the goal is to first plot the second column in the data set (closing values) on the y axis against the first column in the data set. 
#Then, run the real fourier transform on the second column of data and then run the inverse real fourier transform on that again. Compare the 
#results and see if that returns points that are similar to the original data (this is open ended, need to design my own test for this). Then, 
#assuming there are 21 business days per month, calculate how many coefficients in the real fourier transform of the data need to be set to zero 
#for data variations with a period of more than 6 months to be removed. Run the inverse fourier transform on that and then plot it against the data. 

##1a##
#Solution Description#
#using numpy array and textfile functions, load in the closing values and times into separate arrays and plot them. 

fulldata = np.loadtxt(
    "sp500.csv", delimiter=",", dtype=str
)  # load all data from csv file

closingvalues = fulldata[:, 1]  # select the second column of the data only
closingvalues = np.delete(
    closingvalues, 0
)  # remove the column label within the csv file
closingvalues = closingvalues.astype(float)  # convert array into a float

# generate the business day number array for plotting
number_of_days = len(closingvalues)  # determine the total nuber of days
business_day_number = np.arange(0, number_of_days)  # generate xvalues array

# plot closing value
plt.plot(business_day_number, closingvalues)
plt.title("Business Day Number vs S&P 500 Closing Value, 2014-2019")
plt.xlabel("Business Day Number")
plt.ylabel("S&P 500 Closing Value")
plt.show()

##1b##
#Solution Description#
#Run rfft on the closing values data generated in part a, plot that just to visually inspect (should see a plot with distinct spikes)
#Run irfft next and plot that on top of the data. 
#Compute the average percentage difference ((experimential - literature)/literature) of each term in the ifft and data as part of my test. 
#Compute maximum percentage difference to get a better sense of the distribution of the % differences. 
#plot the irfft and the data to visually inspect the difference as another part of my test. 

coefficients = np.fft.rfft(closingvalues)  # calculate the rfft

computed_closing_values = np.fft.irfft(coefficients) # calculate the irfft
numberofterms = len(computed_closing_values)

#plot rfft to ensure it was calculated properly 
plt.plot(np.abs(coefficients), "^")
plt.title("fft of S&P 500 Closing Values, 2014-2019")
plt.show()

# plot the irfft values on top of the original data for visual inspection test
plt.plot(closingvalues, "b")
plt.plot(computed_closing_values, "r--")
plt.title("ifft of fft of S&P 500 Closing Values, 2014-2019")
plt.legend(["original closing values", "ifft computed values"])
plt.xlabel("Business Day Number")
plt.ylabel("S&P 500 Closing Value")
plt.show()


# compute percentage differences of each term in ifft and original closing values
percentage_difference = 0
max_percentage_difference = 0
max_index = 0
for i in range(0, numberofterms):
    percentage_difference += np.abs(
        (computed_closing_values[i] - closingvalues[i]) / closingvalues[i]
    )
    # calculate maximum percentage difference
    max_percentage_difference = max(
        max_percentage_difference,
        np.abs((computed_closing_values[i] - closingvalues[i]) / closingvalues[i]),
    )

    if max_percentage_difference == np.abs((computed_closing_values[i] - closingvalues[i]) / closingvalues[i]):
        max_index = i
percentage_difference /= numberofterms  # calculate average percentage difference

#print % difference values 
print(percentage_difference * 100)
print(max_percentage_difference * 100)
print(max_index)


##1c##
#Solution Description#
#manually calculate how many coefficients a period of 126 business days or less correspond to using w = 2pi/126 and w = 2pi*n/1259 to solve for n. 
#obtain n = 10 as the coefficient that corresponds to a period of 126, so that coefficient and any of the ones that come after it need to be 
#set to 0. 
#set these values to 0 using a for loop, run irfft on the coefficients and then plot that on top of the original data. 

#set appropriate coefficinets to 0
length_of_coefficients = len(coefficients)
low_pass_coefficients = coefficients
for i in range(9, length_of_coefficients):
    low_pass_coefficients[i] = 0.0

#run the irfft on lowpassed fft 
print(low_pass_coefficients)
low_pass_data = np.fft.irfft(low_pass_coefficients)

low_pass_data_length = len(low_pass_data)
xvalues = np.arange(0, low_pass_data_length)

#plot data and irfft output
plt.plot(business_day_number, closingvalues)
plt.plot(xvalues, low_pass_data)
plt.title("Business Day Number vs S&P 500 Closing Value, 2014-2019")
plt.legend(["original data", "6 month low-pass filtered data"])
plt.xlabel("Business Day Number")
plt.ylabel("S&P 500 Closing Value")
plt.show()
