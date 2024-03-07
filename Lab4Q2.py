import numpy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
##Problem Description##
#The goal of the problem is to calculate the amplitude, period and phase of different martian pressure datasets. First, calculate the correction/scaling
#values needed to convert the fft outputs into amplitude, period and phase. Then, using those relations, write a function that returns those values
#for each mode/coefficient returned from a given pressure and time dataset. Test this function using a premade cosine to see if it calculates the values
#correctly. Then, for pressure data columns A, C and D, run the function and use it to help plot the phase of and the amplitude as a function of period
#(on a log scale) on two different plots each. Find the peaks in each (one in A, 4 in C, none needed for D) and record their phases, amplitudes
#and periods. For data column E, run the function and just plot log scale of period vs amplitude. Find 4 peaks below 1 sol in period and 4 peaks
#above 100 sol and record their periods and amplitudes. 

##2a##
#Solution Description#
#run rfft to obtain coefficients (using rfft for increased speed/efficiency)
#determine analytically that the scale factor for amplitude = absolutevalue(coefficient) * 2/length of dataset
#determine that phase = arctan(imag/real) / 2pi
#determine that period = length of dataset/coefficient index * (conversion value)
#determine that period needs to be multiplied by the time difference between each measurement to convert it to real time units
#calculate these three values for each coefficient using a for loop and saving them to a list
#return a list of the amplitudes for all coefficients, list of all periods and a list of all phases. 
#add in a NaN value for the 0th period since it would require dividing by 0. 

#use a modified version of code provided in Lab Manual 4 to generate a cosine wave of known amplitude, period, phase to check if function
#is calculating correctly. 


# load in all data
def marsfft(times, pressures): #takes in arrays of time and pressure, returns 3 lists of amplitude, period, phase of each coefficient
    coefficients = fft.rfft(pressures)
    plt.plot(np.abs(coefficients), "^")
    plt.title("FFT Plot (marsfft() output)")
    plt.xlabel("Mode Number")
    plt.ylabel("Amplitude (unscaled)")
    plt.show()

    # calculate time conversion for period list
    conversion = times[2] - times[1]

    # generate lists to hold values for each mode
    amplitude_list = []
    phase_list = []
    period_list = []

    for i in range(0, len(coefficients)):
        # calculate values based on tested formulas
        amplitude = np.abs(coefficients[i]) * (
            2 / len(pressures)
        )  # scale factor = 2/number of terms in our data
        phase = np.arctan2(np.imag(coefficients[i]), np.real(coefficients[i])) / (
            2 * np.pi
        )  # scale factor to convert it to number of cycles = 1/2*pi
        if i == 0:
            period = float("NaN") # value in reality does not exist since dividing by 0
        else:
            period = (len(pressures) / i) * conversion

        # append values to list
        amplitude_list.append(amplitude)
        period_list.append(period)
        phase_list.append(phase)

    return amplitude_list, phase_list, period_list


# test out function by generating cosine
def cosine_wave(time, amplitude, period, phase):
    """Generate a sine wave with known amplitude, period, and phase"""
    return amplitude * np.cos((time / period + phase) * 2 * np.pi)


t = np.arange(0, 200, 1)
y = cosine_wave(t, 2, 100, 1 / 2)


transform = marsfft(t, y)
cos_amp = transform[0]
cos_phase = transform[1]
cos_period = transform[2]

# function hits maxiumum at n = 2 based on fourier plot. 
print(cos_amp[2], cos_phase[2], cos_period[2])


##2b##
#Solution Description#
#Load Pressure column A and the times into two arrays using method specified in 1a. Call fft function from 2a. 
#Using scipy.find_peaks() with no peakfinding parameters needed, obtain the index value of the peak and find the value corresponding 
#to that index in the amplitude, phase and period arrays. Print values. 
#plot phase and amplitudes as functions of log scale period.  

# load in data
fulldata = np.loadtxt("msl_rems.csv", delimiter=",", dtype=str)
times = fulldata[:, 0]
pressure_a = fulldata[:, 2]

# remove the column label within the csv file
times = np.delete(times, 0)
times = times.astype(float)  # convert array into a float
pressure_a = np.delete(pressure_a, 0)
pressure_a = pressure_a.astype(float)  # convert array into a float

# call fft function
values = marsfft(times, pressure_a)
amplitude_a = values[0]
phase_a = values[1]
period_a = values[2]

#find amplitude, period and phase of the peak and print
max_amp_ind_a , _ = find_peaks(amplitude_a)
max_amp_ind_a = max_amp_ind_a[0]
period_value_a = period_a[max_amp_ind_a]
amp_value_a = amplitude_a[max_amp_ind_a]
phase_value_a = phase_a[max_amp_ind_a]

print("Dataset A: Peak Period, Amp, Phase")
print(period_value_a)
print(amp_value_a)
print(phase_value_a)


#generate plots
plt.plot(period_a, amplitude_a, "r") #log period vs amplitude
plt.xscale("log") 
plt.ylim(0, 1.5)
plt.title("Period vs Amplitude for Dataset A")
plt.xlabel("Period (Sol)")
plt.ylabel("Amplitude (Pa)")
plt.show()

plt.plot(period_a, phase_a, "r") #log period vs phase
plt.xscale("log")
plt.title("Period (log scale) vs Phase for Dataset A")
plt.xlabel("Period (Sol)")
plt.ylabel("Phase (Sol)")
plt.show()


##2c##
#Solution Description#
#Load Pressure column C and the times into two arrays using method specified in 1a. Call fft function from 2a. 
#Using scipy.find_peaks() with no peakfinding parameters needed, obtain the index value of the 4 peaks and find the values corresponding 
#to those indexes in the amplitude, phase and period arrays. Print those three arrays of peak values. 
#plot phase and amplitudes as functions of log scale period. Plot peaks onto those plots. Also plot phase using discrete points instead of the 
#default line for ease of analysis. Keeping amplitude plot as a line helps with peak finding so did not plot amplitude using discrete points. 

# load in data
fulldata = np.loadtxt("msl_rems.csv", delimiter=",", dtype=str)
times = fulldata[:, 0]
pressure_c = fulldata[:, 4]

# remove the column label within the csv file
times = np.delete(times, 0)
times = times.astype(float)  # convert array into a float
pressure_c = np.delete(pressure_c, 0)
pressure_c = pressure_c.astype(float)  # convert array into a float

# call fft function
values = marsfft(times, pressure_c)
amplitude_c = values[0]
phase_c = values[1]
period_c = values[2]

#find amplitude, period and phase and print
period_values_c = []
amp_values_c = []
phase_values_c = []

#Use find_peaks() to find the peaks value
max_amp_ind_c , _ = find_peaks(amplitude_c)
print(max_amp_ind_c)
for i in max_amp_ind_c:
    period_value_c = period_c[i]
    amp_value_c = amplitude_c[i]
    phase_value_c = phase_c[i]

    period_values_c.append(period_value_c)
    amp_values_c.append(amp_value_c)
    phase_values_c.append(phase_value_c)

#print out period, amp, phase of the peaks 
print("Dataset C: Peak Period, Amp, Phase")
print(period_values_c)
print(amp_values_c)
print(phase_values_c)

#generate plots
plt.plot(period_c, amplitude_c, "g") #log period vs amplitude
plt.plot(period_values_c, amp_values_c, "^") #peaks
plt.xscale("log") 
plt.title("Period vs Amplitude for Dataset C")
plt.legend(["amplitude", "peak amplitudes"])
plt.xlabel("Period (Sol)")
plt.ylabel("Amplitude (Pa)")
plt.show()

plt.plot(period_c, phase_c, "g.") #log period vs phase
plt.plot(period_values_c, phase_values_c, "^") #peaks 
plt.xscale("log")
plt.title("Period vs Phase for Dataset C")
plt.legend(["phase", "peak phases"])
plt.xlabel("Period (Sol)")
plt.ylabel("Phase (Sol)")
plt.show()


##2d##
#solution Description#
#Load Pressure column D and the times into two arrays using method specified in 1a. Call fft function from 2a. 
#plot phase and amplitudes as functions of log scale period. Also plot phase using discrete points instead of the 
#default line for ease of analysis. Keeping amplitude plot as a line helps with peak finding so did not plot amplitude using discrete points.

# load in data
fulldata = np.loadtxt("msl_rems.csv", delimiter=",", dtype=str)
times = fulldata[:, 0]
pressure_d = fulldata[:, 5]

# remove the column label within the csv file
times = np.delete(times, 0)
times = times.astype(float)  # convert array into a float
pressure_d = np.delete(pressure_d, 0)
pressure_d = pressure_d.astype(float)  # convert array into a float

# call fft function, convert period
values = marsfft(times, pressure_d)
amplitude_d = values[0]
phase_d = values[1]
period_d = values[2]

#generate plots
plt.plot(period_d, amplitude_d, "b") #log period vs amplitude
plt.xscale("log") 
plt.title("Period vs Amplitude for Dataset D")
plt.xlabel("Period (Sol)")
plt.ylabel("Amplitude (Pa)")
plt.show()

plt.plot(period_d, phase_d, "b^") #log period vs phase
plt.xscale("log")
plt.title("Period vs Phase for Dataset D")
plt.xlabel("Period (Sol)")
plt.ylabel("Phase (Sol)")
plt.show()


##2e##
#Solution Description#
#Load Pressure column E and the times into two arrays using method specified in 1a. Call fft function from 2a. Split amplitude and period output arras
#into two sets of arrays, one containg amplitudes and periods with a period lower than 1.5 sol (Martian Tides) and one
#containing amps and periods with period higher than 100 sol (Seasons) using for loop and if statement. 
#Use find_peaks() again to find the 4 peaks in each of those two sets of data and printed. Needed to add a height parameter for the Seasons dataset
#and height + distance parameter for the Tides data set. Inspected the peakfinding quality by plotting the peak values on top of the data. 
#Due to the Tides data set having a very short peak, it was extremely difficult to prevent noise interference so the approximate peak values 
#were found visually and used to select the real peaks out of the printed array. 
#Plot amplitude on a log period scale. Then plotted peaks on top. 

# load in data
fulldata = np.loadtxt("msl_rems.csv", delimiter=",", dtype=str)
times = fulldata[:, 0]
pressure_e = fulldata[:, 6]

# remove the column label within the csv file
times = np.delete(times, 0)
times = times.astype(float)  # convert array into a float
pressure_e = np.delete(pressure_e, 0)
pressure_e = pressure_e.astype(float)  # convert array into a float

# call fft function
values = marsfft(times, pressure_e)
amplitude_e = values[0]
period_e = values[2]

#calculate mode peaks

#separate data into sections that meet tide and seasonal period criteria
seasonal_amps = []
seasonal_periods = []
tide_amps = []
tide_periods = []

for i in range(len(period_e)):
    if period_e[i] < 1.5:
        tide_amps.append(amplitude_e[i])
        tide_periods.append(period_e[i])
    if period_e[i] > 100.0:
        seasonal_amps.append(amplitude_e[i])
        seasonal_periods.append(period_e[i])



#run find_peaks to determine peaks in tides region
period_values_tide = []
amp_values_tide = []


max_amp_ind_tide , _ = find_peaks(tide_amps, height=1.2, distance= 1) #returns indexes of peaks
print(max_amp_ind_tide)
for i in max_amp_ind_tide:
    period_value_tide = tide_periods[i]
    amp_value_tide = tide_amps[i]

    period_values_tide.append(period_value_tide)
    amp_values_tide.append(amp_value_tide)
   
#print peaks
print("Dataset E (Tides): Period, Amp")
print(period_values_tide)
print(amp_values_tide)


#run find_peaks to determine peaks in seasons array
period_values_seasonal = []
amp_values_seasonal = []

max_amp_ind_seasonal , _ = find_peaks(seasonal_amps, height=5) #returns indexes of peaks
print(max_amp_ind_seasonal)
for i in max_amp_ind_seasonal:
    period_value_seasonal = seasonal_periods[i]
    amp_value_seasonal = seasonal_amps[i]

    period_values_seasonal.append(period_value_seasonal)
    amp_values_seasonal.append(amp_value_seasonal)
   
#print peaks
print("Dataset E (Seasonal Cycles): Period, Amp ")
print(period_values_seasonal)
print(amp_values_seasonal)


#generate plot
plt.plot(period_e, amplitude_e, "m") #log period vs amplitude
plt.plot(period_values_seasonal, amp_values_seasonal, "b^")
plt.plot(period_values_tide, amp_values_tide, "g^")
plt.xscale("log") 
plt.title("Period vs Amplitude for Dataset E")
plt.legend(["amplitude", "seasonal peaks", "tide peaks (with false peaks)"])
plt.xlabel("Period (Sol)")
plt.ylabel("Amplitude (Pa)")
plt.show()
