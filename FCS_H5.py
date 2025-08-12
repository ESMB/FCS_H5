import os
import numpy as np
import h5py
from fretbursts import loader
import matplotlib.pyplot as plt
import csv
import multipletau  # Ensure this package is installed
import scipy.optimize as opt
import math
import pandas as pd


# List to store paths for analysis
pathlist = []

# The root path where results will be stored
path_root = "/Users/Mathew/Documents/Current analysis/PBS_Ab/"

# Filenmane contains

filename='1nM'

# Append the folder path to the list for processing
pathlist.append("/Users/Mathew/Documents/Current analysis/PBS_Ab/")

bin_size=1

k_green=-11.692949698263954
w_green=3.735508891679675e-07 
k_red=1876731.7235220107
w_red=4.149518739269024e-07


# Function to list all HDF5 files in the specified folder
def list_h5_files(folder_path, text):
    """
    List all .h5 files in the specified folder that contain specific text in their filename.

    Parameters:
    folder_path (str): Path to the folder containing HDF5 files.
    text (str): Specific text to look for in the filenames.

    Returns:
    list: A list of file paths for .h5 files containing the specified text.
    """
    return [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.h5') and text in f]


# Function to extract photon data from an HDF5 file
def extract_photon_data(file_name, offset=0):
    with h5py.File(file_name, 'r') as f:
        # Load detector and timestamp data
        detectors = f['photon_data/detectors'][...]
        timestamps = f['photon_data/timestamps'][...] + offset

        # Load spectral channels
        spectral_ch1 = f['photon_data/measurement_specs/detectors_specs/spectral_ch1'][0]
        spectral_ch2 = f['photon_data/measurement_specs/detectors_specs/spectral_ch2'][0]

        # Separate photon data by channel
        channel_1_data = timestamps[detectors == spectral_ch1]
        channel_2_data = timestamps[detectors == spectral_ch2]

  

        return channel_1_data, channel_2_data, np.max(timestamps)
    
    
    # Function to create a time trace from timestamps
def create_time_trace(timestamps, bin_size,min_time,max_time):
    
    
    # Calculate number of bins based on range and bin size
    num_bins = (max_time - min_time) // bin_size + 1
    
    # Use np.bincount for efficient counting
    indices = ((timestamps - min_time) // bin_size).astype(int)
    trace = np.bincount(indices, minlength=num_bins)
    
    # Calculate corresponding time bins
    time_bins = min_time + bin_size * np.arange(num_bins)

    return trace, time_bins




   
                                                 
                                               # Autocorrelated the sample - comment out if you just want to look at the dye.
                                                                           

########################################### Autocorrelate ############################################
    
def autocorrelate():
   global c
   global new_c
   global new_d
  


   
        
   x=np.array(green,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   c=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_c=np.delete(c, (0), axis=0)                                                  # This deletes the first row which is at t = 0. 
   
       
   
   x=np.array(red,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   d=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_d=np.delete(d, (0), axis=0) 



def crosscorrelate():
   global c
   global new_c
   global new_d
   global new_e


   
        
   x=np.array(green,dtype=float)   
   y=np.array(red,dtype=float)                                             # Convert the csv columns to float values - Why does it not know?
   c=multipletau.correlate(x,y,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_e=np.delete(c, (0), axis=0)                                                  # This deletes the first row which is at t = 0. 
   
       
def fungreen(x,n,td):
    k=k_green                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitgreen():
    xdata=new_c[:, 0]  
    ydata=new_c[:,1]
    guesses=np.array([20,6e-5])
    (n_, td_), _ = opt.curve_fit(fungreen, xdata, ydata,guesses)
    params= opt.curve_fit(fungreen, xdata, ydata,guesses)
    
   
    
    y_fit = fungreen(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "-",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_c[:,1]))
    ax1.plot(xdata, y_fit, '-',color='green')
    
    
    print(("Green_N = %r \r") %params[0][0])
    print(("Green_td = %r \r") %params[0][1])

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print(("Green_D = %r \r") %Diff)
    
    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print(("Green_r = %r \r") %Rh)
    
    conf_volume=(math.pi**(3/2)*k_green*w_green**3)*1000 
    concentration=(params[0][0]/6.022e23)/conf_volume
    print(("Green_conc = %r \r") %concentration)
    
def funred(x,n,td):
    k=k_red                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitred():
    xdata=new_d[:, 0]  
    ydata=new_d[:,1]
    guesses=np.array([25,5e-5])
    (n_, td_), _ = opt.curve_fit(funred, xdata, ydata,guesses)
    params= opt.curve_fit(funred, xdata, ydata,guesses)
    
   
    
    y_fit = funred(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "-",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_d[:,1]))
    ax1.plot(xdata, y_fit, '-',color='red')
    
    
    print(("Red_N = %r \r") %params[0][0])
    print(("Red_td = %r \r") %params[0][1])

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print(("Red_D = %r \r") %Diff)

    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print(("Red_r = %r \r") %Rh)
    conf_volume=(math.pi**(3/2)*k_red*w_red**3)*1000 
    concentration=(params[0][0]/6.022e23)/conf_volume
    print(("Red_conc = %r \r") %concentration)

def fitcross():
    xdata=new_e[:, 0]  
    ydata=new_e[:,1]
    guesses=np.array([25,5e-5])
    (n_, td_), _ = opt.curve_fit(funred, xdata, ydata,guesses)
    params= opt.curve_fit(funred, xdata, ydata,guesses)
    
   
    
    y_fit = funred(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "-",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_e[:,1]))
    ax1.plot(xdata, y_fit, '-',color='black')
    
    
    print(("Red_N = %r \r") %params[0][0])
    print(("Red_td = %r \r") %params[0][1])

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print(("Red_D = %r \r") %Diff)

    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print(("Red_r = %r \r") %Rh)
    conf_volume=(math.pi**(3/2)*k_red*w_red**3)*1000 
    concentration=(params[0][0]/6.022e23)/conf_volume
    print(("Red_conc = %r \r") %concentration)
    




for folder_name in pathlist:
    # Get list of all HDF5 files
    file_list = list_h5_files(folder_name,filename)
    num=len(file_list)
    path=folder_name
    print(file_list)
    # Initialize lists to hold all channel data
    all_channel_1_data = []
    all_channel_2_data = []
    offset = 0
    
    # Loop through each file and extract data
    for file_name in file_list:
        channel_1_data, channel_2_data, max_timestamp = extract_photon_data(file_name, offset)
        all_channel_1_data.append(channel_1_data)
        all_channel_2_data.append(channel_2_data)
        
        # Update offset for the next file
        offset = max_timestamp + 1
    
    # Concatenate all data from different files
    all_channel_1_data = np.concatenate(all_channel_1_data)
    all_channel_2_data = np.concatenate(all_channel_2_data)
    
    # Maximum bins
    max_ch1 = np.max(all_channel_1_data)
    max_ch2 = np.max(all_channel_2_data)

    overall_max = max(max_ch1, max_ch2)
    
    # Create time traces for each channel
    trace_ch1, time_bins_ch1 = create_time_trace(all_channel_1_data, bin_size,0,overall_max)
    trace_ch2, time_bins_ch2 = create_time_trace(all_channel_2_data, bin_size,0,overall_max)
    
        
    green = trace_ch1                                                                       # This is where the red and green data will be saved. Makes an array. Data exists as two columns of 10 us bursts (green and red)
    red = trace_ch2               
    
    
    
     
        
    autocorrelate()   
    crosscorrelate()   
    fitgreen()
    fitred()
    fitcross()
  
