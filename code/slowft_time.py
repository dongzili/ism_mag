from PyScint.fourier import *

import numpy as np
import sys
from scipy.fftpack import fftshift, fft2, ifft2, ifftshift, fft, ifft
#import nufft

def main( baseline,freq,f0):
    print 'Band selected: ' + freq

    fnd = '/mnt/scratch-lustre/simard/b0834/data/dynamic_spectrum_'+str(baseline)+'_freq_' + freq + '_f0' + str(f0)+'.bin'

    if baseline=='258':
        dir = '/mnt/scratch-lustre/simard/b0834/data/gb057_1.input_baseline258_freq_'+freq+'_pol_all.rebint.1/'
    elif baseline=='514':
        dir = '/mnt/scratch-lustre/simard/b0834/data/Gb-Gb/freq'+freq+'/'
    elif baseline=='257':
        dir = '/mnt/scratch-lustre/simard/b0834/data/Ar-Ar/freq'+freq+'/'
    else:
        print 'Invalid baseline'

    if freq=='00':
        initial_frequency = 310.5 #MHz
    elif freq=='01':
        initial_frequency = 318.5 #MHz
    elif freq=='02':
        initial_frequency = 326.5 #MHz
    elif freq=='03':
        initial_frequency = 334.5 #MHz
    else:
        print 'Invalid frequency'

    num_rows = 16384
    num_columns = 660
    bandwidth = 8.0
    band_per_channel = bandwidth/num_rows
    frequency = np.arange(num_rows)*band_per_channel + initial_frequency
    timespan = 5700.0 #s
    time_resolution = timespan/num_columns 
    time = np.arange(num_columns)*time_resolution

    filer = 'weighted_average_0_1_2_6.bin'
    filei = 'weighted_average_0_i_1_i_2_i_6_i.bin'
    
    datar = np.fromfile(dir+filer,
                        dtype=np.float32).reshape(num_columns,num_rows).T
    if baseline=='258':
        datai = np.fromfile(dir+filei, 
                            dtype=np.float32).reshape(num_columns,num_rows).T
        I = datar + 1.0j*datai
        del datai
    else:
        I = datar
    del datar

    # Pad the frequency axis
    Ipad = np.pad(I,((0,num_rows),(0,0)),'constant',constant_values=0)
    I = Ipad
    del Ipad

    # Set zeros to the mean - move to after zero padding to include
    # zeros from padding
    mean = np.mean(I)*np.size(I)/(np.count_nonzero(I)) #mean excluding zeros
    for i in range(num_columns):
        if np.sum(I[:,i])==0.:
            I[:,i] = mean
    #for i in range(num_rows):
    #    if np.sum(I[i,:])==0.:
    #        I[i,:] = mean
    
    # Slow FT the time axis
    #Ic = np.ones(I.shape,dtype=complex)
    Iout = np.ones(I.shape,dtype=complex)
    for i in range(len(frequency)):
        x = I[i,:]
        N = x.shape[0]
        n = np.arange(N) - N/2
        k = n.reshape((N,1))*frequency[i]/f0
        M = np.exp(-2.0j*np.pi*k*n/N)
        y = np.dot(M,x)
        
        # Slow FT back
    
        k = n.reshape((N,1))
        M = np.exp(2.0j*np.pi*k*n/N)
        Iout[i,:] = 1./N*np.dot(M,y)

    # Add zero padding
    #I = np.pad(Iout,((0,num_rows),(0,0)),'constant',constant_values=0)
    # Set zeros from zero padding and from the slowft to the mean
    #mean = np.mean(I[nonzero(I)])
    #I = np.where(I==0.,mean,I)
      
    print 'Writing output to ' + fnd
    Iout.tofile(fnd)

    return

if __name__=='__main__':
    if len(sys.argv)>3:
        main( sys.argv[1], sys.argv[2],float(sys.argv[3]))
    else:
        print '\nslowft_time.py requires command line arguments:\n\nconvolution.py baseline freq f0\n\nbaseline can be any of 257 258 or 514\nfreq can be any of 00 01 02 or 03\n'
        # later, buid in ability to use 2 or 1 frequency bands (ie. merge bands)
