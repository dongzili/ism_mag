import numpy as np
import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)
import png

def __init__():
    return

def smooth( array, n0, n1 ):
    array_out = array
    for i in range(n0):
        array_out = array_out + np.roll(array, i, axis=0) + np.roll(array, -1*i, axis=0)
    for i in range(n1):
        array_out = array_out + np.roll(array, i, axis=1) + np.roll(array, -1*i, axis=1)
    array_out = array_out / (1. + 2.*n0 + 2.*n1 )
    return array_out


def plot_secondary(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.real(smooth(I1*np.conjugate(I2),5,0))
                  [num_rows/2:7*num_rows/10,:]),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           #vmin=18,vmax=30,
               cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')
    ax.autoscale(False)



def plot_secondary_full(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.real(smooth(I1*np.conjugate(I2),5,0))
                  ),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[0],delay[-1]],
           #vmin=18,vmax=30,                                                    
               cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')
    ax.autoscale(False)
    
def plot_secondary_zoom(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.real(smooth(I1*np.conjugate(I2),5,0))
                  [num_rows/2:6*num_rows/10,:len(doppler)/2]),origin='lower',
           extent=[doppler[0],doppler[len(doppler)/2],
                   delay[num_rows/2],delay[6*num_rows/10]],
           vmin=18,vmax=30,cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')
    ax.autoscale(False)    

def plot_dynamic(I,time,frequency):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(I,origin='lower',extent=[time[0],time[-1],frequency[0],frequency[-1]],
               #vmin=18,vmax=30,
               cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'time (s)')
    plt.ylabel(r'frequency (MHz)')

def plot_secondary_squared(I1,I2,doppler,delay,num_rows):
    # Same as plot_secondary_abs(I1,I2,doppler,delay,num_rows)
    # Delete this and change everything to plot_secondary_abs
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.absolute(smooth(I1*np.conjugate(I2),5,0))
                  [num_rows/2:7*num_rows/10,:]),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           vmin=18,vmax=30,cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')
    
def plot_secondary_abs(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.absolute(smooth(I1*np.conjugate(I2),5,0))
                  [num_rows/2:7*num_rows/10,:]),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           vmin=18,vmax=30,cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')    

def plot_secondary_imag(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    Ip = np.imag(smooth(I1*np.conjugate(I2),5,0))[num_rows/2:7*num_rows/10,:]
    plt.imshow(np.sqrt(np.sqrt(np.abs(Ip)))*np.sign(Ip),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')

def plot_secondary_real(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    Ip = np.real(smooth(I1*np.conjugate(I2),5,0))[num_rows/2:7*num_rows/10,:]
    plt.imshow(np.sqrt(np.sqrt(np.abs(Ip)))*np.sign(Ip),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')

def plot_secondary_diff(I1,I2,doppler,delay,num_rows):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(np.log(np.absolute(smooth(I1*np.conjugate(I1),5,0))
                  [num_rows/2:7*num_rows/10,:])-np.log(np.absolute(
                      smooth(I2*np.conjugate(I2),5,0))
                  [num_rows/2:7*num_rows/10,:]),origin='lower',
           extent=[doppler[0],doppler[-1],
                   delay[num_rows/2],delay[7*num_rows/10]],
           #vmin=18,vmax=30,
           cmap=plt.get_cmap('gray_r'))
    ax.set_aspect('auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel(r'$f_D$ [mHz]')
    plt.ylabel(r'$\tau$[ms]')    
    
def plot_png(I,fnd):
    I = (I - np.min(I))*255/(np.max(I)-np.min(I))
    f1 = open(fnd,'wb')
    w = png.Writer(I.shape[1],I.shape[0],greyscale=True)
    w.write(f1,I)
    f1.close()

