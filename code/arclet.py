import spec as sp
import numpy as np
from scipy.fftpack import fftshift, fft2, ifft2, ifftshift, fft, ifft
nubox=8000

def __init__():
    return 

def parabola(x,x0,y0):
    return -y0/x0**2.0*(x-x0)**2.0 + y0

def main_parabola_x(y,x1,y1):
    return np.sqrt(x1**2./y1*y)

def main_parabola_y(x,x1,y1):
    return y1/x1**2. * x**2.

def get_conjugate_spectrum(I):
    Ic = fftshift(fft2(I))
    return Ic

def parabola_general_y(x,x0,y0,x1,y1):
    A = (y1-y0)/(x1-x0)**2.
    return A*(x-x0)**2. + y0

def parabola_general_x(y,x0,y0,x1,y1):
    A = (y1-y0)/(x1-x0)**2.
    return x0 + np.sqrt(np.absolute((y-y0)/A)), x0 - np.sqrt(np.absolute((y-y0)/A))

def convolved_conjugate_spectrum(Ic,doppler,delay,pdf,f0=0,tbox=5.7,nubox=800.):

	delay_max = 7*Ic.shape[0]/10 #+ 1177
	delay_min = np.argmin(np.absolute(delay-0.0999)) #+ 692
	doppler_min = 0 #int(660/2 - 206*(f0/314.5))
	doppler_max = Ic.shape[1]/2 #int(660/2 - 41*(f0/314.5))

	#define the thickness of the arclet
	delay1 = 0.148
	doppler1 = -20.52*(f0/338.5)
	#doppler1 = (-20.35-doppler_resolution*1.0)*f0/338.5
	delay2 = 0.138
	doppler2 = doppler1

	delay1=1.41
	delay2=1.36
	doppler1=-19.7
	doppler2=doppler1

	I = ifft2(ifftshift(Ic))
	#shift axis back to 0,0?
	delayi = np.argmin(abs(delay-np.average([delay1,delay2])))
	doppleri = np.argmin(abs(doppler-np.average([doppler1,doppler2])))

	#mask other points besides the moon shape thing
	mask = np.zeros(Ic.shape)
	mask[delay_min:delay_max,doppler_min:doppler_max]=1.0
	for i in range(doppler_min,doppler_max):
		delayu=parabola(doppler[i],doppler1,delay1)
		delayl=parabola(doppler[i],doppler2,delay2)
		mask[:,i] = np.where(delay>delayu, 0, mask[:,i])
		mask[:,i] = np.where(delay<delayl, 0, mask[:,i])
	sp.plot_spec(Ic*mask,pdf,spec_type='sec',tbox=tbox,nubox=nubox)
	Iarc = ifft2(ifftshift(Ic*mask))
	I2 = np.real(I)*np.exp(-1.0j*np.angle(Iarc))
	I2c = fftshift(fft2(I2))/len(I2)/len(I2[0])
	I2c = np.roll(np.roll(I2c,doppleri-len(doppler)/2,1),delayi-len(delay)/2,0)
	return I2c

def fourier_transform_incremental_width(I,frequency,f0):
    If = np.ones(I.shape)*(1.0 + 1.0j)
    next = 0.0
    for i in range(len(frequency)):
        #if i%len(frequency)>next:
        #    print str(next*100.) + ' % through at index ' + str(i)
        #    next += 0.1
        If[i,:] = slow_dft(I[i,:],frequency[i]/f0)
    Ic = fftshift(fft(If,axis=0))
    return Ic
    
def slow_dft(x,scaling):
    #x = np.asarray(x,dtype=float)
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N,1))
    Pk = N #np.max(k)-np.min(k)
    k = np.where(k>=Pk/2,k-Pk,k)*scaling
    k = np.where(k<0,k+Pk,k)
    M = np.exp(-2j*np.pi*k*n/N)
    return np.dot(M,x)

def read_data(dir,num_columns,num_rows,filer,filei):
    datar = np.fromfile(dir+filer,
                        dtype=np.float32).reshape(num_columns,num_rows).T
    datai = np.fromfile(dir+filei, 
                        dtype=np.float32).reshape(num_columns,num_rows).T
    
    I = datar + 1.0j*datai
    
    del datar
    del datai

    return I

def get_mask(num_rows,num_columns,doppler,delay,doppler1,doppler2,delay1,delay2,delaycut=0.07):
	print doppler1,doppler2,delay1,delay2
	mask = np.zeros((int(num_rows),int(num_columns)))
	delay_max = mask.shape[0]//2 + 2500 #np.int(5.6*mask.shape[0]/10 )
	delay_min = mask.shape[0]//2 + int(delaycut*nubox)#np.argmin(np.absolute(delay-0.0999))#16384/2 +559 #cut of at a delay of 0.07 ms                      
	doppler_min = 0 #mask.shape[1]/2 
	doppler_max = mask.shape[1]//2


	mask[delay_min:delay_max,doppler_min:doppler_max]=1.
	for i in range(mask.shape[1]):
		delayh=parabola(doppler[i],doppler1,delay1)
		delayl=parabola(doppler[i],doppler2,delay2)
		mask[:,i] = np.where(delay>delayh, 0, mask[:,i])
		mask[:,i] = np.where(delay<delayl, 0, mask[:,i])

	return mask

def wiener_deconvolution1(y,H):
    signal = np.absolute(y-np.mean(y))+np.absolute(1./H-np.mean(1./H))
    noise = 2.*np.average(np.absolute(y-np.mean(y)))
    F = np.absolute(H)/np.conjugate(H)
    W = np.absolute(F)**2./(np.absolute(F)**2.+noise/signal)
    Ic = fftshift(fft2(y/F*W))
    #Ic = np.roll(np.roll(Ic,doppleri-len(doppler)/2,1),delayi-len(delay)/2,0)
    return Ic

def wiener_deconvolution(y,H):
   signal = np.absolute(y-np.mean(y))
   noise = np.average(signal)
   F = np.absolute(H)/np.conjugate(H)
   W = np.absolute(F)**2./(np.absolute(F)**2.+noise/signal)
   Ic = fftshift(fft2(y/F*W))
   return Ic


def get_main_para(num_rows,num_columns,doppler,delay,doppler1,doppler2,delay1,delay2,delaycut=0.07,delayup=0.5):
	print doppler1,doppler2,delay1,delay2
	mask = np.zeros((int(num_rows),int(num_columns)))
	delay_max = mask.shape[0]//2 + int(delayup*nubox) #np.int(5.6*mask.shape[0]/10 )
	delay_min = mask.shape[0]//2 + int(delaycut*nubox)#np.argmin(np.absolute(delay-0.0999))#16384/2 +559 #cut of at a delay of 0.07 ms                      
	doppler_min = 0 #mask.shape[1]/2 
	doppler_max = mask.shape[1]


	mask[delay_min:delay_max,doppler_min:doppler_max]=1.
	for i in range(mask.shape[1]):
		delayh=main_parabola_y(doppler[i],doppler1,delay1)
		delayl=main_parabola_y(doppler[i],doppler2,delay2)
		mask[:,i] = np.where(delay>delayh, 0, mask[:,i])
		mask[:,i] = np.where(delay<delayl, 0, mask[:,i])

	return mask
