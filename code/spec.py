# a script to read dynamic/conjugate spectra, and plot the spectra 
#=================================================================
# Paths:
#	read in path = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'
#   save to path = '/mnt/raid-cita/dzli/gb057/'

# Functions Included:
#		readdy(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
#			read dynamic spectra, return it with freq*time

#		readconj(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
#			read conjugate spectra: complex, return it with freq*time		

#		conj_spec(data,save = 'no',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
#			calculate conjugate spectra from dynamic spectra,
#			save as float32, time * frequency
#			without correcting unit

#		plot_spec(data, pdf, spec_type = 'dy',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
#			plot multiple spectra into a pdf file
#			Arguments:
#					spec_type: 'dy', 'sec'; secondary spectra are plotted in log format
#					pdf: a pointer pointing to the targeted pdf
#
# First version: Dongzi Li. 2016.10.15
#===============================================================
import numpy as np
import matplotlib.pyplot as plt

def readdy(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):

	path = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'
	filename = path + name
	data = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	print(data.shape)
	# frequency channels
	return data


def readconj(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	path = '/mnt/raid-cita/dzli/gb057/'
	filename = path + 'conj_r_'+name
	datar = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	filename = path + 'conj_i_'+name
	datai = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	print('conj array:',datar.shape)
	data = datar + datai*1.0j
	del datar
	del datai
	return data

#========================================================
def conj_spec(data,save = 'no',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	from scipy.fftpack import fft2, fftshift
	Ic = fftshift(fft2(data))
	if save == 'no':
		return Ic
	else:
		path = '/mnt/raid-cita/dzli/gb057/'
		np.real(Ic).astype('float32').T.tofile(path+'conj_r_'+name)
		np.imag(Ic).astype('float32').T.tofile(path+'conj_i_'+name)
		print ('save conj_spec to: ')
		print (path+'conj_r/i_'+name)
		del Ic
		return

#========================================================
def plot_spec(data, pdf, spec_type = 'dy',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	if spec_type == 'dy':
		plt.figure(1)
		plt.imshow(data,origin='lower',interpolation='none',aspect='auto')
	elif spec_type == 'sec':
		plt.figure(1)
		plt.imshow(np.log10(np.absolute(data)**2.), origin='lower',interpolation='none',aspect='auto',cmap='gray')
		plt.xlim(300,1000)
		plt.ylim(5000,12000)
		plt.title(name)
	
	plt.colorbar()

	#	plt.savefig(spec_type+name+'.pdf')
	pdf.savefig(1)
	plt.clf() # clear figure for next plot
	return 

