import numpy as np
import matplotlib.pyplot as plt

def read():
	path = '/mnt/scratch-lustre/pen/VLBI/Gate0/'
	name='gb057_1.input_baseline257_freq_03_pol_LL.bin.all'

	path = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'
	name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'	
	filename = path + name
	data = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	print(data.shape)
	# frequency channels
	return data

def plot_spec(data,spec_type = 'dy',save = 'no',name = 'll'):
	if spec_type == 'dy':
		plt.figure(1)
		plt.imshow(data,origin='lower',interpolation='none',aspect='auto')
	elif spec_type == 'sec':
		plt.figure(2)
		plt.imshow(np.log10(np.absolute(data)**2.), origin='lower',interpolation='none',aspect='auto',cmap='gray')
		plt.xlim(300,1000)
		plt.ylim(5000,12000)
	
	plt.colorbar()

	if save == 'no':
		plt.show()
	else:
		plt.savefig(spec_type+name+'.pdf')
	return 

def sec_spec(data):
	from scipy.fftpack import fft2, fftshift

	Ic = fftshift(fft2(data))

	return Ic
