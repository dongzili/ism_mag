# plot LL-RR of dynamic and secondary spectrum to see if there is position displacements
#secondary spectra is plotted with log(LL)-log(RR)

#first version Dongzi Li 2016.10.16
#====================================================
import spec as sp
import numpy as np
import matplotlib.pyplot as plt

def check_lr(pdf, spec_type = 'dy',namel =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint', namer =  'gb057_1.input_baseline257_freq_00_pol_RR.rebint'):

	if spec_type == 'dy':
		datal = sp.readdy(name=namel)
		datar = sp.readdy(name=namer)
		diff = datal-datar
		del datal,datar
		plt.figure(1)
		plt.imshow(diff,origin='lower',interpolation='none',aspect='auto')
		plt.title(namel)

	elif spec_type == 'sec':
		datal = sp.readconj(name=namel)
		datar = sp.readconj(name=namer)
		diff = np.log10(np.absolute(datal)**2.) - np.log10(np.absolute(datar)**2.)
		del datal,datar
		plt.figure(1)
		plt.imshow(diff, origin='lower',interpolation='none',aspect='auto',cmap='gray')
		plt.xlim(300,1000)
		plt.ylim(5000,12000)
		plt.title(namel)
	
	plt.colorbar()

	#	plt.savefig(spec_type+name+'.pdf')
	pdf.savefig(1)
	plt.clf() # clear figure for next plot
	return 
