import numpy as np
import matplotlib.pyplot as plt

#saving paths
saving_path = '/mnt/raid-cita/dzli/gb057/'

def phase0(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint', save = 'no',mask_lowI=0.0):
	import spec as sp
	conj = sp.readconj(name=name)
	phase = np.angle(conj)
#apply mask
	mask_indices = np.log10(np.absolute(conj)**2.)< mask_lowI
	print 'number of unmasked points',len(phase)*len(phase[0])-np.count_nonzero(mask_indices)
	phase[mask_indices] = 100.0

	if save !='no':
		phase.astype('float32').T.tofile(saving_path+'phase_mask'+str(mask_lowI)+'_'+name)
		print 'save phase to: '
		print (saving_path+'phase_mask'+str(mask_lowI)+'_'+name)
	return phase

def readphase0(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',phase_type='phase',mask_lowI = 0.0):
	mask=''
	if mask_lowI != 0:
		mask = 'mask'+str(mask_lowI)+'_'

	if phase_type == 'phase':
		filename = saving_path + 'phase_'+mask + name
	else:
		filename = saving_path + 'phase_diff_'+mask + name

	print 'reading: ',filename	
	phase = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	return phase

def phase0_LminusR (phasel,phaser,save='no',name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',mask_lowI=0.0):
	phase_diff = phasel - phaser
	if save != 'no':
		phase_diff.astype('float32').T.tofile(saving_path+'phase_diff_'+name)
		print 'save phase LL-RR to: '
		print (saving_path+'phase_diff_mask'+str(mask_lowI)+'_'+name)
	return phase_diff

def plot_phase (phase,pdf,name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint' ):
#	import math
#	phase*=180.0*3600.0/math.pi #transform into arcsec
	plt.figure(1)
	plt.imshow(phase, origin='lower',interpolation='none',aspect='auto',cmap='gray',vmin = -3.14, vmax=3.14)
	plt.xlim(300,1000)
	plt.ylim(5000,12000)
	plt.title(name)
	plt.colorbar()
	
	pdf.savefig(1)
	plt.clf() # clear figure for next plot
	return 		
