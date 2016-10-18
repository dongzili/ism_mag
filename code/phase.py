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
	print 'number of unmasked points',np.count_nonzero(np.logical_not(mask_indices))
	phase[mask_indices] = 0.0

	if save !='no':
		phase.astype('float32').T.tofile(saving_path+'phase_mask'+str(mask_lowI)+'_'+name)
		print 'save phase to: '
		print (saving_path+'phase_mask'+str(mask_lowI)+'_'+name)
	return phase

def readphase0(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',phase_type='phase',mask_lowI = 0.0):
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

#to make sure delta psi is in the range of -pi to pi
	renorm = phase_diff > np.pi
	phase_diff[renorm]-=2.0*np.pi
	renorm = phase_diff < -np.pi
	phase_diff[renorm]+=2.0*np.pi

	if save != 'no':
		phase_diff.astype('float32').T.tofile(saving_path+'phase_diff_mask'+str(mask_lowI)+'_'+name)
		print 'save phase LL-RR to: '
		print (saving_path+'phase_diff_mask'+str(mask_lowI)+'_'+name)
	return phase_diff

def plot_phase (phase,pdf,name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint',vrange = np.pi ):
	if vrange == 0:
		print 'vrange = 0, auto set vrange = pi'
		vrange = np.pi
#	import math
#	phase*=180.0*3600.0/math.pi #transform into arcsec

# mask data points with abnormal large delta psi
	mask_indices = abs(phase) > vrange
	phase[mask_indices] = 0.0
#	phase[np.logical_not(mask_indices)] = -100	

	plt.figure(1)
	plt.imshow(phase, origin='lower',interpolation='none',aspect='auto',cmap='seismic',vmin=-vrange,vmax=vrange)
	plt.xlim(300,1000)
	plt.ylim(5000,12000)
	plt.title(name)
	plt.colorbar()
	
	pdf.savefig(1)
	plt.clf() # clear figure for next plot
	return 		
