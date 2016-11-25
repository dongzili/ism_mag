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
#		mean_spec(
#
# First version: Dongzi Li. 2016.10.15
#===============================================================
import numpy as np
import matplotlib.pyplot as plt

#rebin spectrum
def rebin(a, shape):
	sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
	return a.reshape(sh).mean(-1).mean(1)

def readdy(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):

	path = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'
	filename = path + name
	data = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	print(data.shape)
	# frequency channels
	return data


def readconj(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',cross = 'no'):
	path = '/mnt/raid-cita/dzli/gb057/'
	if cross != 'no':
		name = name.replace('_pol_LL','')
		filenamer = path + 'conj_cross_r'+name
		filenamei = path + 'conj_cross_i'+name
	else:
		filenamer = path + 'conj_r_'+name
		filenamei = path + 'conj_i_'+name
		
	datar = np.fromfile(filenamer,dtype=np.float32) .reshape(-1,16384).T
	datai = np.fromfile(filenamei,dtype=np.float32) .reshape(-1,16384).T
	print('conj array:',datar.shape)
	data = datar + datai*1.0j
	del datar
	del datai
	return data

#========================================================
def conj_spec(data,save = 'no',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	from scipy.fftpack import fft2, fftshift
	Ic = fftshift(fft2(data))
	Ic = Ic/len(Ic)/len(Ic[0])
	if save == 'no':
		return Ic
	else:
		save_conj(Ic,name)
		return
#========================================================
def save_conj(Ic,name):
	print ('save conj_spec: ')
	path = '/mnt/raid-cita/dzli/gb057/'
	np.real(Ic).astype('float32').T.tofile(path+'conj_r_'+name)
	np.imag(Ic).astype('float32').T.tofile(path+'conj_i_'+name)
	print (path+'conj_r/i_'+name)
	del Ic
	return
	
#========================================================
def cross_sec(namel,namer,save = 'no',outname=''):
	print 'start cross conj'
	path = '/mnt/raid-cita/dzli/gb057/'
	icr = readconj(name=namer)
	crosssec = readconj(name=namel)*np.conj(icr)
	del icr
	if save != 'no':
		np.real(crosssec).astype('float32').T.tofile(path+'sec_cross_r_'+outname)
		np.imag(crosssec).astype('float32').T.tofile(path+'sec_cross_i_'+outname)
		print 'save to', path+'conj_cross_r/i'+outname
		crosssec = None
	return crosssec
#========================================================
def plot_spec(data, pdf, spec_type = 'dy',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint',cmap = 'Greys',show = 'plot',tbox=5.7,nubox=8000.,bins=1):
	print tbox,nubox
	print 'spec_type', spec_type
	fig=plt.figure(1)
	if spec_type == 'dy':
		fig=plt.figure(1)
		plt.imshow(data,origin='lower',interpolation='none',aspect='auto')
	else:
		# x: fs; y:tau
		dpix,dpiy = [300.,5000]
		cenx,ceny = [len(data[0])//2,len(data)//2]
		exten = [cenx-dpix,cenx+dpix,ceny,ceny+dpiy]
		data=data[exten[2]:exten[3],exten[0]:exten[1]]
		rebin_dpix,rebin_dpiy=[dpix,dpiy//bins]
		print data.shape
		print rebin_dpiy,rebin_dpix
		if bins!=1:
			data = rebin(data,(rebin_dpiy*2,rebin_dpix))
#		exten = [(exten[0]-cenx)/tbox,(exten[1]-cenx)/tbox,0,(ceny-exten[2])/nubox] 
		exten = [-dpix/tbox,(dpix-1)/tbox,0,(dpiy-1)/nubox] 
		print exten
		dpi = 100
		fig=plt.figure(1,figsize=(rebin_dpix*2/dpi,rebin_dpiy/dpi),dpi=dpi)

		
		if spec_type == 'sec':
			plt.imshow(np.log10(np.absolute(data)**2.), origin='lower',interpolation='none',aspect='auto',cmap= cmap,extent=exten)
		elif spec_type == 'conj':
			plt.imshow(np.log10(np.absolute(np.imag(data))), origin='lower',interpolation='none',aspect='auto',cmap= cmap,extent=exten)
			name = 'cross_sec_i'+name.replace('_pol_LL.rebint','')
		elif spec_type == 'phase':
			plt.imshow(np.sign(data)*abs(data)**0.5,origin='lower',interpolation='none',aspect='auto',cmap= 'Dark2',extent=exten)

#		plt.xlim(-50,0)
#		plt.ylim(0,0.5)

	plt.colorbar()
	plt.title(name)
	plt.xlabel(r'differential doppler frequency $f_D$ mHz')
	plt.ylabel(r'time delay $\tau$ ms')
	#for interactive mode	
	if show == 'interact':	
		def onclick(event):
   	 		print('xdata=%f, ydata=%f' %(event.xdata, event.ydata))
		cid = fig.canvas.mpl_connect('button_press_event', onclick)
		plt.show()
#	elif show=='plot':
	#	plt.show()
	else:
	#	plt.savefig(spec_type+name+'.pdf')
		pdf.savefig(1,figsize=(rebin_dpix*2/dpi,rebin_dpiy/dpi),dpi=dpi)
		print fig.get_size_inches()
#		plt.clf() # clear figure for next plot
	return 
#=============================================
#x is an array denote x axis
def plot_parabola(x,x0,y0,pdf):
	import arclet as arc
	y=np.zeros(len(x))
	for i in range(len(x)):
		y[i]=arc.parabola(x[i],x0,y0)
	posi = y>0
	x=x[posi]
	y=y[posi]
	plt.plot(x,y)
#	pdf.savefig(1)
#	plt.show()
	return 
#================================================
def mean_spec(data,box,cover='no'):
	print 'start mean_spec' 
	num_points = np.count_nonzero(data[box[2]:box[3]+1,box[0]:box[1]+1])
	[average] = sum(data[box[2]:box[3]+1,box[0]:box[1]+1].reshape(-1,1))*1.0/num_points
	if cover !='no':
		data[box[2]:box[3]+1,box[0]:box[1]+1] = average
	print 'number of points averaged: ', num_points, 'mean: %i log10 %.1f' %(np.sign(average), np.log10(np.abs(average)))
	return average,num_points

