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

def readdy(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',path = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'):
	filename = path + name
	data = np.fromfile(filename,dtype=np.float32) .reshape(-1,16384).T
	print(data.shape)
	# frequency channels
	return data


def readconj(name = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint',cross = 'no',path = '/mnt/raid-cita/dzli/gb057/',num_t=1321,num_nu=16384):
	if cross != 'no':
		name = name.replace('_pol_LL','')
		filenamer = path + 'conj_cross_r'+name
		filenamei = path + 'conj_cross_i'+name
	else:
		filenamer = path + 'conj_r_'+name
		filenamei = path + 'conj_i_'+name
		
	datar = np.fromfile(filenamer,dtype=np.float32) .reshape(num_t,num_nu).T
	datai = np.fromfile(filenamei,dtype=np.float32) .reshape(num_t,num_nu).T
	print('conj array:',datar.shape)
	data = datar + datai*1.0j
	del datar
	del datai
	return data

#========================================================
def conj_spec(data,save = 'no',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	'start calculate conj'
	from scipy.fftpack import fft2, fftshift
	Ic = fftshift(fft2(data))
	if save == 'no':
		return Ic
	else:
		save_conj(Ic,name)
		return
#========================================================
def save_conj(Ic,name,type='float32',path = '/mnt/raid-cita/dzli/gb057/Gate0/'):
	print ('save conj_spec: ')
	np.float32(np.real(Ic)).astype('float32').T.tofile(path+'conj_r_'+name)
	np.float32(np.imag(Ic)).astype('float32').T.tofile(path+'conj_i_'+name)
	print (path+'conj_r/i_'+name)
	return
#=======================================================
def save_dy(Ic,name,type='float32',path = '/mnt/raid-cita/dzli/gb057/Gate0/'):
	print ('save dy_spec: ')
	if np.mean(np.imag(Ic))!=0:
		np.float32(np.real(Ic)).astype('float32').T.tofile(path+'r_'+name)
		np.float32(np.imag(Ic)).astype('float32').T.tofile(path+'i_'+name)
		print (path+'r/i_'+name)
	else:
		np.float32(Ic).astype('float32').T.tofile(path+name)
		print (path+name)
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
def select_rebin(data,dpix,dpiy,xbins,ybins,ypositive=0):
        cenx,ceny = [len(data[0])//2,len(data)//2]
        exten = [cenx-dpix,cenx+dpix,ceny-dpiy,ceny+dpiy]
        print 'select range: ',exten
        print 'original data size: ',data.shape
        data=data[exten[2]:exten[3],exten[0]:exten[1]]
        rebin_dpix,rebin_dpiy=[dpix//xbins,dpiy//ybins]
        print 'truncated data size: ',data.shape
        print 'rebint data size: ',rebin_dpiy*2,rebin_dpix*2
        if xbins*ybins!=1:
            data = rebin(data,(rebin_dpiy*2,rebin_dpix*2))
        if ypositive==1:
            data = data[len(data)//2:len(data),:]
        return data
#========================================================
def plot_spec(data, pdf=0, spec_type = 'dy',name =  'gb057_1.input_baseline257_freq_00_pol_LL.rebint',cmap = 'jet',show = 'plot',tbox=5.7,nubox=8000.,xbins=1,ybins=1,dpix=350,dpiy=5000,dpi=100,vmin=0,vmax=0):
    print 'total time, frequency',tbox,nubox
    print 'spec_type', spec_type
    if show!='interact':
        plt.ioff()
    if spec_type == 'dy':
        fig=plt.figure(1)
        exten = [0,tbox*1000.,-nubox/2,nubox/2]
        plt.imshow(data,origin='lower',interpolation='none',aspect='auto',extent=exten)
        plt.xlabel(r'time s')
        plt.ylabel(r'frequency kHz')

        rebin_dpiy,rebin_dpix=data.shape
    else:
        data=select_rebin(data,dpix,dpiy,xbins,ybins,ypositive=1)
        rebin_dpiy,rebin_dpix=data.shape
        # x: fs; y:tau
        exten = [-dpix/tbox,dpix/tbox,0,dpiy/nubox] 
        print exten
        figsize = (rebin_dpix/dpi,rebin_dpiy/dpi)
        fig=plt.figure(1,figsize=figsize,dpi=dpi)

        
        if spec_type == 'sec':
            data = np.log10(np.absolute(data)**2.)
        elif spec_type == 'conj':
            data = np.log10(np.absolute(np.imag(data)))
            name = 'cross_sec_i'+name.replace('_pol_LL.rebint','')
        elif spec_type == 'phase':
            data = np.sign(data)*abs(data)**(1./8)		

        if vmin!=0 and vmax!=0:
            plt.imshow(data,origin='lower',interpolation='none',aspect='auto',cmap= cmap,extent=exten,vmin=vmin,vmax=vmax)
        else:
            plt.imshow(data,origin='lower',interpolation='none',aspect='auto',cmap= cmap,extent=exten)
        plt.xlabel(r'differential doppler frequency $f_D$ mHz')
        plt.ylabel(r'time delay $\tau$ ms')
        plt.colorbar()

    #	plt.xlim(-60,60)
    #	plt.ylim(0,1.2)

    plt.title(name)
    #for interactive mode	
    if show == 'interact':	
        def onclick(event):
            print('xdata=%f, ydata=%f' %(event.xdata, event.ydata))
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
    elif show=='plot':
        plt.show()
    elif show=='save':
    #	plt.savefig(spec_type+name+'.pdf')
        figsize = (rebin_dpix/dpi,rebin_dpiy/dpi)
        pdf.savefig(1,figsize=figsize,dpi=dpi)
        print fig.get_size_inches()
        print 'figsize',figsize,'dpi',dpi
    #		plt.clf() # clear figure for next plot
    return 
#=============================================
#x is an array denote x axis
def plot_parabola(x,x0,y0,type='main',delaycut=0,delayup=0.5):
	import arclet as arc
	y=np.zeros(len(x))
	for i in range(len(x)):
		if type=='main':
			y[i]=arc.main_parabola_y(x[i],x0,y0)
		else:
			y[i]=arc.parabola(x[i],x0,y0)

	posi1 = y>delaycut
	posi2 = y<delayup
	x=x[posi1*posi2]
	y=y[posi1*posi2]
	plt.plot(x,y)
#	plt.xlim(-60,60)
#	plt.ylim(0,1)
#	pdf.savefig(1)
#	plt.show()
	return 
#================================================
def plot_noisebox(noisebox):
	xline=[noisebox[2],noisebox[3]]
	yline=[noisebox[0],noisebox[1]]
	plt.plot(xline,[yline[0],yline[0]])
	plt.plot(xline,[yline[1],yline[1]])
	plt.plot([xline[0],xline[0]],yline)
	plt.plot([xline[1],xline[1]],yline)
	return

#================================================
def plot_fd_arg(fd,arg,title='fd_arg',xerr=0,yerr=0):
    plt.xlabel(r'differential doppler frequency $f_D$ mHz')
    #	plt.ylabel(r'arg LR* rad')
    plt.ylabel(r'differential Faraday rotation')
    plt.title(title)
    if np.sum(xerr)==0 and np.sum(yerr)==0:
       plt.plot(fd,arg,'b.',markersize=1)
    else:
        plt.errorbar(fd, arg, xerr=xerr, yerr=yerr,fmt='.',markersize=1,ecolor='g', elinewidth=1,capsize=1)
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

def dy2conj(save=1,plot=1,inpath = '/mnt/scratch-lustre/simard/B0834_2012/Gate6/',inname = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'):
	figure_path='/home/dzli/ism/figures/'
	outname='gb057_baseline257_freq_00_pol_LL'
#output file
	outpath = '/mnt/raid-cita/dzli/gb057/Gate6/'

	dy=readdy(name=inname,path=inpath)
	ic=sp.conj_spec(dy)
	del dy
	if save==1:
		sp.save_conj(ic,outname+'.rebint',path=outpath)
	if plot==1:
		cmap='jet'
		figurename = outname+'.pdf'
		pdf = pgs.PdfPages(figure_path+figurename)
		sp.plot_spec(ic,spec_type='sec',name=outname, cmap = cmap,show= 'save',tbox=tbox,bins=4,pdf=pdf,dpix=500,dpiy=8100,dpi=250)
		pdf.close()

def get_noise(ic,box,tbox=5.7,nubox=8000.,lenx=1321,leny=16384):
	noisebox=[int(box[0]*nubox)+leny//2,int(box[1]*nubox)+leny//2,int(box[2]*tbox)+lenx//2,int(box[3]*tbox)+lenx//2]
	print noisebox
	icb=ic[noisebox[0]:noisebox[1],noisebox[2]:noisebox[3]]
	mean=np.mean(icb)
	noise=np.mean(np.absolute(icb-mean))
	print 'noise ',noise
	return noise
