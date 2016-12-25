task = 19
plotphase=0
mask_lowI = 6.7
vrange = 0.0
filenumber = 1

lenx=1321 #time
leny=16384 #freq
tbox=6.729 # e^3 s (from .time)
nubox=8000. #kHz
f0=322.5
fcenters=[314.5,322.5,330.5,338.5]
inpath = '/mnt/scratch-lustre/simard/B0834_2012/Gate0/'
inname = 'gb057_1.input_baseline257_freq_00_pol_LL.rebint'
#output demo
figure_path='/home/dzli/ism/figures/'
data_path='/home/dzli/ism/data/'
outname='gb057_baseline257_freq_00_pol_'
#output file
outpath = '/mnt/raid-cita/dzli/gb057/'
slowft_path = '/mnt/raid-cita/dzli/gb057/Gate0/'
slowft_name = 'gb057_baseline257_freq_00_pol_LL_f0_322.5'
multigate_path = '/mnt/raid-cita/dzli/gb057/multigate/'
multigate_name_old = 'LL_gate012_gb057_baseline257_freq_00.rebint'
multigate_name='gb057_baseline257_gate012_snvary_'
pol=['LL','RR']
freq=['00','01','02','03']
freq_count,pol_count=4,2

# task0 add gates
# task1 plot secondary specs of 4 freq
# task2 plot the arclet cut
# task3 arclet deconvolution
# 4 calculate cross sec_spec: (sec_LL sec_RR*)
# 5 plot phase difference
# 6 check LL+RR and LL-RR
# 7 check LL+RR and LL-RR in log plot

# 8 average phase of 4 whole area to reduce statistical error
# 9 average Ei of 4 area before calc phase, with mask 6.5

#10 arclet deconvolve
#task16 calculate cc, cut the main parabola
#task17 fit fd-psi

#=============================================================

import arclet as arc
import spec as sp
import glob,os
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as pgs 
from scipy.fftpack import ifft2,ifftshift,fft,fftshift
import numpy as np

os.chdir('/mnt/scratch-lustre/simard/B0834_2012/Gate0/')
[namer] = glob.glob('*freq_00_pol_RR.rebint')
[namel] = glob.glob('*freq_00_pol_LL.rebint')
#for task 0,1,2: names = namers
#add gates012
#go0
if task == 0:
	import add_gates as gate
	gate.sum_dy()
	#task+=1

#plot added sec spec
#go1
if task==1:
    cmap,plot_small_arc='jet',1
    xbins,ybins,vmin,vmax,dpix,dpiy=(1,10,6,10.5,330,6400)
    hat='ybins'+str(ybins)+'_'
    if plot_small_arc==1:
        hat='small_arc_'+hat
    for num_pol in xrange(0,pol_count):
        figurename = multigate_name+'pol_'+pol[num_pol]+'.pdf'
        pdf = pgs.PdfPages(figure_path+hat+figurename)
        for num_freq in xrange(0,freq_count):
    #        sn = sn_dy(num_freq=num_freq)
            name=multigate_name+'freq_'+freq[num_freq]+'_pol_'+pol[num_pol]
            ic=sp.readconj(name=name+'.rebint',path=multigate_path)
            sp.plot_spec(ic,spec_type='sec',name=name, cmap = cmap,show= 'continue',tbox=tbox,xbins=xbins,ybins=ybins,pdf=pdf,dpix=dpix,dpiy=dpiy,vmin=vmin,vmax=vmax)

            if plot_small_arc == 1:
                print 'plot small arc'
                doppler =(-np.arange(dpix))/tbox			
                f = fcenters[num_freq]
                doppler1,doppler2,delay1,delay2 = [-16.8*f/f0,-16.8*f/f0,0.1463,0.1363] 
                sp.plot_parabola(doppler,doppler1,delay1,type='sec')
                sp.plot_parabola(doppler,doppler2,delay2,type='sec')
        #plt.show()
            pdf.savefig(1)
            plt.clf()
       
        pdf.close()


#arclet deconvolution
#go3
if task == 3:
    cmap='jet'
    pol_count,freq_count=2,1
    saveic,savedy=1,1
    xbins,ybins,vmin,vmax,dpi,dpix,dpiy=(1,8,6,10.5,50,400,6400)
    delaycut,doppler_keep,dopplercut,delay_keep=0.006,20.,3./tbox,0.02
    hat='process_ybins%.0f'%ybins+'_'#+str(delaycut)+'_'+str(dopplercut)
    for num_pol in xrange(0,pol_count):
        hat='process_ybins%.0f'%ybins+'_delaycut%.0f'%delaycut+'dopplercut%.0f'%dopplercut
        figurename = multigate_name+'pol_'+pol[num_pol]+'.pdf'
        pdf = pgs.PdfPages(figure_path+hat+figurename)
        for num_freq in xrange(0,freq_count):
    #        sn = sn_dy(num_freq=num_freq)
            name=multigate_name+'freq_'+freq[num_freq]+'_pol_'+pol[num_pol]
            ic=sp.readconj(name=name+'.rebint',path=multigate_path)
            if delaycut!=0:
                exclude=int(delaycut*nubox)
                ic[leny//2-exclude:leny//2+exclude+1,0:lenx//2-np.ceil(doppler_keep*tbox)+1]=0
                ic[leny//2-exclude:leny//2+exclude+1,lenx//2+np.ceil(doppler_keep*tbox):lenx+1]=0

            if dopplercut!=0:
                exclude=int(dopplercut*tbox)
                ic[0:leny//2-np.ceil(delay_keep*nubox)+1,lenx//2-exclude:lenx//2+exclude+1]=0
                ic[leny//2+np.ceil(delay_keep*nubox):leny,lenx//2-exclude:lenx//2+exclude+1]=0

            sp.plot_spec(ic,spec_type='sec',name=outname, cmap = cmap,show= 'save',tbox=tbox,xbins=xbins,ybins=ybins,pdf=pdf,dpix=dpix,dpiy=dpiy,dpi=dpi,vmin=vmin,vmax=vmax)
            plt.clf()

        print 'original shape',ic.shape
        ic=sp.select_rebin(ic,dpix,dpiy,xbins,ybins)
        print 'original shape',ic.shape
        hat='_ybins%.0f'%ybins+'_'+'array_f%.0f'%ic.shape[0]+'t%.0f'%ic.shape[1]+'_taufcutted'
        name=name+hat+'.rebint'
        if saveic==1:
            print np.sum(ic.imag)
            sp.save_conj(ic,name,path=multigate_path)
        if savedy==1:
            from scipy.fftpack import ifft2,ifftshift
            dy=ifft2(ifftshift(ic))
            print 'r,i dy',np.sum(dy.real),np.sum(dy.imag)
            hat='_ybins%.0f'%ybins+'_'+'array_f%.0f'%dy.shape[0]+'t%.0f'%dy.shape[1]+'_taufcutted'
            sp.save_dy(dy.real,name,path=multigate_path)
        check=1
        if check==1:
            from scipy.fftpack import fft2,fftshift
            ic=fftshift(fft2(dy.real))
            sp.plot_spec(ic,spec_type='sec',name='real dy', cmap = cmap,show= 'save',tbox=tbox,ybins=ybins,pdf=pdf,dpix=dpix,dpiy=dpiy,dpi=dpi,vmin=vmin,vmax=vmax)
            plt.clf()

        pdf.close()

#go4
if task == 4:
	for i in xrange(0,filenumber):
		namerr = namer.replace('freq_00','freq_0'+str(i))
		namell = namel.replace('freq_00','freq_0'+str(i))
		savename = namell.replace('_pol_LL','')
		crosssec = sp.cross_sec(namell,namerr,save='y',outname = savename)

#go8
if task == 8:
	import phase as psi

	mask = ''
	if mask_lowI !=0:
		mask = '_mask'+str(mask_lowI)
	if vrange != 0:
		mask = mask + '_vrange'+str(vrange)
# select box
	box_lu = [400,600,9000,11000] #lrbt
	print 'box: ',box_lu

#plot phase LL minus RR
	pdf = pgs.PdfPages(figure_path+'phase_LLminusRR'+mask+'_average'+str(box_lu)+str(filenumber)+'.pdf')
	for i in xrange(0,filenumber):
		phase = psi.readphase0(name = namel,phase_type = 'diff',mask_lowI=mask_lowI)
		lenx = len(phase[0])
		leny = len(phase)
		
		mean = np.arange(4)
		mean[0] = sp.mean_spec(phase,box_lu,cover='yes')
	
		box_ru = [lenx-box_lu[1],lenx-box_lu[0],box_lu[2],box_lu[3]]
		print 'box: ',box_ru
		mean[1] = sp.mean_spec(phase,box_ru,cover='yes')

		box_ld = [box_lu[0],box_lu[1],leny-box_lu[3],leny-box_lu[2]]
		print 'box: ',box_ld
		mean[2] = sp.mean_spec(phase,box_ld)

		box_rd = [box_ru[0],box_ru[0],leny-box_ru[3],leny-box_ru[2]]
		print 'box: ',box_rd
		mean[3] = sp.mean_spec(phase,box_ru)

		name = namel.replace('freq_00_pol_LL','freq_0'+str(i)+'LLminusRR')
		name = namel.replace('pol_LL','average_phase'+str(box_lu))
		name = 'freq_00'+mask+repr(mean)
		psi.plot_phase(phase,pdf,name=name,vrange = vrange)
		del phase
	pdf.close()


#go9
if task == 9:
	import phase as psi

# select box
	boxnumber = 2
	box = np.empty([boxnumber,4],dtype=int)
	box[0,0:4] = [300,625,8400,11000] #lrbt
	box[1,0:4] = [lenx-box[0][1],lenx-box[0][0],box[0][2],box[0][3]]

#plot phase LL minus RR
	dat = open (data_path+'mean_phase.dat','a')
	dat.write('#################\n')
	for i in xrange(0,filenumber):
		namell = namel.replace('freq_00','freq_0'+str(i))
		dat.write('#==============\n'+namell+'\n')
		for j in xrange(0,boxnumber):
			phasel,num_points= psi.mean_phase(namell,mask_lowI=mask_lowI,box=box[j])
			dat.write(str(box[j])+'points averaged:'+str(num_points)+' ')
			dat.write('phase diff: %.5f \n' %phasel)
	dat.close()
#go10
delaycut,doppler_keep,dopplercut,delay_keep=0.006,0.,3./tbox,0.0
arcdelaycut=0.006
hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
if task == 10:
    wiener='y'
    from scipy.fftpack import fftshift,fft2,ifft2,ifftshift
    import arclet as arc
    path=multigate_path
    for num_freq in xrange(0,freq_count):
        for num_pol in xrange(0,pol_count):
            f = fcenters[num_freq]
            delay = (np.arange(leny)-np.floor(leny/2.))/nubox
            doppler =(np.arange(lenx)-np.floor(lenx/2.))/tbox			
            doppler1,doppler2,delay1,delay2 = [-16.8*f/f0,-16.8*f/f0,0.1463,0.1363] 
            print 'get mask'
            mask = arc.get_mask(leny,lenx,doppler,delay,doppler1,doppler2,delay1,delay2,delaycut=arcdelaycut)
            print 'non zeron mask: ',np.count_nonzero(mask)
            print 'deconvolve'
            name=multigate_name+'freq_'+freq[num_freq]+'_pol_'+pol[num_pol]+'.rebint'
            ic = sp.readconj(name=name,path=path)
            dyarc = ifft2(ifftshift(ic*mask))
            del mask
            if delaycut*dopplercut!=0:
                if delaycut!=0:
                    exclude=int(delaycut*nubox)
                    ic[leny//2-exclude:leny//2+exclude+1,0:lenx//2-np.ceil(doppler_keep*tbox)+1]=0
                    ic[leny//2-exclude:leny//2+exclude+1,lenx//2+np.ceil(doppler_keep*tbox):lenx+1]=0

                if dopplercut!=0:
                    exclude=int(dopplercut*tbox)
                    ic[0:leny//2-np.ceil(delay_keep*nubox)+1,lenx//2-exclude:lenx//2+exclude+1]=0
                    ic[leny//2+np.ceil(delay_keep*nubox):leny,lenx//2-exclude:lenx//2+exclude+1]=0
                from scipy.fftpack import ifft2,ifftshift
                dy=ifft2(ifftshift(ic))

            else:
                del ic
                dy = sp.readdy(name = name,path=path)
    #			dy = np.fromfile(slowft_path+slowft_name+'_nopad.rebint',dtype=complex).reshape((leny,lenx))

            if wiener == 'y':
                icd = arc.wiener_deconvolution(dy,dyarc)
            else:
                ang = np.angle(dyarc)
                icd = fftshift(fft2(dy*np.conjugate(dyarc)/abs(dyarc)))
                print 'average arclet angle', np.mean(ang),np.std(ang)
                del ang
    #after deconvolution, move back to center
            delayi = np.argmin(abs(delay-np.average([delay1,delay2])))
            doppleri = np.argmin(abs(doppler-np.average([doppler1,doppler2])))
            icd = np.roll(np.roll(icd,doppleri - len(doppler)/2,1),delayi-len(delay)/2,0)
            del dy,dyarc
            sp.save_conj(icd,hat+name,path=path)
    task+=1
#go11
if task == 11:
    hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
    xbins,ybins,vmin,vmax,dpi,cmap=1,10,6.5,11,50,'jet'
    for num_pol in xrange(0,pol_count):
        figurename = multigate_name+'pol_'+pol[num_pol]
        pdf = pgs.PdfPages(figure_path+hat+figurename+'.pdf')
        for num_freq in xrange(0,freq_count):
            path=multigate_path
            name=multigate_name+'freq_'+freq[num_freq]+'_pol_'+pol[num_pol]+'.rebint'
            icd = sp.readconj(name=hat+name,path=path)
            sp.plot_spec(icd,pdf,spec_type='sec',name=figurename, tbox=tbox,nubox=nubox,show='save',ybins=10,cmap = cmap,dpix=350,dpiy=5000,vmin=vmin,vmax=vmax)
            plt.clf()
        pdf.close()
#    task+=2

#go12
if task == 12:
	cmap = 'Greys'
	name = 'wiener_recon_fd-17_' 
	icl = sp.readconj(name=name+namel)
	icr = sp.readconj(name=name+namer)
	name = 'sec_grey_v5-20'+name 
	dpix,dpiy=[300,4000]
	cc = icl*np.conjugate(icr)
#	phase = np.angle(cc) 
	mask1 = np.real(cc)>0
	mask2 = np.imag(cc)>0
	ybins=20
	pdf = pgs.PdfPages(figure_path+'cc-1'+name+'.pdf')
	sp.plot_spec(np.real(cc)*(mask1-1),pdf,spec_type='sec',name='cc-1_real_'+name, tbox=tbox,nubox=nubox,show='save',cmap=cmap,ybins=ybins,dpix=dpix,dpiy=dpiy)
	plt.clf()
	sp.plot_spec(np.imag(cc)*(mask2-1),pdf,spec_type='sec',name='cc-1_ima'+name, tbox=tbox,nubox=nubox,show='save',cmap=cmap,ybins=ybins,dpix=dpix,dpiy=dpiy)
	plt.clf()
	pdf.close()
	pdf = pgs.PdfPages(figure_path+'cc+1'+name+'.pdf')
	sp.plot_spec(np.real(cc)*mask1,pdf,spec_type='sec',name='cc+1_real'+name, tbox=tbox,nubox=nubox,show='save',cmap=cmap,dpix=dpix,dpiy=dpiy,ybins=ybins)
	plt.clf()
	sp.plot_spec(np.imag(cc)*mask2,pdf,spec_type='sec',name='cc+1_ima'+name, tbox=tbox,nubox=nubox,show='save',cmap=cmap,dpix=dpix,dpiy=dpiy,ybins=ybins)
	del cc
	pdf.close()

#go13
if task == 13:
	cmap = 'jet'
	path = multigate_path
	namell,namerr= [multigate_name,multigate_name.replace('LL','RR')]
	icl = sp.readconj(name=hat+namell,path=path)
	icr = sp.readconj(name=hat+namerr,path=path)
	dpix,dpiy=[300,4000]
	cc = icl*np.conjugate(icr)
#	phase = np.angle(cc) 
	dpi,xbins,ybins=[50,2,20]
	pdf = pgs.PdfPages(figure_path+'cc_'+hat+str(ybins)+'.pdf')
	sp.plot_spec(np.real(cc),pdf,spec_type='phase',name='cc_real_'+hat, tbox=tbox,nubox=nubox,show='save',cmap=cmap,ybins=ybins,dpix=dpix,dpiy=dpiy,dpi=dpi,xbins=xbins)
	plt.clf()
	sp.plot_spec(np.imag(cc),pdf,spec_type='phase',name='cc_ima'+hat, tbox=tbox,nubox=nubox,show='save',cmap=cmap,ybins=ybins,dpix=dpix,dpiy=dpiy,dpi=dpi,xbins=xbins)
	plt.clf()
	pdf.close()
#go14
if task == 14:
#	import slowft as ft
#	ft.main()
	if 1==0:
		dy1 = np.fromfile(slowft_path+slowft_name+'_nopad.rebint',dtype=complex).reshape((leny,lenx))
		print dy1.shape
		pdf = pgs.PdfPages(figure_path+'dy_'+slowft_name+'LL.pdf')
		sp.plot_spec(np.real(dy1[0:leny-1,:]),pdf=pdf,spec_type='dy',name='slowft', tbox=tbox,nubox=nubox,show='plot',ybins=4)
		pdf.close()
	else:
		task+=1
#go15
if task == 15:
	path = multigate_path
	name = 'gb057_baseline257_freq_00_f0_322.5'
	if 1 == 0:
		from scipy.fftpack import fft2, fftshift
		dy1 = np.fromfile(path+name+'_nopad.rebint',dtype=complex).reshape((leny,lenx))
		ic1 = fftshift(fft2(dy1))
		del dy1
		sp.save_conj(ic1,name,path=path)
	else:
		ic1 = sp.readconj(path=path,name=name,num_t=lenx,num_nu=leny)
	pdf = pgs.PdfPages(figure_path+'totslowft.pdf')
	f=314.5
	f=f0
	sp.plot_spec(ic1,pdf,spec_type='sec',name='', tbox=tbox/f*f0,nubox=nubox,show='save',ybins=6,dpix=500,dpiy=8100,dpi=100,cmap='jet')
	pdf.close()


#go16
if task == 16:
    delaycut,dopplercut,dpix,dpiy,dpi=0.006,3./tbox,400,6400,50
    plot,save=[1,1]
    delayup=0.8
    noisebox=[delaycut,1.0,-80,-40]
    xbins,ybins,vmin,vmax,dpi,cmap=1,10,13,22,50,'jet'
    hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
    inname = hat+multigate_name
    figurename='cc_main_para_'+inname
    path = multigate_path
    if plot==1:
        pdf = pgs.PdfPages(figure_path+figurename+'.pdf')
    for num_freq in xrange(0,freq_count):
        name=inname+'freq_'+freq[num_freq]+'_pol_'+pol[0]+'.rebint'
        outname='cc_main_para_'+inname+'freq_'+freq[num_freq]
        icl = sp.readconj(name=name,path=path)
        icr = sp.readconj(name=name.replace('LL','RR'),path=path)
        ic=icl*np.conjugate(icr)
        del icl,icr
        delay = (np.arange(leny)-leny//2)/nubox
        doppler =(np.arange(lenx)-lenx//2)/tbox			
        f = fcenters[num_freq]
        doppler1,doppler2,delay1,delay2 = [-16.8*f/f0,-16.8*f/f0,0.18,0.13] 
        if plot==1:
            sp.plot_spec(ic,pdf,spec_type='sec',cmap = cmap,show= 'continue',tbox=tbox,ybins=ybins,dpix=dpix,dpiy=dpiy,dpi=dpi,vmin=vmin,vmax=vmax,name=outname)
            sp.plot_parabola(doppler,doppler1,delay1,delaycut=delaycut)
            sp.plot_parabola(doppler,doppler2,delay2,delaycut=delaycut)
            #plt.show()
            #plot noise box
            pdf.savefig(1)
            plt.clf()

        if save==1:
            mask=arc.get_main_para(leny,lenx,doppler,delay,doppler1,doppler2,delay1,delay2,delaycut=delaycut,delayup=delayup)	
            print 'non_zero_mask: ',np.count_nonzero(mask)
            sp.save_conj(ic*mask,outname+'.rebint',path=path)

        noise=sp.get_noise(ic,noisebox,tbox=tbox,nubox=nubox,lenx=lenx,leny=leny)

    if plot==1:
        pdf.close()

#go17
if task==17:
    delaycut,dopplercut,dpix,dpiy,dpi=0.006,3./tbox,400,6400,50
    delaylow,sn=0.1,2
    noise=2e7
    xbins,ybins,vmin,vmax,dpi,cmap=1,4,13,22,50,'jet'
    hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
    inname = 'cc_main_para_'+hat+multigate_name
    figurename='fd-psi_sn%.0f_ybins%.0f'%(sn,ybins)+inname
    path = multigate_path
    plot=1
    if plot==1:
        pdf = pgs.PdfPages(figure_path+figurename+'.pdf')

    dat = open (data_path+'ratio_tau%.1e.dat'%delaycut,'a')
    dat.write(inname+'\n')
    dat.write('sn %.0f'%sn+' ybins %.0f'%ybins+' xbins %.0f'%xbins)
    dat.write(' delaylow%.1e'%delaylow)
    for num_freq in xrange(0,freq_count):
        name=inname+'freq_'+freq[num_freq]+'.rebint'
        cc = sp.readconj(name=name,path=path)
        cc=sp.select_rebin(cc,dpix,dpiy,xbins,ybins,ypositive=1)
        if delaylow!=0:
            exclude=int(delaycut*nubox)
            cc=cc[exclude:cc.shape[0]+1,:]
        if sn!=0:
            select=cc.real>noise*sn
            print 'num of points selected',np.count_nonzero(select)
        #find x index
            
        (drop,select_fd)=np.where(select)
        select_fd=(select_fd-cc.shape[1]//2*xbins)/tbox
        del drop
        #find arg
        arg=np.angle(cc[select])
        sigma = abs(noise/cc[select].real)
    #		arg=cc[select].imag
    #===========================
    #plot x,y
        plotri=0
        if plotri==1:
            pdf2 = pgs.PdfPages(figure_path+'ri_'+inname+'.pdf')
            sp.plot_spec(np.real(cc),pdf2,spec_type='phase',name='cc_real_'+hat, show='save',cmap=cmap,dpiy=dpiy,dpix=dpix,dpi=dpi,tbox=tbox/xbins,nubox=nubox/ybins)
            plt.clf()
            sp.plot_spec(np.imag(cc),pdf2,spec_type='phase',name='cc_ima'+hat,show='save',cmap=cmap,dpix=dpix,dpiy=dpiy,dpi=dpi,tbox=tbox/xbins,nubox=nubox/ybins)
            plt.clf()
            pdf2.close()


        fit=1
        if fit==1:
            from scipy.optimize import curve_fit
            def func(x,a,b):
                return a*x+b
            popt,pcov=curve_fit(func,select_fd,arg,sigma=sigma,absolute_sigma=1)
            minfd,maxfd=[min(select_fd),max(select_fd)]
            mintau,maxtau=func(minfd,popt[0],popt[1]),func(maxfd,popt[0],popt[1])
            print 'para',popt,'error',pcov, 'maxfd %.1e maxtau %.1e'%(maxfd,maxtau)
    #			ratio=sum(cc[select].imag)/sum(cc[select].real)

        if plot==1:
            title='sn'+str(sn)+'_ybins'+str(ybins)+'_ratio: %.1e'%popt[0]
            sp.plot_fd_arg(select_fd,arg,title=title)
            pdf.savefig(1)
            plt.clf()
            sp.plot_fd_arg(select_fd,arg,title=title,yerr=sigma,xerr=xbins/tbox)
            pdf.savefig(1)
            plt.clf()

        dat.write('\n freq'+freq[num_freq]) 
        dat.write('\n %.1e  x+  %.1e'%(popt[0],popt[1]))
        dat.write('\n pcov %.1e %.1e %.1e %.1e'%(pcov[0,0],pcov[0,1],pcov[1,0],pcov[1,1]))
        dat.write('\n maxfd %.1e maxtau %.1e'%(maxfd,maxtau))
        dat.write('\n')
    
    dat.write('ybins%.0f \n \n'%ybins)
    dat.close()
    if plot==1:
        pdf.close()

#see phase fourier
#go18
if task==18:
    delaycut,dopplercut,dpix,dpiy,dpi=0.006,3./tbox,lenx//2,leny//2,100
    xbin,ybins=1,10
    vmin,vmax=0,0
    plot,save=[1,1]
    freq_count=1
    hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
    inname = hat+multigate_name
    figurename='cc_phase_fourier_'+inname
    path = multigate_path
    if plot==1:
        pdf = pgs.PdfPages(figure_path+figurename+'.pdf')
    for num_freq in xrange(0,freq_count):
        name=inname+'freq_'+freq[num_freq]+'_pol_'+pol[0]+'.rebint'
        outname='cc_phase_fourier_'+inname+'freq_'+freq[num_freq]
        icl = sp.readconj(name=name,path=path)
        icr = sp.readconj(name=name.replace('LL','RR'),path=path)
        ic=icl*np.conjugate(icr)
        ic=np.angle(ic)
        dy=ifft2(ifftshift(ic))
        del ic
        if plot==1:
            sp.plot_spec(np.absolute(dy)**1./16,pdf,spec_type='dy',show= 'continue',tbox=tbox,ybins=ybins,dpix=dpix,dpiy=dpiy,dpi=dpi,vmin=vmin,vmax=vmax,name=outname)
            pdf.savefig(1)
            plt.clf()
        if save==1:
            sp.save_dy(np.absolute(dy)**2,'pow_'+outname+'.rebint',path=path)

    if plot==1:
        pdf.close()

#go19
#plot noise fourier
if task==19:
    delaycut,dopplercut,dpix,dpiy,dpi=0.006,3./tbox,400,6400,50
    delaylow,sn=0.1,2
    noise=2e7
    xbins,ybins,vmin,vmax,dpi,cmap=1,1,13,22,50,'jet'
    hat='exclude_tau'+str(delaycut)+'_fd'+str(dopplercut)+'_wiener_recon_fd-17_'
    inname = 'cc_main_para_'+hat+multigate_name
    figurename='psi_div_fd_fourier_sn%.0f_ybins%.0f'%(sn,ybins)+inname
    path = multigate_path
    plot=1
    if plot==1:
        pdf = pgs.PdfPages(figure_path+figurename+'.pdf')

    for num_freq in xrange(0,freq_count):
        name=inname+'freq_'+freq[num_freq]+'.rebint'
        cc = sp.readconj(name=name,path=path)
        cc=sp.select_rebin(cc,dpix,dpiy,xbins,ybins,ypositive=1)
        if delaylow!=0:
            exclude=int(delaycut*nubox)
            cc=cc[exclude:cc.shape[0]+1,:]

        if sn!=0:
            select=cc.real>noise*sn
            cc[select]=0
            print 'num of points selected',np.count_nonzero(select)
        #find x index
            
        count=cc.real!=0
        print count.shape
        lenx=cc.shape[1]
        cc=sum(cc)
        count=sum(count)
        #print cc.shape()#,count.shape()
        select,=np.where(cc.real)
        #print select.shape()
        cc[select]/=count[select]
        select_fd=(select-lenx//2)*xbins/tbox
        #find arg
        arg=np.angle(cc[select])/select_fd
        dy=fftshift(fft(arg))
        del arg
        if plot==1:
            plt.plot(np.arange(len(dy)),np.absolute(dy))
            pdf.savefig(1)
            plt.clf()
    if plot==1:
        pdf.close()
