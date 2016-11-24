task = 11
plotphase=0
mask_lowI = 6.7
vrange = 0.0
filenumber = 1

figure_path='/home/dzli/ism/figures/'
data_path='/home/dzli/ism/data/'
lenx=1321 #time
leny=16384 #freq
tbox=5.095424 # e3 s
nubox=8000. #kHz
# 0 plot secondary specs from dy_spec 
# 1 save conj_spec data
# 2 plot secondary specs from conj_spec
# 3 see LL minus RR sec_spec to check if there is displacement
# 4 calculate cross sec_spec: (sec_LL sec_RR*)
# 5 plot phase difference
# 6 check LL+RR and LL-RR
# 7 check LL+RR and LL-RR in log plot

# 8 average phase of 4 whole area to reduce statistical error
# 9 average Ei of 4 area before calc phase, with mask 6.5

#11 arclet deconvolve

# 1 auto goes to 2; 
# 4 auto goes to 5
#=============================================================

import spec as sp
import glob,os
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as pgs 
import numpy as np

os.chdir('/mnt/scratch-lustre/simard/B0834_2012/Gate0/')
[namer] = glob.glob('*freq_00_pol_RR.rebint')
[namel] = glob.glob('*freq_00_pol_LL.rebint')
#for task 0,1,2: names = namers
#go1
if task == 1:
	for i in xrange(0,filenumber):
		name = namer.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name = name)
		sp.conj_spec(data,save = 'y',name = name)

	for i in xrange(0,filenumber):
		name = namel.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name = name)
		sp.conj_spec(data,save = 'y',name = name)
	task=2
#go2
if task == 2:
	cmap = 'Dark2'
	cmap = 'Greys'
	namell = namel
	for j in xrange(0,1):
		figurename = 'sec_pol_LL_grey.pdf'
		pdf = pgs.PdfPages(figure_path+figurename)
		for i in xrange(0,filenumber):
			name = namell.replace('freq_00','freq_0'+str(i))
			ic = sp.readconj(name=name)
			sp.plot_spec(ic,pdf,spec_type='sec',name=name, cmap = cmap,show= 'insave')
			para=1
			if para == 1:	
				doppler =(-np.arange(lenx/2.))/tbox			
				doppler1,doppler2,delay1,delay2 = [-6.0,-6.0,0.116,0.088] 
				sp.plot_parabola(doppler,doppler1,delay1,pdf)
				sp.plot_parabola(doppler,doppler2,delay2,pdf)
				doppler1,doppler2,delay1,delay2 = [-11.3,-11.3,0.4,0.37] 
				sp.plot_parabola(doppler,doppler1,delay1,pdf)
				sp.plot_parabola(doppler,doppler2,delay2,pdf)
				doppler1,doppler2,delay1,delay2 = [-10.1,-10.1,0.3,0.27] 
				sp.plot_parabola(doppler,doppler1,delay1,pdf)
				sp.plot_parabola(doppler,doppler2,delay2,pdf)
				doppler1,doppler2,delay1,delay2 = [-8.2,-8.0,0.2,0.16] 
				sp.plot_parabola(doppler,doppler1,delay1,pdf)
				sp.plot_parabola(doppler,doppler2,delay2,pdf)
			#plt.show()
			pdf.savefig(1)
			plt.clf()
		pdf.close()
		figurename = figurename.replace('LL','RR')
		namell = namell.replace('LL','RR')

#go3
if task == 3:
	import check as ck
	pdf = pgs.PdfPages(figure_path+'sec_LLminusRR_nolog.pdf')
	for i in xrange(0,filenumber):
		namell = namel.replace('freq_00','freq_0'+str(i))
		namerr = namer.replace('freq_00','freq_0'+str(i))
		ck.check_lr(pdf,spec_type='sec',namel=namell,namer=namerr)
	pdf.close()
#go4
if task == 4:
	for i in xrange(0,filenumber):
		namerr = namer.replace('freq_00','freq_0'+str(i))
		namell = namel.replace('freq_00','freq_0'+str(i))
		savename = namell.replace('_pol_LL','')
		crosssec = sp.cross_sec(namell,namerr,save='y',outname = savename)
#go5
if task == 5:
	import plot
	name = ''
	icd = sp.readconj(name=name+namer)
	extent = [300,1000,0,leny/2]
	icd = icd[extent[2]:extent[3],extent[0]:extent[1]]
	filename = figure_path+name+'.png'
	plot.plot_png(icd,filename)

#go6
if task == 6:
	import phase as psi

	mask = ''
	if mask_lowI !=0:
		mask = '_mask'+str(mask_lowI)
	if vrange != 0:
		mask = mask + '_vrange'+str(vrange)
#plot phase LL minus RR
	pdf = pgs.PdfPages(figure_path+'phase_LLminusRR'+mask+'_'+str(filenumber)+'.pdf')
		
	# set a filter, only consider phasel+phaser(ideally = 0) < vrange
	phasel = psi.readphase0(name = namel,mask_lowI=mask_lowI)
	phaser = psi.readphase0(name = namer,mask_lowI=mask_lowI)
	check_ful = phasel - phaser
	del phasel,phaser	
	renorm = check_ful>np.pi
	check_ful[renorm]-=2.0*np.pi
	renorm = check_ful<-np.pi
	check_ful[renorm]+=2.0*np.pi
	del renorm

	for i in xrange(0,filenumber):
		check = check_ful
		filter_indices = abs(check) < vrange	
		print np.count_nonzero(filter_indices)
		check[np.logical_not(filter_indices)]=0.0
		del filter_indices

		name = namel.replace('pol_LL','vrange'+str(vrange))
		psi.plot_phase(check,pdf,name=name,vrange = vrange)
		del check
		vrange /=10
	pdf.close()
	del check_ful

#go7
if task == 7:
	sign = 1.0
	cmap = 'Dark2'
	import phase as psi
	mask = 'sign'+str(int(sign))+'_log_'+cmap
#plot phase LL minus RR
	pdf = pgs.PdfPages(figure_path+'phase_LLminusRR_'+mask+'_'+str(filenumber)+'.pdf')
		
	# set a filter, only consider phasel+phaser(ideally = 0) < vrange
	phasel = psi.readphase0(name = namel,mask_lowI=mask_lowI)
	phaser = psi.readphase0(name = namer,mask_lowI=mask_lowI)
	check = phasel - phaser
	del phasel,phaser	
	renorm = check>np.pi
	check[renorm]-=2.0*np.pi
	renorm = check<-np.pi
	check[renorm]+=2.0*np.pi
	del renorm

	for i in xrange(0,filenumber):
		filter_indices = sign*check >0	
		print np.count_nonzero(filter_indices)
		check[np.logical_not(filter_indices)]=0.0
		del filter_indices

		check = np.log10(abs(check))

		name = namel.replace('pol_LL',mask)
		psi.plot_phase(check,pdf,name=name,vrange = vrange,cmap = cmap)
	pdf.close()
	del check


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
if task == 10:
	wiener = 'y'
	from scipy.fftpack import fftshift,fft2,ifft2,ifftshift
	import arclet as arc
	cmap = 'Greys'
	namell = namel
	for j in xrange(0,2):
		for i in xrange(0,filenumber):
			delay = (np.arange(leny)-np.floor(leny/2.))/nubox
			doppler =(np.arange(lenx)-np.floor(lenx/2.))/tbox			
			doppler1,doppler2,delay1,delay2 = [-8.2,-8.0,0.2,0.16] 
			print 'get mask'
			mask = arc.get_mask(leny,lenx,doppler,delay,doppler1,doppler2,delay1,delay2)
			print 'non zeron mask: ',np.count_nonzero(mask)
			print 'deconvolve'
			name = namell.replace('freq_00','freq_0'+str(i))
			ic = sp.readconj(name=name)
			dyarc = ifft2(ifftshift(ic*mask))
			del ic,mask
			dy = sp.readdy(name = name)
			if wiener == 'y':
				icd = arc.wiener_deconvolution(dy,dyarc)
			else:
				ang = np.angle(dyarc)
				icd = fftshift(fft2(dy*np.exp(-1.0j*ang)))
		#		print 'average arclet angle', np.mean(ang),np.std(ang)
				del ang
#after deconvolution, move back to center
			delayi = np.argmin(abs(delay-np.average([delay1,delay2])))
			doppleri = np.argmin(abs(doppler-np.average([doppler1,doppler2])))
			icd = np.roll(np.roll(icd,doppleri - len(doppler)/2,1),delayi-len(delay)/2,0)
			del dy,dyarc
			sp.save_conj(icd,'wiener_recon_fd-8_'+name)
		namell = namell.replace('LL','RR')
	task+=1
#go11
if task == 11:
	name = 'recon_fd-8_' 
	name = ''
	icd = sp.readconj(name=name+namer)
	pdf = pgs.PdfPages(figure_path+name+'RR.pdf')
	sp.plot_spec(icd,pdf,spec_type='sec',name=name, tbox=tbox,nubox=nubox,show='save')
	pdf.close()

#go12
if task == 12:
	name = 'recon_fd-8_' 
	content = 'real'
	icl = sp.readconj(name=name+namel)
	icr = sp.readconj(name=name+namer)
	cc = icl*np.conjugate(icr)
#	phase = np.angle(cc) 
	pdf = pgs.PdfPages(figure_path+name+'ri.pdf')
	sp.plot_spec(np.real(cc),pdf,spec_type='phase',name=content+name, tbox=tbox,nubox=nubox,show='save')
	plt.clf()
	content = 'imag'
	sp.plot_spec(np.imag(cc),pdf,spec_type='phase',name=content+name, tbox=tbox,nubox=nubox,show='save')
	del cc
	pdf.close()
