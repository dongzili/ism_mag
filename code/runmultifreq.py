task = 5
plotphase=0
mask_lowI = 0
vrange = 1e-3
filenumber = 3

figure_path='/home/dzli/ism/figures/'
# 0 plot secondary specs from dy_spec 
# 1 save conj_spec data
# 2 plot secondary specs from conj_spec
# 3 see LL minus RR sec_spec to check if there is displacement
# 4 calculate phase0 and phase0_LL-RR
# 5 plot phase difference
	#plotphase = 1 also plot phase RR, phase LL


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
if task == 0:
	pdf = pgs.PdfPages(figure_path+'sec_pol_RR.pdf')
	for i in xrange(0,filenumber):
		name = namer.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name=name)
	#sp.plot_spec(data,save='y')
		ic = sp.conj_spec(data)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

	pdf = pgs.PdfPages(figure_path+'sec_pol_LL.pdf')
	for i in xrange(0,filenumber):
		name = namel.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name=name)
	#sp.plot_spec(data,save='y')
		ic = sp.conj_spec(data)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

elif task == 1:
	for i in xrange(0,filenumber):
		name = namer.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name = name)
		sp.conj_spec(data,save = 'y',name = name)

	for i in xrange(0,filenumber):
		name = namer.replace('freq_00','freq_0'+str(i))
		data = sp.readdy(name = name)
		sp.conj_spec(data,save = 'y',name = name)
	task=2

if task == 2:
	pdf = pgs.PdfPages(figure_path+'sec_pol_RR.pdf')
	for i in xrange(0,filenumber):
		name = namel.replace('freq_00','freq_0'+str(i))
		ic = sp.readconj(name=name)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

	pdf = pgs.PdfPages(figure_path+'sec_pol_LL.pdf')
	for i in xrange(0,filenumber):
		name = namel.replace('freq_00','freq_0'+str(i))
		ic = sp.readconj(name=name)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

if task == 3:
	import check as ck
	pdf = pgs.PdfPages(figure_path+'sec_LLminusRR.pdf')
	for i in xrange(0,filenumber):
		namell = namel.replace('freq_00','freq_0'+str(i))
		namerr = namer.replace('freq_00','freq_0'+str(i))
		ck.check_lr(pdf,spec_type='sec',namel=namell,namer=namerr)
	pdf.close()

elif task == 4:
	import phase as psi 
	for i in xrange(0,filenumber):
		namell = namel.replace('freq_00','freq_0'+str(i))
		namerr = namer.replace('freq_00','freq_0'+str(i))
		phasel = psi.phase0(name = namell,save ='y',mask_lowI = mask_lowI)
		phaser = psi.phase0(name = namerr,save ='y',mask_lowI = mask_lowI)
		phase_diff = psi.phase0_LminusR(phasel,phaser,save='y',name = namell,mask_lowI=mask_lowI)
		del phasel,phaser,phase_diff
	task = 5

if task == 5:
	import phase as psi

	mask = ''
	if mask_lowI !=0:
		mask = '_mask'+str(mask_lowI)
	if vrange != 0:
		mask = mask + '_vrange'+str(vrange)
#plot phase LL minus RR
	pdf = pgs.PdfPages(figure_path+'phase_LLminusRR'+mask+'_'+str(filenumber)+'.pdf')
	for i in xrange(0,filenumber):
	#	namell = namel.replace('freq_00','freq_0'+str(i))
		namell = namel
		vrange *=10
		phase = psi.readphase0(name = namell,phase_type = 'diff',mask_lowI=mask_lowI)
		name = namel.replace('freq_00_pol_LL','freq_0'+str(i)+'LLminusRR')
		name = namel.replace('pol_LL','LLminusRR')
		psi.plot_phase(phase,pdf,name=name,vrange = vrange)
		del phase
	pdf.close()

	if plotphase == 1:
	#plot phase LL and phase RR
		pdf = pgs.PdfPages(figure_path+'phase_LL'+mask+'_'+str(filenumber)+'.pdf')
		for i in xrange(0,filenumber):
			namell = namel.replace('freq_00','freq_0'+str(i))
			phase = psi.readphase0(name = namell,mask_lowI=mask_lowI)
			psi.plot_phase(phase,pdf,name=namell,vrange=vrange)
			del phase
		pdf.close()
		
		pdf = pgs.PdfPages(figure_path+'phase_RR'+mask+'_'+str(filenumber)+'.pdf')
		for i in xrange(0,filenumber):
			namerr = namer.replace('freq_00','freq_0'+str(i))
			phase = psi.readphase0(name = namerr,mask_lowI=mask_lowI)
			psi.plot_phase(phase,pdf,name=namerr,vrange = vrange)
			del phase
		pdf.close()

