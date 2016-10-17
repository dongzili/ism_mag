task = 4
plotphase=1
mask_lowI = 8.0

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
namers = glob.glob('*pol_RR.rebint')
namels = glob.glob('*pol_LL.rebint')
#for task 0,1,2: names = namers
if task == 0:
	pdf = pgs.PdfPages(figure_path+'sec_pol_RR.pdf')
	for name in names:
		data = sp.readdy(name=name)
	#sp.plot_spec(data,save='y')
		ic = sp.conj_spec(data)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

elif task == 1:
	for name in names:
		data = sp.readdy(name = name)
		sp.conj_spec(data,save = 'y',name = name)
	task=2

if task == 2:
	pdf = pgs.PdfPages(figure_path+'sec_pol_RR.pdf')
	for name in names:
		ic = sp.readconj(name=name)
		sp.plot_spec(ic,pdf,spec_type='sec',name=name)
	pdf.close()

if task == 3:
	import check as ck
	pdf = pgs.PdfPages(figure_path+'sec_LLminusRR.pdf')
	for i in xrange(0,len(namels)):
		ck.check_lr(pdf,spec_type='sec',namel=namels[i],namer=namers[i])
	pdf.close()

elif task == 4:
	import phase as psi 
	for i in xrange(0,len(namels)):
		phasel = psi.phase0(name = namels[i],save ='y',mask_lowI = mask_lowI)
		phaser = psi.phase0(name = namers[i],save ='y',mask_lowI = mask_lowI)
		phase_diff = psi.phase0_LminusR(phasel,phaser,save='y',name = namels[i],mask_lowI=mask_lowI)
		del phasel,phaser,phase_diff
	task = 5

if task == 5:
	import phase as psi

	if mask_lowI !=0:
		mask = '_mask'+str(mask_lowI)+'_'
#plot phase LL minus RR
	pdf = pgs.PdfPages(figure_path+'phase_LLminusRR'+mask+'.pdf')
	for i in xrange(0,len(namels)):
		phase = psi.readphase0(name = namels[i],phase_type = 'diff',mask_lowI=mask_lowI)
		psi.plot_phase(phase,pdf,name=namels[i])
		del phase
	pdf.close()

	if plotphase == 1:
	#plot phase LL and phase RR
		pdf = pgs.PdfPages(figure_path+'phase_LL'+mask+'.pdf')
		for i in xrange(0,len(namels)):
			phase = psi.readphase0(name = namels[i],mask_lowI=mask_lowI)
			psi.plot_phase(phase,pdf,name=namels[i])
			del phase
		pdf.close()
		
		pdf = pgs.PdfPages(figure_path+'phase_RR'+mask+'.pdf')
		for i in xrange(0,len(namels)):
			phase = psi.readphase0(name = namers[i],mask_lowI=mask_lowI)
			psi.plot_phase(phase,pdf,name=namers[i])
			del phase
		pdf.close()

