import multispec as sp
import glob,os
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as pgs 
import numpy as np

os.chdir('/mnt/scratch-lustre/simard/B0834_2012/Gate0/')
names = glob.glob('*pol_LL.rebint')
pdf = pgs.PdfPages('/home/dzli/ism/sec_pol_LL.pdf')
for name in names:
	data = sp.read(name=name)
#sp.plot_spec(data,save='y')
	ic = sp.sec_spec(data)
	sp.plot_spec(ic,pdf,spec_type='sec',name=name)

pdf.close()
