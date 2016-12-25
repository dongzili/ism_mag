import numpy as np
from scipy.fftpack import fftshift, fft2
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pgs 
#for write output
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')#set effective value

import spec as sp
import arclet as arc
figure_path='/home/dzli/ism/figures/'

leny=16384 #freq
tbox=6.729 # e^3 s (from .time)
nubox=8000. #kHz
f0=322.5


#f0: the average frequency you scaled
def gate_weight( baseline=257,freq='00',pol='LL',Gate=0,pdf=0):
    print 'Band selected: ' + freq
    inpath = '/mnt/scratch-lustre/simard/B0834_2012/Gate'+str(Gate)+'/'
    infile = 'gb057_1.input_baseline'+str(baseline)+'_freq_'+str(freq)+'_pol_'+pol+'.rebint'
    outpath = '/mnt/raid-cita/dzli/gb057/multigate/'
    outname='Gate'+str(Gate)+'freq_'+str(freq)+'_pol_'+str(pol)
    data_path='/home/dzli/ism/data/'

    if freq=='00':
        f = 314.5
    elif freq=='01':
        f = 322.5
    elif freq=='02':
        f = 330.5
    else:
        f = 338.5

    dy = sp.readdy(name=infile,path=inpath)
    print 'read in dy: ',dy.shape
    lenx = dy.shape[1]

    ic = fftshift(fft2(dy))
    del dy
    sp.save_conj(ic,outname+'.rebint',path=outpath)

    #calculate average S/N
    print 'get arclet'
    delay = (np.arange(leny)-np.floor(leny/2.))/nubox
    doppler =(np.arange(lenx)-np.floor(lenx/2.))/tbox			
    doppler1,doppler2,delay1,delay2 = [-16.8*f/f0,-16.8*f/f0,0.146,0.136] 
    mask = arc.get_mask(leny,lenx,doppler,delay,doppler1,doppler2,delay1,delay2)
    print 'non zeron mask: ',np.count_nonzero(mask)
    mean = np.mean(ic)
    signal =  np.sum(np.absolute(ic*mask-mean))/np.count_nonzero(mask)
    noise = np.average(np.absolute(ic-mean))
    print 'noise total',noise
    delaycut,dopplercut=[0.02,0.5]
    noisebox=[delaycut,1.0,-80,-40]
    noise=sp.get_noise(ic,noisebox,tbox=tbox,nubox=nubox,lenx=lenx,leny=leny)
    sn = signal/noise
    print 'mean ic, signal, noise, S/N \n'
    print mean, signal, noise, sn

    plot=1
    if plot==1:
        dpix,dpiy=300,4000
        sp.plot_spec(ic,spec_type='sec',name=outname, show= 'continue',tbox=tbox,ybins=10,dpix=dpix,dpiy=dpiy,vmin=5,vmax=9.5)
        #plot arclet to see whether the position matches for all gates
        doppler =(np.arange(dpix*2)-dpix)/tbox			
        sp.plot_parabola(doppler,doppler1,delay1,type='sec')
        sp.plot_parabola(doppler,doppler2,delay2,type='sec')
        pdf.savefig(1)
        plt.clf()

    savedata=0
    if savedata==1:
        print 'write to: ',data_path+'gate_weights.dat'
        dat = open (data_path+'gate_weights.dat','a')
        dat.write(outname+'\n')
        dat.write('mean_real,mean_imag,signal,noise, S/N: \n')
        json.dump([mean.real.float32,mean.imag.float32,signal,noise,sn],dat)
        dat.write('\n')
        dat.close()
    return sn



def readdy( baseline=257,freq='00',pol='LL',Gate=0):
	print 'Band selected: ' + freq, 'polar: '+pol 
	inpath = '/mnt/scratch-lustre/simard/B0834_2012/Gate'+str(Gate)+'/'
	infile = 'gb057_1.input_baseline'+str(baseline)+'_freq_'+str(freq)+'_pol_'+pol+'.rebint'
	dy = sp.readdy(name=infile,path=inpath)
	print 'read in dy: ',dy.shape
	return dy

def sn_dy(num_freq=0,pol_count=2,gate_count=3):
    pol=['LL','RR']
    freq=['00','01','02','03']
    Gate=[0,1,2]
    sn = np.zeros((pol_count,gate_count))
    for num_pol in xrange(0,pol_count):

        figurename = 'multigates_arclet_freq_'+freq[num_freq]+'pol_'+pol[num_pol]+'.pdf'
        pdf = pgs.PdfPages(figure_path+figurename)

        for num_Gate in xrange(0,gate_count):
            sn[num_pol,num_Gate]=gate_weight(freq=freq[num_freq],pol=pol[num_pol],Gate=Gate[num_Gate],pdf=pdf)
        pdf.close()
        print '============'
        print 'freq'+freq[num_freq]+'pol[L,R]'+' Gate012'
        print sn
    return sn

def sum_dy(freq_count=4,pol_count=2,gate_count=3,lenx=1321,leny=16384):
    outpath = '/mnt/raid-cita/dzli/gb057/multigate/'
    outname='gb057_baseline257_gate012_snvary_'
    pol=['LL','RR']
    freq=['00','01','02','03']
    Gate=[0,1,2]
    sn=np.zeros((2,3)) #polar,gates
    sn[0,:]=[4.3,3.6,2.8] #left
    sn[1,:]=[3.9,3.5,2.6]
    sn[0,:]=[4.1,3.5,2.7]
    sn[1,:]=sn[0,:]
    for num_freq in xrange(0,freq_count):
        sn = sn_dy(num_freq=num_freq)
	sn[0,:]=sum(sn)/sn.shape[0]
	sn[1,:]=sn[0,:] 
        for num_pol in xrange(0,pol_count):
            dytotal=np.zeros((leny,lenx))
            for num_Gate in xrange(0,gate_count):
                dy=readdy(freq=freq[num_freq],pol=pol[num_pol],Gate=Gate[num_Gate])
                dytotal=dytotal+dy*sn[num_pol,num_Gate]
            name=outname+'freq_'+freq[num_freq]+'_pol_'+pol[num_pol]+'.rebint'
            sp.save_dy(dytotal,name,path=outpath)
            ic=sp.conj_spec(dytotal)
            sp.save_conj(ic,name,path=outpath)
            del dytotal,ic
    return	
