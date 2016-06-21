from numpy import *
from pandas import Series, DataFrame
import pandas as pd
import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import colors
import pylab 
import scipy.stats as stats
from scipy.special import erf

def gaussian(x,x0,s):
	return exp(-(x-x0)**2/(2*s**2))
	
points=1000
x = linspace(-1,1,points)
ys = zeros(points)
ys+=gaussian(x,0,0.1)
yn=0.2*random.normal(size=len(x))
ysn=ys+yn

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ys,c='r')
ax1.set_ylim([1.2*min(ys),1.2*max(ys)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian',dpi=300,bbox_inches='tight')
fig.clf()

fftys=fft.rfft(ys)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftys),c='gray')
ax1.set_xlim([0,50])
fig.savefig('gaussian_fft',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,yn,c='r')
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('noise',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
stats.probplot(yn, dist="norm", plot=pylab)
fig.savefig('noise_qqplot',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.hist(yn,bins=100,histtype='step',color='black',normed=False)
ax1.set_xlim([1.2*min(ysn),1.2*max(ysn)])
fig.savefig('noise_hist',dpi=300,bbox_inches='tight')
fig.clf()
	
fftyn=fft.rfft(yn)
fftysn=fft.rfft(ysn)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftysn),c='white')
ax1.plot(abs(fftyn),c='gray')
ax1.set_xlim([0,50])
fig.savefig('noise_fft',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn,c='r')
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise',dpi=300,bbox_inches='tight')
fig.clf()

fftysn=fft.rfft(ysn)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftysn),c='gray')
ax1.set_xlim([0,50])
fig.savefig('gaussian_noise_fft',dpi=300,bbox_inches='tight')
fig.clf()

#Smoothing
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<1:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	if window == 'flat': 
		w=ones(2*window_len+1,'d')
	else:
		w=eval(window+'(2*window_len+1)')
	y=convolve(w/w.sum(),x,'valid')
	return r_[[y[0]]*int(window_len),y,[y[-1]]*int(window_len)]

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn,c='black')
ax1.plot(x,ys,c='r',lw=2)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise_smooth_no',dpi=300,bbox_inches='tight')
fig.clf()

ysn_smooth_flat_3=smooth(ysn,3,'flat')
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn_smooth_flat_3,c='black')
ax1.plot(x,ys,c='r',lw=2)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise_smooth_flat_3',dpi=300,bbox_inches='tight')
fig.clf()

fftysn_smooth_flat_3=fft.rfft(ysn_smooth_flat_3)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftysn_smooth_flat_3),c='gray')
ax1.plot(abs(fftysn),c='magenta')
ax1.set_xlim([0,50])
fig.savefig('gaussian_noise_smooth_flat_3_fft',dpi=300,bbox_inches='tight')
fig.clf()

ysn_smooth_flat_11=smooth(ysn,11,'flat')
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn_smooth_flat_11,c='black')
ax1.plot(x,ys,c='r',lw=2)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise_smooth_flat_11',dpi=300,bbox_inches='tight')
fig.clf()

ysn_smooth_flat_31=smooth(ysn,31,'flat')
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn_smooth_flat_31,c='black')
ax1.plot(x,ys,c='r',lw=2)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise_smooth_flat_31',dpi=300,bbox_inches='tight')
fig.clf()

fftysn_smooth_flat_31=fft.rfft(ysn_smooth_flat_31)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftysn_smooth_flat_31),c='gray')
ax1.plot(abs(fftysn),c='magenta')
ax1.set_xlim([0,50])
fig.savefig('gaussian_noise_smooth_flat_31_fft',dpi=300,bbox_inches='tight')
fig.clf()

ysn_smooth_flat_91=smooth(ysn,91,'flat')
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn_smooth_flat_91,c='black')
ax1.plot(x,ys,c='r',lw=2)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('gaussian_noise_smooth_flat_91',dpi=300,bbox_inches='tight')
fig.clf()

fftysn_smooth_flat_91=fft.rfft(ysn_smooth_flat_91)
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(abs(fftysn_smooth_flat_91),c='gray')
ax1.plot(abs(fftysn),c='magenta')
ax1.set_xlim([0,50])
fig.savefig('gaussian_noise_smooth_flat_91_fft',dpi=300,bbox_inches='tight')
fig.clf()

