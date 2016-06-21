from numpy import *
import numpy as np
from pandas import Series, DataFrame
import pandas as pd
import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import colors
import pylab 
import scipy.stats as stats
from scipy.special import erf
import sys

def gaussian(x,x0,s):
	return exp(-(x-x0)**2/(2*s**2))

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

if sys.argv[0]:
	window_size_signal=0
	if len(sys.argv)>1:
		window_size_signal=int(sys.argv[1])
	window_size_background=3
	if len(sys.argv)>2:
		window_size_background=int(sys.argv[2])
	filter_type='flat'
	if len(sys.argv)>3:
		filter_type=sys.argv[3]
	print str(window_size_signal),str(window_size_background),filter_type

points=1000
x = linspace(-1,1,points)
ys = zeros(points)
ys+=gaussian(x,0,0.1)
ys+=0.75*gaussian(x,0.3,0.01)
ys+=1.2*gaussian(x,0.7,0.01)
ys+=0.9*gaussian(x,0.8,0.01)
yn=0.2*random.normal(size=len(x))
ysn=ys+yn

for ws in range(0, 20):
	for wb in range(0, 20):
		for ft in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
			ysn_smooth_signal=smooth(ysn,ws,ft)
			ysn_smooth_background=smooth(ysn,wb,ft)
			corrected_signal = ysn_smooth_signal - ysn_smooth_background
			maxinds = corrected_signal.argsort()[-3:][::-1]
			xvals = [x[n] for n in maxinds]
			if any([t>0.21 and t<0.31 for t in xvals]) and any([t>0.61 and t<0.71 for t in xvals]) and any([t>0.71 and t<0.81 for t in xvals]):
				print (ws, wb, ft)


"""

ysn_smooth_signal=smooth(ysn,window_size_signal,filter_type)
ysn_smooth_background=smooth(ysn,window_size_background,filter_type)

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ys,c='r')
ax1.set_ylim([1.2*min(ys),1.2*max(ys)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,yn,c='black')
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('noise',dpi=300,bbox_inches='tight')
fig.clf()
	
fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn,c='black')
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks_plus_noise',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn,c='black',lw=1)
ax1.plot(x,ysn_smooth_signal,c='red',lw=3,alpha=0.5)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks_plus_noise_smooth_signal',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn,c='black',lw=1)
ax1.plot(x,ysn_smooth_background,c='r',lw=3,alpha=0.5)
ax1.set_ylim([1.2*min(ysn),1.2*max(ysn)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks_plus_noise_smooth_background',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn-ysn_smooth_background,c='black')
ax1.set_ylim([1.2*min(ysn-ysn_smooth_background),1.2*max(ysn-ysn_smooth_background)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks_plus_noise_background_subtracted',dpi=300,bbox_inches='tight')
fig.clf()

fig, (ax1) = plt.subplots(1,figsize=(6,6))
ax1.plot(x,ysn_smooth_signal-ysn_smooth_background,c='black')
ax1.set_ylim([1.2*min(ysn_smooth_signal-ysn_smooth_background),1.2*max(ysn_smooth_signal-ysn_smooth_background)])
ax1.set_xlim([min(x),max(x)])
fig.savefig('peaks_plus_noise_smooth_background_subtracted',dpi=300,bbox_inches='tight')
fig.clf()

"""
