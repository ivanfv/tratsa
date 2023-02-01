#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np #Fundamental package for scientific computing with Python
import matplotlib.pyplot as plt #Librería gráfica
#Importo la clase StatsFile de mi librería de parsing
import parseStatsTransPrecision as pst
from subprocess import call #Para llamar a un comando del shell
from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import pylab


files=[["Power(6,20) r=16K", "../../results/transcampfpga/power_outputs/MP_power_scaled.txt_1536_5_16_16384.txt"],
       ["Power(6,20) r=64K", "../../results/transcampfpga/power_outputs/MP_power_scaled.txt_1536_5_16_65536.txt"],
       ["Power(6,20) r=256K", "../../results/transcampfpga/power_outputs/MP_power_scaled.txt_1536_5_16_262144.txt"],
       ["Power(6,20) r=off",  "../../results/transcampfpga/power_outputs/MP_power_scaled.txt_1536_5_16_2100000.txt"]]

#Exponentes
exp = [7]
#Mantisas
man = [16]
#Recalcs
recal = [16384]


#Line width
lw=1
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
flexfloatColor = [191./256,63./256,63./256]



fig = plt.figure(figsize=(7,5))
for i in range(len(files)):
  print(files[i][1])

  mp_d = pylab.loadtxt("../../mp_tseries_ref/power/MP_power_scaled_w1536_t2100000_double.txt");
  #mp_f = pylab.loadtxt("../../mp_tseries_ref/epg/MP_epg_scaled_w16384_t2100000_double.txt");
  mp_ff = pylab.loadtxt(files[i][1]);

  
  X = np.arange(0.0, len(mp_d), 1)
  ax = fig.add_subplot(2,2,i+1)
  ax.set_aspect("auto")
  ax.plot(X[1700000:1700050],mp_d[1700000:1700050,0], color=motifColor,   linewidth=lw, alpha = 0.8, label='MP Double')
  ax.plot(X[1700000:1700050],mp_ff[1700000:1700050,0], '--',color=flexfloatColor, linewidth=lw, alpha = 0.8, label='MP Trans.')

  ax.set_title(files[i][0], y=0.98)

  ax.legend(fontsize=10,ncol=1, columnspacing=1, loc='upper left')

	  


plt.savefig("./power_MP_profiles_recalc.pdf", format='pdf')


