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


files=[["Speech(7,20)","../../results/transcampfpga/audio/accuracy_audio_w16384_e6_m20_t%d.txt"],
       ["ECG(7,13)", "../../results/transcampfpga/e0103/accuracy_e0103_w512_e7_m13_t%d.txt"],
       ["Power(5,16)", "../../results/transcampfpga/power/accuracy_power_w1536_e5_m16_t%d.txt"],
       ["Seismo.(6,16)", "../../results/transcampfpga/seismology/accuracy_seismology_w64_e6_m16_t%d.txt"],
       ["IMU(5,13)", "../../results/transcampfpga/imu/accuracy_imu_w256_e5_m13_t%d.txt"],
       ["EPG(6,20)", "../../results/transcampfpga/epg/accuracy_epg_w16384_e6_m20_t%d.txt"]]

#Exponentes
exp = [7]
#Mantisas
man = [16]
#Recalcs
recal = [16384, 65536, 262144, 2100000]


#Line width
lw=2
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]

X = [16384, 710922, 1405460, 2100000]
ZMff = np.empty_like(X)   #Motifs ff
ZDff = np.empty_like(X)   #Discords ff


fig = plt.figure(figsize=(7,5))
for i in range(len(files)):
  for l in range(len(recal)):
    tmp = pst.StatsFile(files[i][1] % (recal[l]), False)
    tmp.ParseMotifsAccu()
    tmp.ParseDiscordsAccu()
    # 4 porcentajes: accu w.r.t float, float+-10, ff, ff+-10
    m = tmp.GetMotifsAccu()
    d = tmp.GetDiscordsAccu()
    ZMff[l] = m[2]
    ZDff[l] = d[2]
  
  ax = fig.add_subplot(2,3,i+1)
  ax.set_aspect("auto")
  ax.plot(X,ZMff, color=motifColor,   linewidth=lw, alpha = 0.8, label='Motifs')
  ax.plot(X,ZDff, color=discordColor, linewidth=lw, alpha = 0.8, label='Discords')


  print("%s motifs and discords." % files[i][0])
  print(ZMff)
  print(ZDff)
  ax.set_title(files[i][0], y=0.98)

  ax.grid()

  ax.set_yticks([0,25,50,75,100])
  
  labels = [item.get_text() for item in ax.get_xticklabels()]
  labels[0] = '16K'
  labels[1] = '64K'
  labels[2] = '256K'
  labels[3] = 'off'
  #ax.set_xticklabels(labels)
  plt.xticks(X, labels[0:4], rotation=45)
  #ax.set_xticks(man)
  plt.tight_layout()

  if i == 0:
	  ax.set_ylabel("\% Top-1000 Accur.")
	  labels = [item.get_text() for item in ax.get_xticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_xticklabels(empty_string_labels)
	  ax.legend(frameon=True, loc='lower center', fontsize = 'small')

 
  if i == 1:

	  labels = [item.get_text() for item in ax.get_xticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_xticklabels(empty_string_labels)
	  
	  labels = [item.get_text() for item in ax.get_yticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_yticklabels(empty_string_labels)
	  
  if i == 2:
	  labels = [item.get_text() for item in ax.get_xticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_xticklabels(empty_string_labels)
	
	  labels = [item.get_text() for item in ax.get_yticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_yticklabels(empty_string_labels)

  if i == 3:
	  ax.set_ylabel("\% Top-1000 Accur.")
	  #ax.set_yticks([0,25,50,75,100])
	  ax.set_xlabel("Recalc. Factor")
	  #ax.set_xticks(man) 

  if i == 4:

	  ax.set_xlabel("Recalc. Factor")
	  labels = [item.get_text() for item in ax.get_yticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_yticklabels(empty_string_labels)

  if i == 5:
	  ax.set_xlabel("Recalc. Factor")
	  labels = [item.get_text() for item in ax.get_yticklabels()]
	  empty_string_labels = ['']*len(labels)
	  ax.set_yticklabels(empty_string_labels) 
	  


plt.savefig("./accuracyTranSCAMPfpga_tiling.pdf", format='pdf')


