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

import matplotlib.cm as cm


files=[["Speech","../../results/transcampfpga/audio/accuracy_audio_w16384_e%d_m%d_t65536.txt"],
       ["ECG", "../../results/transcampfpga/e0103/accuracy_e0103_w512_e%d_m%d_t262144.txt"],
       ["Power", "../../results/transcampfpga/power/accuracy_power_w1536_e%d_m%d_t262144.txt" ],
       ["Seismology", "../../results/transcampfpga/seismology/accuracy_seismology_w64_e%d_m%d_t65536.txt"],
       ["IMU", "../../results/transcampfpga/imu/accuracy_imu_w256_e%d_m%d_t65536.txt"],
       ["EPG", "../../results/transcampfpga/epg/accuracy_epg_w16384_e%d_m%d_t65536.txt"]]

#Exponentes
exp = [2,3,4,5,6,7,8]
#Mantisas
#man = [2,5,7,10,12,14,16,18,20,22]
man = [2,5,7,10,13,16,20,23]

#Line width
lw=2
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]

X,Y = np.meshgrid(exp,man)
ZMff = np.empty_like(X)   #Motifs ff
ZMff10 = np.empty_like(X) #Motifs ff+-10
ZDff = np.empty_like(X)   #Discords ff
ZDff10 = np.empty_like(X) #Discords ff+-10
ZMfloat = 0. #Motif float
ZDfloat = 0. #Discords float
#fig = plt.figure(figsize=(9,7))
fig = plt.figure(figsize=(7,10))
for i in range(len(files)):
  for j in range(len(exp)):
    for k in range(len(man)):
      tmp = pst.StatsFile(files[i][1] % (exp[j],man[k]), False)
      tmp.ParseMotifsAccu()
      tmp.ParseDiscordsAccu()
      # 4 porcentajes: accu w.r.t float, float+-10, ff, ff+-10
      m = tmp.GetMotifsAccu()
      d = tmp.GetDiscordsAccu()
      ZMfloat = m[0]
      ZMff[k,j] = m[2]
      ZMff10[k,j] = m[3]
      ZDfloat = d[0]
      ZDff[k,j] = d[2]
      ZDff10[k,j] = d[3]
  
  ax = fig.add_subplot(3,2,i+1)
  
  
#  ax.set_aspect(0.25)
  ax.set_aspect("auto")
  CS =ax.contourf(X,Y,ZDff,[0,10,20,30,40,50,60,70,80,90,95,100]) 
  manual_locations = [(7, 12), (7, 13), (7, 14), (7, 15), (7, 16)]
  #ax.clabel(CS, inline=1, fontsize=10)
   
  CB = plt.colorbar(CS, shrink=0.8, extend='both')
  

#, manual=manual_locations
  plt.xlim(8, 4)
  plt.ylim(10, 23)
  #ax.plot_wireframe(X,Y,ZMff, color=motifColor, linewidth=lw, alpha = 0.8, label='Motifs')
  #ax.plot_wireframe(X,Y,ZDff, color=discordColor, linewidth=lw, alpha = 0.8, label='Discords')
#  ax.plot_surface(X,Y,ZMff, cmap=plt.get_cmap('Blues'), linewidth=lw, alpha = 0.3,antialiased=True)
#  ax.plot_surface(X,Y,ZDff, cmap=plt.get_cmap('Reds'), linewidth=lw, alpha = 0.3,antialiased=True)
  #ax.scatter(8,23,ZMfloat,c=motifColor,edgecolors=motifColor,marker='o')
  #ax.text(8-8./50,23,ZMfloat,"Single", ha="right", va="bottom")
  #ax.scatter(8,23,ZDfloat,c=discordColor,edgecolors=discordColor,marker='o')
  print("%s motifs and discords." % files[i][0])
  print(ZMff)
  print(ZDff)
  #ax.view_init(30,-70) #elevación, azimut
  ax.set_title(files[i][0], y=0.98)
  ax.set_xlabel("Exponent")
  ax.set_ylabel("Mantissa")
  #ax.set_zlabel("Top-1000 Accuracy")
  #ax.set_xticks(exp)
  ax.set_yticks([10,13,16,20,23])
  #ax.set_zticks([0,20,40,60,80,100])
  ax.invert_xaxis()
  if i == 0:
    ax.legend(frameon=False, loc='upper right', fontsize = 'small')


plt.savefig("./accuracyMotTranSCAMPfpga.pdf", format='pdf')
