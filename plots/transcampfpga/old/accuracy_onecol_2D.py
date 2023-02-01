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
import matplotlib


files=[["Speech","../../results/transcampfpga/audio/accuracy_audio_w16384_e%d_m%d_t65536.txt"],
       ["ECG", "../../results/transcampfpga/e0103/accuracy_e0103_w512_e%d_m%d_t262144.txt"],
       ["Power", "../../results/transcampfpga/power/accuracy_power_w1536_e%d_m%d_t262144.txt" ],
       ["Seismology", "../../results/transcampfpga/seismology/accuracy_seismology_w64_e%d_m%d_t65536.txt"],
       ["IMU", "../../results/transcampfpga/imu/accuracy_imu_w256_e%d_m%d_t65536.txt"],
       ["EPG", "../../results/transcampfpga/epg/accuracy_epg_w16384_e%d_m%d_t65536.txt"]]

#Exponentes
exp = [2,3,4,5,6,7,8]
#Mantisas
man = [2,5,7,10,13,16,20,23]


cmap = matplotlib.cm.get_cmap('CMRmap')
colors = (cmap(0),cmap(20),cmap(70),cmap(110),cmap(140),cmap(160),cmap(200))

#Line width
lw=2
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]

#X,Y = np.meshgrid(exp,man)
#X = np.arange(0, len, 1)
ZMff = np.empty_like(man)   #Motifs ff
ZMff10 = np.empty_like(man) #Motifs ff+-10
ZDff = np.empty_like(man)   #Discords ff
ZDff10 = np.empty_like(man) #Discords ff+-10
ZMfloat = 0. #Motif float
ZDfloat = 0. #Discords float

fig = plt.figure(figsize=(7,10))
for i in range(len(files)):
  ax = fig.add_subplot(3,2,i+1)
  ax.set_aspect("auto")
  
  for j in range(len(exp)):
    for k in range(len(man)):
      tmp = pst.StatsFile(files[i][1] % (exp[j],man[k]), False)
      tmp.ParseMotifsAccu()
      tmp.ParseDiscordsAccu()
      # 4 porcentajes: accu w.r.t float, float+-10, ff, ff+-10
      m = tmp.GetMotifsAccu()
      d = tmp.GetDiscordsAccu()
      ZMfloat = m[0]
      ZMff[k] = m[2]
      ZMff10[k] = m[3]
      ZDfloat = d[0]
      ZDff[k] = d[2]
      ZDff10[k] = d[3]
  

    ax.plot(man,ZMff, color=colors[j], linewidth=lw, alpha = 0.8, label=('Motifs exp %i' % exp[j]))
  
  #ax.plot_wireframe(X,Y,ZMff, color=motifColor, linewidth=lw, alpha = 0.8, label='Motifs')
  #ax.plot_wireframe(X,Y,ZDff, color=discordColor, linewidth=lw, alpha = 0.8, label='Discords')

  #ax.scatter(8,23,ZMfloat,c=motifColor,edgecolors=motifColor,marker='o')
  #ax.text(8-8./50,23,ZMfloat,"Single", ha="right", va="bottom")
  #ax.scatter(8,23,ZDfloat,c=discordColor,edgecolors=discordColor,marker='o')
  print("%s motifs and discords." % files[i][0])
  print(ZMff)
  print(ZDff)
  #ax.view_init(30,-70) #elevación, azimut
  ax.set_title(files[i][0], y=0.98)
  ax.set_xlabel("Mantissa")
  ax.set_ylabel("Accuracy")
  #ax.set_zlabel("Top-1000 Accuracy")
  ax.set_xticks(man)
  #ax.set_yticks(man+[23])
  #ax.set_zticks([0,20,40,60,80,100])
  ax.invert_xaxis()
  if i == 0:
    ax.legend(frameon=False, loc='upper right', fontsize = 'small')

plt.savefig("./accuracyTranSCAMPfpga.pdf", format='pdf')


