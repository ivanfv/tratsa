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
import matplotlib.colors as colors
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
#man = [2,5,7,10,13,16,20,23]
man = [7,10,13,16,20,23]

barx_size = 0.3
bary_size = 1

#Line width
lw=2
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]

X,Y     = np.meshgrid(exp,man, copy=False)
ZMff    = np.empty_like(X)   #Motifs ff
ZMff10  = np.empty_like(X) #Motifs ff+-10
ZDff    = np.empty_like(X)   #Discords ff
ZDff10  = np.empty_like(X) #Discords ff+-10
ZMfloat = 0. #Motif float
ZDfloat = 0. #Discords float
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
      ZMfloat     = m[0]
      ZMff[k,j]   = m[2]
      ZMff10[k,j] = m[3]
      ZDfloat     = d[0]
      ZDff[k,j]   = d[2]
      ZDff10[k,j] = d[3]
  
  ax = fig.add_subplot(3, 2, i + 1,projection='3d')
  
  ax.set_aspect("auto")

  x_coords = X.ravel() - barx_size / 2
  y_coords = Y.ravel() - bary_size / 2
  z_coords = ZMff.ravel() 
  bottom   = np.zeros_like(z_coords)
 
  r_color = np.zeros_like(z_coords)
  g_color = z_coords
  b_color = np.zeros_like(z_coords)
  print(z_coords)
  print("here")


  dz = ZMff
  offset  = dz + np.abs(dz.min())
  fracs   = offset.astype(float)/offset.max()
  norm    = colors.Normalize(fracs.min(), fracs.max())
  colorss = cm.Greens(norm(fracs))

  CS = ax.bar3d(x_coords, y_coords,bottom, dx=barx_size, dy=bary_size ,dz=z_coords, shade=True, edgecolor = "black", color=colorss.reshape(-1,4)) 
  
  plt.xlim(8, 2)
  ax.set_xticks([2,3,4,5,6,7,8])
  plt.ylim(7, 23)
  ax.set_yticks([7,10,13,16,20,23])  
 
  ax.set_zlim(0, 100)
 

  print("%s motifs." % files[i][0])
  print(ZMff)


  ax.set_title(files[i][0], y=0.98)
  ax.set_xlabel("Exponent")
  ax.set_ylabel("Mantissa")
  if((i + 1) % 2 == 0):
    ax.set_zlabel("\% Top-1000 Mot.")


#fig.subplots_adjust(bottom=0.8)
#cbar_ax = fig.add_axes([0.1, -0.02, 0.865, 0.03])
#plt.colorbar(CS, cax=cbar_ax,orientation="horizontal")

#cbar_ax.text(3, 30, "\% Top-1000 Motif Accuracy")

plt.savefig("./accuracyMotTranSCAMPfpga.pdf", format='pdf',bbox_inches='tight')
