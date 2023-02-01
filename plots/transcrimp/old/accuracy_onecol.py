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

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif', size=14)

#files=[["Audio","../results/result_audio-MPIII-SVD_w200_s1_d0_t256_he%d_hm%d_le%d_lm%d_*.csv"],
#       ["ECG", "../results/result_e0103_n180000_w500_s1_d0_t256_he%d_hm%d_le%d_lm%d_*.csv"],
#       ["Power", "../results/result_power-MPIII-SVF_n180000_w1325_s0.1_d0_t256_he%d_hm%d_le%d_lm%d_*.csv" ],
#       ["Seismology", "../results/result_seismology-MPIII-SVE_n180000_w50_s0.01_d0_t256_he%d_hm%d_le%d_lm%d_*.csv"],
#       ["Human Activity", "../results/result_human_activity-MPIII-SVC_w120_s1_d0_t256_he%d_hm%d_le%d_lm%d_*.csv"],
#       ["Penguin Behavior", "../results/result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he%d_hm%d_le%d_lm%d_*.csv"]]

files=[["Audio","../../../transcrimp/results/result_audio-MPIII-SVD_w200_s1.00_t256_he%d_hm%d_le%d_lm%d_*.csv"],
       ["ECG", "../../../transcrimp/results/result_e0103_n180000_w500_s1.00_t256_he%d_hm%d_le%d_lm%d_*.csv"],
       ["Power", "../../../transcrimp/results/result_power-MPIII-SVF_n180000_w1325_s0.10_t256_he%d_hm%d_le%d_lm%d_*.csv" ],
       ["Seismology", "../../../transcrimp/results/result_seismology-MPIII-SVE_n180000_w50_s0.01_t256_he%d_hm%d_le%d_lm%d_*.csv"],
       ["Human Activity", "../../../transcrimp/results/result_human_activity-MPIII-SVC_w120_s1.00_t256_he%d_hm%d_le%d_lm%d_*.csv"],
       ["Penguin Behavior", "../../../transcrimp/results/result_penguin_sample_TutorialMPweb_w800_s1.00_t256_he%d_hm%d_le%d_lm%d_*.csv"]]

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
      tmp = pst.StatsFile(files[i][1] % (exp[j],man[k],exp[j],man[k]), False)
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
  
  ax = fig.add_subplot(3,2,i+1, projection='3d')
#  ax.set_aspect(0.25)
  ax.set_aspect("auto")
  ax.plot_wireframe(X,Y,ZMff, color=motifColor, linewidth=lw, alpha = 0.8, label='Motifs')
  ax.plot_wireframe(X,Y,ZDff, color=discordColor, linewidth=lw, alpha = 0.8, label='Discords')
#  ax.plot_surface(X,Y,ZMff, cmap=plt.get_cmap('Blues'), linewidth=lw, alpha = 0.3,antialiased=True)
#  ax.plot_surface(X,Y,ZDff, cmap=plt.get_cmap('Reds'), linewidth=lw, alpha = 0.3,antialiased=True)
  ax.scatter(8,23,ZMfloat,c=motifColor,edgecolors=motifColor,marker='o')
  ax.text(8-8./50,23,ZMfloat,"Single", ha="right", va="bottom")
  ax.scatter(8,23,ZDfloat,c=discordColor,edgecolors=discordColor,marker='o')
  print("%s motifs and discords." % files[i][0])
  print(ZMff)
  print(ZDff)
  ax.view_init(30,-70) #elevación, azimut
  ax.set_title(files[i][0], y=0.98)
  ax.set_xlabel("Exponent")
  ax.set_ylabel("Mantissa")
  ax.set_zlabel("Top-100 Accuracy")
  ax.set_xticks(exp)
  ax.set_yticks(man+[23])
  ax.set_zticks([0,20,40,60,80,100])
  ax.invert_xaxis()
  if i == 0:
    ax.legend(frameon=False, loc='upper right', fontsize = 'small')
    #ax.text2D(-0.11, -0.05, 'TranSCRIMP', rotation=90 , fontsize=20)
  #ax.spines["right"].set_visible(False)
  #ax.spines["top"].set_visible(False)
  ##ax.spines["left"].set_linewidth(1)
  ##ax.spines["bottom"].set_linewidth(1)
  #ax.get_xaxis().tick_bottom()
  #ax.get_yaxis().tick_left()
  #ax.yaxis.grid(color=lgc, linewidth=1, linestyle="-")
  #ax.set_axisbelow(True) #Para que el grid no se muestre por encima de las barras
  #print ax.yaxis.label.get_fontname()
#plt.text(0, 0, 'SCAMP', fontdict=None, withdash=False, fontsize=12)
#plt.tight_layout(pad=5., w_pad=5., h_pad=5.)
#fig = plt.gcf()
#si = fig.get_size_inches()
#fig.set_size_inches(si[0]*3.5, si[1])
#bbox_extra_artists=(lgd,)
plt.savefig("./accuracyTranSCRIMP.pdf", format='pdf')
#plt.savefig("./accuracy.png", format='png', bbox_inches='tight')
#call("cp ../z_figs/SU.pdf ../../paper-TPDS/figs/", shell=True)
#plt.show()

