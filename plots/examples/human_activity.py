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


files=[["Human Activity(7,13)", "../../results/transcamp/result_human_activity-MPIII-SVC_w120_s1_d0_t256_cov7-13_corr7-13_stat7-13_prof7-13_28702.csv"],
       ["Human Activity(7,10)", "../../results/transcamp/result_human_activity-MPIII-SVC_w120_s1_d0_t256_cov7-10_corr7-10_stat7-10_prof7-10_28443.csv"],
       ["Human Activity(5,13)", "../../results/transcamp/result_human_activity-MPIII-SVC_w120_s1_d0_t256_cov5-13_corr5-13_stat5-13_prof5-13_24560.csv"],
       ["Human Activity(5,7)",  "../../results/transcamp/result_human_activity-MPIII-SVC_w120_s1_d0_t256_cov5-7_corr5-7_stat5-7_prof5-7_24043.csv"]]

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
  tmp = pst.StatsFile(files[i][1], False)
  tmp.ParseProfLen()
  tmp.ParseSeries()

  mp_d = np.array(tmp.GetDoubleMP())
  mp_f = np.array(tmp.GetFloatMP())
  mp_ff = np.array(tmp.GetFlexMP())

  
  X = np.arange(0.0, len(mp_d), 1)
  ax = fig.add_subplot(2,2,i+1)
  ax.set_aspect("auto")
  ax.plot(X[1000:1200],mp_d[1000:1200,0], color=motifColor,   linewidth=lw, alpha = 0.8, label='MP Double')
  ax.plot(X[1000:1200],mp_ff[1000:1200,0], '--',color=flexfloatColor, linewidth=lw, alpha = 0.8, label='MP Trans.')

  ax.set_title(files[i][0], y=0.98)

  ax.legend(fontsize=10,ncol=1, columnspacing=1, loc='upper left')

	  


plt.savefig("./human_activity_MP_profiles.pdf", format='pdf')


