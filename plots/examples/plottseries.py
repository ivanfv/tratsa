#!/usr/bin/python
# -*- coding: utf-8 -*-
import csv     #Para leer ficheros de comma separated values
import sys
import os      #Paquete para acceder a comandos del OS
import fnmatch #Paquete para Unix filename pattern matching
# Importo la clase que lee los ficheros de estadísticas de gems
#sys.path.append('/home/quislant/Ricardo/Research/myPyLib/')
import numpy as np #Fundamental package for scientific computing with Python
import matplotlib.pyplot as plt #Librería gráfica
import matplotlib
#Importo la clase StatsFile de mi librería de parsing
#import parseStatsHaswell as psf
from subprocess import call #Para llamar a un comando del shell

tseries = "./e0103.txt"
mprofDouble = "./double.txt"
mprofFloat = "./float.txt"

tseries_2 = "./seismology.txt"
mprofDouble_2 = "./double_2.txt"
mprofFloat_2 = "./float_2.txt"

print("#### Leyendo la serie temporal 1")
tseriesVal = []
with open(tseries, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = 0
  for row in csv_reader:
    tseriesVal.append(float(row[-1])) #Hay dos espacios en la serie
    line_count += 1
  print("Elementos: %d" % line_count)
tseriesCount = line_count

print("#### Leyendo la serie temporal 2")
tseriesVal_2 = []
with open(tseries_2, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = 0
  for row in csv_reader:
    tseriesVal_2.append(float(row[-1])) #Hay dos espacios en la serie
    line_count += 1
  print("Elementos: %d" % line_count)
tseriesCount_2 = line_count

print("#### Leyendo el vector matrix double 1")
mprofDVal = []
with open(mprofDouble, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = -1
  for row in csv_reader:
    if line_count == -1: # Obvio el tiempo de ejecución
      line_count += 1
    else:
      mprofDVal.append(float(row[0])) # El tercer valor es el mprofile
      line_count += 1
  print("Elementos: %d" % line_count)
mprofDCount = line_count


print("#### Leyendo el vector matrix double 2")
mprofDVal_2 = []
with open(mprofDouble_2, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = -1
  for row in csv_reader:
    if line_count == -1: # Obvio el tiempo de ejecución
      line_count += 1
    else:
      mprofDVal_2.append(float(row[0])) # El tercer valor es el mprofile
      line_count += 1
  print("Elementos: %d" % line_count)
mprofDCount_2 = line_count


print("#### Leyendo el vector matrix profile float")
mprofFVal = []
with open(mprofFloat, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = -1
  for row in csv_reader:
    if line_count == -1: # Obvio el tiempo de ejecución
      line_count += 1
    else:
      mprofFVal.append(float(row[0])) # El tercer valor es el mprofile
      line_count += 1
  print("Elementos: %d" % line_count)
mprofFCount = line_count


print("#### Leyendo el vector matrix profile float 2")
mprofFVal_2 = []
with open(mprofFloat_2, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=' ')
  line_count = -1
  for row in csv_reader:
    if line_count == -1: # Obvio el tiempo de ejecución
      line_count += 1
    else:
      mprofFVal_2.append(float(row[0])) # El tercer valor es el mprofile
      line_count += 1
  print("Elementos: %d" % line_count)
mprofFCount_2 = line_count

#Texto en latex y fuente
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

xini = 100500
nele = 1500

xini_2 = 103150
nele_2 = 600

cmap = matplotlib.cm.get_cmap('CMRmap')
colors = (cmap(0),cmap(30),cmap(100),cmap(200))

x = [i for i in range(xini,xini+nele)]

x_2 = [i for i in range(xini_2,xini_2+nele_2)]
plt.figure(figsize=(9,4))
#En el primer subplot va la tseries
ax = plt.subplot(2,2,1)
plt.plot(x, tseriesVal[xini:xini+nele], '-', linewidth=1.5, color=colors[1])
plt.ylabel("Amplitude")
#ax.yaxis.set_label_coords(-0.05,0.5)
plt.grid(axis='x')
plt.yticks(range(0, int(max(tseriesVal[xini:xini+nele]))+2 ))
plt.xticks(range(xini,xini+nele+1,500),())

ax = plt.subplot(2,2,3)
#plt.subplot(2,2,3)
plt.plot(x, mprofDVal[xini:xini+nele], '-', linewidth=1.5, color=colors[3], label='double')
plt.plot(x, mprofFVal[xini:xini+nele], ':', linewidth=1.5, color=colors[2], label='single')
plt.ylabel("Profile")
ax.yaxis.set_label_coords(-0.08,0.5)
plt.grid(axis='x')
plt.yticks(range(0, int(max(mprofDVal[xini:xini+nele]))+2, 4 ))
plt.xticks(range(xini,xini+nele+1,500),('100K','100.5K','101K','101.5K','102K','102.5K'))
plt.xlabel("Data Points (ECG)")
plt.legend(frameon=False, labelspacing=0.5, loc='upper left', fontsize=15)


ax = plt.subplot(2,2,2)
plt.plot(x_2, [x+40 for x in tseriesVal_2[xini_2:xini_2+nele_2]], '-', linewidth=1.5, color=colors[1])
#plt.ylabel("T")
#ax.yaxis.set_label_coords(-0.05,0.5)
plt.grid(axis='x')
plt.yticks(np.arange(0, 65, step=20))
plt.xticks(range(xini_2,xini_2+nele_2+1,200),())

plt.subplot(2,2,4)
plt.plot(x_2, mprofDVal_2[xini_2:xini_2+nele_2], '-', linewidth=1.5, color=colors[3], label='double')
plt.plot(x_2, mprofFVal_2[xini_2:xini_2+nele_2], ':', linewidth=1.5, color=colors[2], label='single')
#plt.ylabel("P")
plt.grid(axis='x')
plt.yticks(np.arange(5, 7.5, step=0.8))
plt.xticks(range(xini_2,xini_2+nele_2+1,200),('103K','103.2K','103.4K','103.6K','102K','102.5K'))

plt.xlabel("Data Points (Seismology)")
#plt.legend(frameon=False, labelspacing=0.5, loc='lower right', fontsize=15)


plt.savefig("ECGe0103DoubleFloat.pdf", format='pdf', bbox_inches='tight')
#plt.savefig("ECGe0103DoubleFloat.png", format='png', bbox_inches='tight')
#plt.show()
