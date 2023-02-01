#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
# Importo la clase que lee los ficheros de estadisticas de gems
sys.path.append('/home/quislant/Ricardo/Research/myPyLib/')
import numpy as np #Fundamental package for scientific computing with Python
import matplotlib.pyplot as plt #Librería gráfica
#Importo la clase StatsFile de mi librería de parsing
import parseStatsTransPrecision as pst
from subprocess import call #Para llamar a un comando del shell

#Lo hago sólo para uno. Los resultados son proporcionales independientemente de la longitud de latime series.
#files=[["Audio","../results/result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le8_lm23_*.csv"],
#       ["ECG", "../results/result_e0103_n180000_w500_s1_d0_t256_he8_hm23_le8_lm23_*.csv"],
#       ["Power", "../results/result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm23_le8_lm23_*.csv" ],
#       ["Seis.", "../results/result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm23_le8_lm23_*.csv"],
#       ["Human", "../results/result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm23_le8_lm23_*.csv"],
#       ["Peng.", "../results/result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm23_le8_lm23_*.csv"]]

files=[["../results/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le8_lm23_*.csv", 8, 23],
       ["../results/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le8_lm7_*.csv", 8, 7],
       ["../results/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le5_lm10_*.csv", 5, 10]]
#Datos de energía del paper:
#Mach, S., Schuiki, F., Zaruba, F., Benini, L.: A 0.80 pj/flop, 1.24 tflop/sw 8-to-64bit transprecision floating-point unit for a 64 bit risc-v processor in 22nm fd-soi. In:2019 IFIP/IEEE 27th International Conference on Very Large Scale Integration(VLSI-SoC). pp. 95–98. IEEE (2019)
#####columnas### Double, single, half(8,7), half(5,10), halfhalf
#fila 1: FMA
#fila 2: Mult
#fila 3: Add
#fila 4: Comp
eScalar=[[26.7, 9.4, 4.4,   5,  2.5],
         [24.3, 8.5,   4, 4.5, 2.4],
         [14.1, 6.7, 2.9, 3.5,  1.9],
         [ 3.9, 2.5, 1.6, 1.6,  1.3]]
#Para SIMD no hay double por eso relleno con 0
eSIMD  =[[   0,   10, 3.425, 4.025,  1.6125],
         [   0, 8.85, 2.975, 3.475,   1.375],
         [   0,  6.9, 2.375,  2.75,  1.1125],
         [   0,    2, 0.775,  0.85,     0.4]]
#Line width
lw=2
bw=0.5
plt.rc('text', usetex=True) #Le digo que el procesador de textos sea latex (se pueden usar fórmulas)
plt.rc('font', family='sans-serif') #Le digo que la fuente por defecto sea sans-serif
plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automáticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]

#La energía de double w.r.t. double es 1
energy = [1]
for i in range(len(files)):
  if i==0:
    ######Primero saco los valores para todo a float
    tmp = pst.StatsFile(files[i][0], False)
    # (INV,ADD,SUB,MUL,DIV,FMA,CMP)
    #considero INV,ADD,SUB como sumas (FlexFloat considera acumulaciones como INV)
    #MUL y DIV como multiplicaciones
    tmp.ParseFlexFloatBreak(files[i][1],files[i][2])
    inst = tmp.GetFlexFloatBreak()
    print("FlexFloat<%d,%d> %s" % (files[i][1], files[i][2], inst))
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    doubleE = fmas*eScalar[0][0] + muls*eScalar[1][0] + adds*eScalar[2][0] + comps*eScalar[3][0]
    floatE = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    floatESIMD = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]
    energy.append(float(floatE)/doubleE)
    energy.append(float(floatESIMD)/doubleE)
  else:
    ######Segundo los valores para mixed con 8,7
    tmp = pst.StatsFile(files[i][0], False)
    tmp.ParseFlexFloatBreak(8,23) #Primero los float
    inst = tmp.GetFlexFloatBreak()
    print("FlexFloat<%d,%d>-<%d,%d> %s " % (8,23,files[i][1], files[i][2], inst)),
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    mixESIMD = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]
    tmp.ParseFlexFloatBreak(files[i][1],files[i][2]) #luego los 8,7
    inst = tmp.GetFlexFloatBreak()
    print(inst)
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE += fmas*eScalar[0][i+1] + muls*eScalar[1][i+1] + adds*eScalar[2][i+1] + comps*eScalar[3][i+1]
    mixESIMD += fmas*eSIMD[0][i+1] + muls*eSIMD[1][i+1] + adds*eSIMD[2][i+1] + comps*eSIMD[3][i+1]
    energy.append(float(mixE)/doubleE)
    energy.append(float(mixESIMD)/doubleE)

print (energy)
x=range(len(energy))
xshift=[i+0.25 for i in x]
plt.bar(xshift, energy, bw, color=motifColor, edgecolor=[i*0.8 for i in motifColor])
plt.axes().set_aspect(3)
plt.title('SCRIMP', fontsize=26, loc='right')
plt.ylabel('FPU Energy', fontsize=26)
yticks = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
plt.yticks([float(i) for i in yticks], yticks, fontsize=15)
plt.xticks([i+0.5 for i in x], ['Double', 'Single', 'Single\nSIMD', 'Mixed 8,7', 'Mixed 8,7\nSIMD','Mixed 5,10', 'Mixed 5,10\nSIMD'], rotation=32, fontsize=20)
plt.grid(True, axis='y')

plt.savefig("./energySCRIMP.pdf", format='pdf', bbox_inches='tight')
#plt.savefig("./accuracy.png", format='png', bbox_inches='tight')
#call("cp ../z_figs/SU.pdf ../../paper-TPDS/figs/", shell=True)
#plt.show()

