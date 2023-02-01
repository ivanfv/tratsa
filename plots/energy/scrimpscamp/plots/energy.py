#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
# Importo la clase que lee los ficheros de estadisticas de gems
sys.path.append('/home/ifernandez/repos/tratsa-private/plots/energy/scrimpscamp/plots')
import numpy as np #Fundamental package for scientific computing with Python
import matplotlib.pyplot as plt 
import matplotlib
import parseStatsTransPrecision as pst
from subprocess import call #Para llamar a un comando del shell



files_scamp=[["../results/scamp/result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le8_lm23_*.csv", 8, 23],
       ["../results/scamp/result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le8_lm7_*.csv", 8, 7],
       ["../results/scamp/result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le5_lm10_*.csv", 5, 10]]


files_scrimp=[["../results/scrimp/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le8_lm23_*.csv", 8, 23],
              ["../results/scrimp/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le8_lm7_*.csv", 8, 7],
              ["../results/scrimp/result_audio-MPIII-SVD_w200_s1.00_t256_he8_hm23_le5_lm10_*.csv", 5, 10]]

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
#plt.rc('text', usetex=True) #Le digo que el procesador de textos sea latex (se pueden usar fórmulas)
#plt.rc('font', family='sans-serif') #Le digo que la fuente por defecto sea sans-serif

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

plt.rc('figure', autolayout=True) #Para que no se corten las etiquetas de los ejes (calcula el bouding box automaticamente)
motifColor = [63./256,127./256,191./256]
discordColor = [191./256,63./256,63./256]
otherColor = [63./256,191./256,63./256]

cmap = matplotlib.cm.get_cmap('CMRmap')
colors = (cmap(0),cmap(30),cmap(100),cmap(200))

#La energia de double w.r.t. double es 1
energy_scamp = [1]
energy_scrimp = [1]

energy_comp = []

for i in range(len(files_scamp)):
  if i==0:
    ######Primero saco los valores para todo a float SCAMP
    tmp = pst.StatsFile(files_scamp[i][0], False)
    # (INV,ADD,SUB,MUL,DIV,FMA,CMP)
    #considero INV,ADD,SUB como sumas (FlexFloat considera acumulaciones como INV)
    #MUL y DIV como multiplicaciones
    tmp.ParseFlexFloatBreak(files_scamp[i][1],files_scamp[i][2])
    inst = tmp.GetFlexFloatBreak()
    print("[SCAMP ] FlexFloat<%d,%d> %s" % (files_scamp[i][1], files_scamp[i][2], inst))
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    doubleE = fmas*eScalar[0][0] + muls*eScalar[1][0] + adds*eScalar[2][0] + comps*eScalar[3][0]
    floatE = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    floatESIMD = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]

    ######Primero saco los valores para todo a float SCRIMP
    tmp2 = pst.StatsFile(files_scrimp[i][0], False)
    # (INV,ADD,SUB,MUL,DIV,FMA,CMP)
    #considero INV,ADD,SUB como sumas (FlexFloat considera acumulaciones como INV)
    #MUL y DIV como multiplicaciones
    tmp2.ParseFlexFloatBreak(files_scrimp[i][1],files_scrimp[i][2])
    inst = tmp2.GetFlexFloatBreak()
    print("[SCRIMP] FlexFloat<%d,%d> %s" % (files_scrimp[i][1], files_scrimp[i][2], inst))
    print("\n")
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    doubleE_scrimp = fmas*eScalar[0][0] + muls*eScalar[1][0] + adds*eScalar[2][0] + comps*eScalar[3][0]
    floatE_scrimp = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    floatESIMD_scrimp = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]

    energy_comp.append(float(doubleE)/float(doubleE_scrimp))
    energy_comp.append(float(floatE)/float(floatE_scrimp))
    energy_comp.append(float(floatESIMD)/float(floatESIMD_scrimp))   
    energy_scrimp.append(float(floatE) / float(doubleE))
    energy_scrimp.append(float(floatESIMD) / float(doubleE))
    energy_scamp.append(float(floatE) / float(doubleE))
    energy_scamp.append(float(floatESIMD) / float(doubleE))

  else:
    ######Segundo los valores para mixed con 8,7 SCAMP
    tmp = pst.StatsFile(files_scamp[i][0], False)
    tmp.ParseFlexFloatBreak(8,23) #Primero los float
    inst = tmp.GetFlexFloatBreak()
    print("[SCAMP ] FlexFloat<%d,%d>-<%d/%d> %s " % (8,23,files_scamp[i][1], files_scamp[i][2], inst)),
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    mixESIMD = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]
    tmp.ParseFlexFloatBreak(files_scamp[i][1],files_scamp[i][2]) #luego los 8,7
    inst = tmp.GetFlexFloatBreak()
    print(inst)
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE += fmas*eScalar[0][i+1] + muls*eScalar[1][i+1] + adds*eScalar[2][i+1] + comps*eScalar[3][i+1]
    mixESIMD += fmas*eSIMD[0][i+1] + muls*eSIMD[1][i+1] + adds*eSIMD[2][i+1] + comps*eSIMD[3][i+1]

    ######Segundo los valores para mixed con 8,7 SCRIMP
    tmp = pst.StatsFile(files_scrimp[i][0], False)
    tmp.ParseFlexFloatBreak(8,23) #Primero los float
    inst = tmp.GetFlexFloatBreak()
    print("[SCRIMP] FlexFloat<%d,%d>-<%d,%d> %s " % (8,23,files_scrimp[i][1], files_scrimp[i][2], inst)),
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE_scrimp = fmas*eScalar[0][1] + muls*eScalar[1][1] + adds*eScalar[2][1] + comps*eScalar[3][1]
    mixESIMD_scrimp = fmas*eSIMD[0][1] + muls*eSIMD[1][1] + adds*eSIMD[2][1] + comps*eSIMD[3][1]
    tmp.ParseFlexFloatBreak(files_scrimp[i][1],files_scrimp[i][2]) #luego los 8,7
    inst = tmp.GetFlexFloatBreak()
    print(inst)
    adds = inst[0]+inst[1]+inst[2]
    muls = inst[3]+inst[4]
    fmas = inst[5]
    comps = inst[6]
    mixE_scrimp += fmas*eScalar[0][i+1] + muls*eScalar[1][i+1] + adds*eScalar[2][i+1] + comps*eScalar[3][i+1]
    mixESIMD_scrimp += fmas*eSIMD[0][i+1] + muls*eSIMD[1][i+1] + adds*eSIMD[2][i+1] + comps*eSIMD[3][i+1]
    print("\n")
    energy_comp.append(float(mixE)/float(mixE_scrimp))
    energy_comp.append(float(mixESIMD)/float(mixESIMD_scrimp))
    energy_scrimp.append(float(mixE_scrimp)/float(doubleE_scrimp))
    energy_scrimp.append(float(mixESIMD_scrimp)/float(doubleE_scrimp))
    energy_scamp.append(float(mixE)/float(doubleE))
    energy_scamp.append(float(mixESIMD)/float(doubleE))

print (energy_scrimp)
print (energy_scamp)
print (energy_comp)
fig, axs = plt.subplots(3)
x=range(len(energy_comp))
xshift=[i for i in x]
axs[0].bar(xshift, energy_scrimp, bw, color=colors[1], edgecolor=[i*0.6 for i in colors[1]])
axs[1].bar(xshift, energy_scamp, bw, color=colors[2], edgecolor=[i*0.6 for i in colors[2]])
axs[2].bar(xshift, energy_comp, bw, color=colors[3], edgecolor=[i*0.6 for i in colors[3]])

axs[0].set_title('a) \\texttt{TranSCRIMP}', loc='left')
axs[1].set_title('b) \\texttt{TranSCAMP}', loc='left')
axs[2].set_title('c) \\texttt{TranSCAMP} w.r.t \\texttt{TranSCRIMP}', loc='left')

axs[0].yaxis.set_ticks(np.arange(0, 1.25, 0.25))
axs[1].yaxis.set_ticks(np.arange(0, 1.25, 0.25))
axs[2].yaxis.set_ticks(np.arange(0, 1.25, 0.25))

axs[0].set_xticks([])
axs[1].set_xticks([])

axs[0].set_ylim(0, 1)
axs[1].set_ylim(0, 1)
axs[2].set_ylim(0, 1)

plt.xticks([i for i in x], ['Double', 'Single', 'Single\nSIMD', 'Mixed 8/7', 'Mixed 8/7\nSIMD','Mixed 5/10', 'Mixed 5/10\nSIMD'], rotation=32, fontsize=22)

axs[0].grid(True, axis='y')
axs[1].grid(True, axis='y')
axs[2].grid(True, axis='y')

fig.set_figheight(10)
fig.set_figwidth(10)

plt.savefig("./energy.pdf", format='pdf', bbox_inches='tight')
#plt.savefig("./accuracy.png", format='png', bbox_inches='tight')
#call("cp ../z_figs/SU.pdf ../../paper-TPDS/figs/", shell=True)
#plt.show()

