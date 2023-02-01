#!/usr/bin/python
# -*- coding: utf-8 -*-
import re #Paquete para expresiones regulares
#import fnmatch #Paquete para Unix filename pattern matching
import os      #Paquete para acceder a comandos del OS
import sys
import math    #Paquete para pow y sqrt
import glob    #Para obtener una lista de archivos de un directorio

################################################################################
# class StatsFile()
#
# Clase que parsea un solo archivo
################################################################################
class StatsFile():
  MotifsAccuRe = re.compile("#Top-100 accuracy w\.r\.t double: Motifs float,float±10,flexfloat,flexfloat±10")
  DiscordsAccuRe = re.compile("#Top-100 accuracy w\.r\.t double: Discords float,float±10,flexfloat,flexfloat±10")
  ProfLenRe = re.compile("#Profile Length")
  SeriesRe = re.compile("#i,tseries,doubleMP,j,floatMP,j,errorf\(%\),flexfloatMP,j,errorff\(%\)")
  FlexFloatBreakRe = re.compile("#Flexfloat<(\d+),(\d+)>.*INV,ADD,SUB,MUL,DIV,FMA,CMP")
  
  def __init__(self, FileName, ParseAll):
    self.FileName = FileName
    self.ProfLen = 0
    self.Tseries = []
    self.MPdouble = []
    self.MPfloat = []
    self.MPff = []
    self.ErrorFloat = []
    self.ErrorFF = []
    self.MotifsAccu = []
    self.DiscordsAccu = []
    self.FlexFloatBreak = [] #almacena el desglose de instrucciones (INV,ADD,SUB,MUL,DIV,FMA,CMP) para una exp y man dados

    #El nombre del archivo puede venir con wildcards
    fileList = glob.glob(self.FileName)
    if len(fileList) > 1:
      raise Exception("El número de archivos es mayor que 1: %s"%fileList)
    elif len(fileList) < 1:
      self.FileName = ""
      self.MotifsAccu = [0,0,0,0]
      self.DiscordsAccu = [0,0,0,0]
      self.Tseries = [0]
      self.MPdouble = [0]
      self.MPfloat = [0]
      self.MPff = [0]
      self.ErrorFloat = [0]
      self.ErrorFF = [0]
      self.FlexFloatBreak = [0]
    else:
      self.FileName = fileList[0] #Me quedo con el nombre del archivo
      if ParseAll:
        self.ParseMotifsAccu()
        self.ParseDiscordsAccu()
        self.ParseProfLen()
        self.ParseSeries()
        self.ParseFlexFloatBreak(8,23)

  def ParseMotifsAccu(self):
    f = open(self.FileName, 'r') #Por defecto es 'r'
    for line in f:
      m = self.MotifsAccuRe.match(line)
      if m:
        break
    for line in f:
      self.MotifsAccu = [float(i) for i in line.split(",")]
      break
    f.close()

  def ParseDiscordsAccu(self):
    f = open(self.FileName, 'r') #Por defecto es 'r'
    for line in f:
      m = self.DiscordsAccuRe.match(line)
      if m:
        break
    for line in f:
      self.DiscordsAccu = [float(i) for i in line.split(",")]
      break
    f.close()
    
  def ParseProfLen(self):
    f = open(self.FileName, 'r') #Por defecto es 'r'
    for line in f:
      m = self.ProfLenRe.match(line)
      if m:
        break
    for line in f:
      self.ProfLen = int(line)
      break
    f.close()

  def ParseSeries(self):
    f = open(self.FileName, 'r') #Por defecto es 'r'
    for line in f:
      m = self.SeriesRe.match(line)
      if m:
        break
    i = 0
    for line in f:
      if i >= self.ProfLen:
        break
      tmp = line.split(",")
      self.Tseries.append(float(tmp[1]))
      self.MPdouble.append([float(tmp[2]), int(tmp[3])])
      self.MPfloat.append([float(tmp[4]), int(tmp[5])])
      self.ErrorFloat.append(float(tmp[6]))
      self.MPff.append([float(tmp[7]), int(tmp[8])])
      self.ErrorFF.append(float(tmp[9]))
      i += 1
    f.close()
    
  def ParseFlexFloatBreak(self,exp,man):
    f = open(self.FileName, 'r') #Por defecto es 'r'
    for line in f:
      m = self.FlexFloatBreakRe.match(line)
      if m:
        if int(m.group(m.lastindex - 1)) == exp and int(m.group(m.lastindex)) == man:
          break
    for line in f:
      tmp = line.split(",")
      self.FlexFloatBreak = [int(i) for i in tmp]
      break
    f.close()

  def GetMotifsAccu(self):
    return self.MotifsAccu
  def GetDiscordsAccu(self):
    return self.DiscordsAccu
  def GetProfLen(self):
    return self.ProfLen
  def GetTSeries(self):
    return self.Tseries
  def GetDoubleMP(self):
    return self.MPdouble
  def GetFloatMP(self):
    return self.MPfloat
  def GetFlexMP(self):
    return self.MPff
  def GetErrorFloat(self):
    return self.ErrorFloat
  def GetErrorFF(self):
    return self.ErrorFF
  def GetFlexFloatBreak(self):
    return self.FlexFloatBreak


####################################################################
# Pruebas
####################################################################
if __name__ == "__main__":
   a=StatsFile("/home/quislant/Ricardo/Research/PostTesis/TransPrecision/scrimp_ff/results/result_audio-MPIII-SVD_w200_s1.00_t256_e5_m10_232255.csv", True)
   tmp = a.GetMotifsAccu()
   print("Motifs Accuracy: %.2f,%.2f,%.2f,%.2f" % (tmp[0], tmp[1], tmp[2], tmp[3]))
   tmp = a.GetDiscordsAccu()
   print("Discords Accuracy: %.2f,%.2f,%.2f,%.2f" % (tmp[0], tmp[1], tmp[2], tmp[3]))

