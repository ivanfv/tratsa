import numpy as np
import pylab
import sys

inputFile = open(str(sys.argv[1]), 'r')
    
num_list = [float(num) for num in inputFile.read().split()]

inputFile.close()
print("Time Series: " + str(sys.argv[1]))
print("Time Series length: " +  str(len(num_list)))
print("Time Series MAX val: " + str(max(num_list)))
print("Time Series MIN val: " + str(min(num_list)))
