import numpy as np
import pylab
import matplotlib.pyplot as plt
import sys

original_data = pylab.loadtxt(str(sys.argv[1]));
mat_profile   = pylab.loadtxt(str(sys.argv[2]));

ind_org_data = np.arange(0, len(original_data), 1);
ind_mat_prof = np.arange(0, len(mat_profile), 1);

print("Original data length: " + str(len(original_data)))

plt.subplot(2, 1, 1)
plt.plot(ind_org_data, original_data, 'C3')

plt.xlabel('Sample')
plt.ylabel('Value')
plt.title(sys.argv[1])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(ind_mat_prof, mat_profile, 'C3')
plt.xlabel('Sample')
plt.ylabel('Value')
plt.title(sys.argv[2])
plt.grid(True)

plt.show()
