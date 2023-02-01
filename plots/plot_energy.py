import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as mticker


# data to plot
n_groups = 5


total = (1,0.3931,0.40035,0.335802,0.33011)


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=26)

# create plot
fig, ax = plt.subplots(nrows=1, ncols=1)
index = np.arange(n_groups)
bar_width = 0.5
opacity = 1



rects2 = plt.bar(index, total, bar_width,
alpha=opacity,
color='b',
#hatch='OO',
edgecolor = "black",label='Memory')



#plt.xlabel('Configuration')
plt.ylabel('Energy consumption',fontsize=26 )
#plt.title('Normalized energy consumption')
plt.xticks(index , ('bin64', 'bin32', 'SIMD\n bin32','SIMD \nbin32/bin16', 'SIMD \nbin32/bin16alt'), rotation=0,fontsize=20 )


plt.yticks(np.arange(0, 1, step=0.5),fontsize=26)


plt.yscale("log")
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())


fig.set_size_inches(10, 7,forward=True)
#plt.savefig('energy_instant.pdf')
plt.tight_layout()

plt.show()
#plt.savefig('gpc_speedup.pdf')
