#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# data to plot
n_groups = 6
exec_times = (10387.043, 5889.1089, 2376.0204, 990.7457, 428.900575, 18.79785)

							
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

index = np.arange(n_groups)
bw = 0.6
cmap = matplotlib.cm.get_cmap('CMRmap')
colors = (cmap(0),cmap(30),cmap(100),cmap(200))

# create plot
#fig, ax = plt.subplots(nrows=1, ncols=1)
#
# create plot
#ax =plt.figure(figsize=(8,4))

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

f.set_size_inches(8, 2.5,forward=True)


ax1.bar(index[0:len(index) -1], exec_times[0:len(index) -1], bw, color='#0071C5', edgecolor='#0071C5', label='Xeon Phi KNL')
ax1.bar(index[len(index) -1], exec_times[len(exec_times) -1], bw, color='#e52719', edgecolor='#0071C5', label='Alveo U50')

ax2.bar(index[0:len(index) -1], exec_times[0:len(index) -1], bw, color='#0071C5', edgecolor='#0071C5', label='Xeon Phi KNL')
ax2.bar(index[len(index) -1], exec_times[len(exec_times) -1], bw, color='#e52719', edgecolor='#e52719', label='Alveo U50')
#plt.xlabel('Time Series')

#ax1.text(index[0] - 0.1, 1800, '10387', fontsize=18, rotation=90, color='w')


# zoom-in / limit the view to different portions of the data
ax1.set_ylim(1100, 11000)  # outliers only
ax2.set_ylim(0, 1000)  # most of the data

# hide the spines between ax and ax2
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .01  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


#plt.ylabel('Execution Time (s)')
ax2.text(index[0] - 0.1, 70, '10387', fontsize=18, rotation=90, color='w')
ax2.text(index[1] - 0.1, 70, '5889', fontsize=18, rotation=90, color='w')
ax2.text(index[2] - 0.1, 70, '2376', fontsize=18, rotation=90, color='w')
ax2.text(index[3] - 0.1, 70, '990', fontsize=18, rotation=90, color='w')
ax2.text(index[4] - 0.1, 70, '428', fontsize=18, rotation=90, color='w')
ax2.text(index[5] - 0.1, 70, '18.8', fontsize=18, rotation=90, color='k')


ax1.grid(axis='y')
ax2.grid(axis='y')
plt.xticks(index, ('Core2Quad\nQ9400', 'i5\n4570', 'i7\n8700', 'XeonPhi\n7210', 'XeonGold\n6154', 'Alveo\nU50'))
#plt.xlim((-3*bw, n_groups - 2*bw))
#plt.legend(fontsize=16, ncol=2, columnspacing=1, loc='upper left',frameon=False)
#plt.yticks(np.arange(0, 16, step=5))
plt.yticks((0,500,1000))

f.text(0.01, 0.5, 'Execution Time (s)', va='center', rotation='vertical')
#plt.ylim((0,2500))



#plt.tight_layout()
plt.savefig("Speedup.pdf", format='pdf', bbox_inches='tight')#, bbox_extra_artist=lgd)
#plt.show()
