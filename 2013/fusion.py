with open('fusion.txt', 'r') as f:
	fusion = [line.split('\t') for line in f.read().splitlines()]
	fuscol = [list(row) for row in zip(*fusion)]

import matplotlib.pyplot as plt
import numpy as np
from itertools import groupby

chromos = [line[2] for line in fusion]
ind = np.arange(22)
y = [len(list(group)) for key, group in groupby(sorted(chromos))]
chrlist = [line[0] for line in [list(group) for key, group in groupby(sorted(chromos))]]
width = 0.5

plt.bar(ind, y, width, color='#003366')
plt.ylabel('Frequency of Chromosome')
plt.title('Chromosomes')
plt.xticks(ind+width*1.0/2, chrlist)
plt.xlim(0,len(chrlist))
plt.show()
