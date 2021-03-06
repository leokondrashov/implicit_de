#!/usr/bin/python

import sys
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(open(sys.argv[1], "rb"), delimiter=",").T
plt.scatter(data[1], data[2], s=0.1, c=data[0])
plt.xlabel('theta')
plt.ylabel('phi')
plt.grid()
plt.savefig(sys.argv[1] + '.png', dpi=300)
