import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("energy.dat")
#data = np.loadtxt("initialVelocities.out")
plt.hist(data[:,0], bins=50)
#for i in range(0, 1):
#    plt.hist(data[:, i], bins=20)
#plt.hist(data[:, 1], bins=20)
#plt.hist(data[:, 2], bins=20)

plt.show()
