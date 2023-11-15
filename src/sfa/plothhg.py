import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt('hhg_At.out')
frequencies = data[:, 0]
real_part = data[:, 1]
imag_part = data[:, 2]

frequencies = frequencies/.057
norm = np.sqrt((real_part**2 + imag_part**2))
#print the dimensions of data
print(np.shape(data))

#plot data on log scale

fig = plt.figure()
plt.semilogy(frequencies, norm)
plt.xlabel('Frequency')
plt.ylabel('HHG')
#plt.xlim([0, 20])
fig.savefig('hhg_At.png')

data = np.loadtxt('dipole.out')
time = data[:, 0]
real_part = data[:, 1]
imag_part = data[:, 2]

norm = np.sqrt((real_part**2 + imag_part**2))

fig = plt.figure()
plt.plot(time, norm)
plt.xlabel('Time')
plt.ylabel('Dipole')

fig.savefig('dipole.png')