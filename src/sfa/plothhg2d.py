import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt('hhg.out')
frequencies = data[:, 0]
real_partx = data[:, 1]
imag_partx = data[:, 2]
real_party = data[:, 3]
imag_party = data[:, 4]

xdat = np.sqrt((real_partx**2 + imag_partx**2))
ydat = np.sqrt((real_party**2 + imag_party**2))

frequencies = frequencies/.057
#plot data on log scale

print(np.shape(xdat),np.shape(ydat))

fig = plt.figure()
plt.semilogy(frequencies, xdat, label='x')
plt.semilogy(frequencies, ydat, label='y')
plt.legend()
#plt.xlim([0, 20])
plt.xlabel('Frequency')
plt.ylabel('HHG')
fig.savefig('hhg.png')

data = np.loadtxt('dipole2d.out')
time = data[:, 0]
real_partx = data[:, 1]
imag_partx = data[:, 2]
real_party = data[:, 3]
imag_party = data[:, 4]

xdat = np.sqrt((real_partx**2 + imag_partx**2))
ydat = np.sqrt((real_party**2 + imag_party**2))

fig = plt.figure()
plt.plot(time, xdat, label='x')
plt.plot(time, ydat, label='y')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Dipole')

fig.savefig('dipole.png')
