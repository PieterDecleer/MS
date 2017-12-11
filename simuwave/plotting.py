import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from numpy import sqrt, exp, pi
from scipy.constants import mu_0, epsilon_0, c, m_e, hbar, e

fi = open("file.txt","r")
A = [line.strip("\n").split("\t") for line in fi]
fi.close()

for row in A:
	del row[-1]

for a in range(len(A)):
	for b in range(len(A[0])):
		A[a][b] = A[a][b].replace(',','.')
C = np.array(A,float)
plt.imshow(C,aspect = 'auto',cmap = 'jet')
plt.colorbar()
#plt.clim(0,1e26)
plt.savefig("file.png")
plt.close()

fi = open("time_evolution_0.txt","r")
A = [line.strip("\n").split("\t") for line in fi]
fi.close()


for a in range(len(A)):
	for b in range(len(A[0])):
		A[a][b] = A[a][b].replace(',','.')
C = np.array(A,float)
dt = C[0,0]
dx = C[0,1]
Pr = C[1:,0]
Pi = C[1:,1]
frequencies = np.fft.fftfreq(len(Pr),dt)
Spectrum1 = np.fft.ifft(np.array(Pr+1j*Pi,dtype=complex))
plt.plot(hbar*2*pi*frequencies/e,np.abs(Spectrum1))
plt.xlim(0,10)
plt.ylim(0,8*10**12)
plt.xlabel('Electron energy in eV')
plt.ylabel('Amplitude (arbitrary units)')
plt.savefig('Spectrum_no_field.png')
plt.show()
plt.close()


for frame in range(300):
	fi = open("snapshot_%d.txt" %frame,"r")
	A = [line.strip("\n").split("\t") for line in fi]
	fi.close()
	for row in A:
		del row[-1]

	for a in range(len(A)):
		for b in range(len(A[0])):
			A[a][b] = A[a][b].replace(',','.')
	C = np.array(A,float)

	plt.imshow(C,aspect = 'auto',cmap = 'jet')
	plt.colorbar()
	plt.clim(0,4.5e25)
	plt.savefig("snapshot_%d.png" %frame)
	plt.close()




