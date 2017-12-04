import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
for frame in range(10):
	fi = open("snapshot_%d.txt" %frame,"r")
	A = [line.strip("\n").split("\t") for line in fi]
	fi.close()
	print(len(A))
	print(len(A[0]))
	for row in A:
		del row[-1]
	print(len(A))
	print(len(A[0]))
	for a in range(len(A)):
		for b in range(len(A[0])):
			A[a][b] = A[a][b].replace(',','.')
	C = np.array(A,float)

	plt.imshow(C,aspect = 'auto',cmap = 'jet')
	plt.colorbar()
	#plt.clim(0,5e25)
	plt.savefig("snapshot_%d.png" %frame)
	plt.close()




