# REMARKS: 
# (1) the stl-file and mesh objects use millimeter units.
# (2) the constraint points should also take into account the curvature of the objects
#       --> triangular-mesh splitting along coordinate planes 
# (3) the number of samples per skin depth should be applied to the boundaries of the conductor,
#       allowing for a coarser nonuniform discretization inside the conductor (CST does this too)
#       --> for a grading ratio r, the first skin depth of a conductor with n samples per skin depth at its boundary 
#           is resolved by log(n)/log(r) samples
#       --> e.g. r=1.5 and n=20 resolves the skin depth by 7.4 nonuniform samples     
# (4) nonuniformity and material interfaces should not coincide to guarantee sufficient accuracy
#       => constraint points pose a havier burden than in Berens' paper as the steps left and right of it need to be identical
# (5) Newton-Raphson achieves quadratic convergence if no zero derivative occurs

from math import ceil, floor, sqrt, log, exp, pi
import numpy as np


c0 = 299792458.0      # vacuum speed of light 
Z0 = 376.730313461    # vacuum wave impedance
tol = 1e-9            # precision of the spatial discretization (mm)


#--------------------------------------------------------------------------------------------------------------------------------------------------


def newton(r0, n, c, higher=False):

	'''
	Newton-Raphson root-finding algorithm for the truncated geometric series with variable r (grading ratio).

	Input:
	   r0 = grading ratio (variable, guess) 
	   n  = number of cells (fixed)
           c  = L/d = interval length divided by the smallest step (fixed)
	   higher = True if we are looking for r>r0, which is usually not the case 
	
	Output:
	   r  = grading ratio (root)
	   warning = True if something went wrong  
	'''

	warning = False

	if n<=1: return 0.0  , warning                  # resolve trivial cases
	if n==2: return c-1.0, warning

	Niter = 50

	if higher:
		a = exp(log(2)/Niter)                   # the guess for r ranges from r0 to 2*r0
	else:
		a = exp(-log(r0)/Niter)                 # the guess for r ranges from r0 to 1
	r = r0/a

	converged = False 
	i = 0           

	while not converged and i<Niter:

		r = r0*a**i
		j = 0   

		while not converged and j<Niter:
			
			f  = (1.0-r**n)-(1.0-r)*c       # condition: the sum of all steps should equal L (mostly, r=1 is a fake root)
			df = c-n*r**(n-1)               # derivative (c=1, n=1 --> df=0)
			
			if abs(df)<1e-10: break			
			
			r_err = f/df
			r     = r-r_err
			
			converged = abs(r_err)<1e-10
			
			j = j+1      

		i = i+1   

	if not converged: 
		print('WARNING: Newton-Raphson method did not converge.'); warning = True
	elif abs(r-1.0)<1e-10:
		print('WARNING: Newton-Raphson method converged to the possibly fake root 1.'); warning = True
		
	return r, warning                                       


#-------------------------------------------------------------------------------------------------------------------------------------------------- 


def discr1D(cp, d_seg, d_cp, r):

	'''
	segment discretization

	Input:
	   cp    = array of constraint points, delimiting segments
	   d_seg = array of desired step sizes in each segment (bulk discretization)
	   d_cp  = array of desired step sizes at each constraint point (edge discretization)
           r     = grading ratio 
	
	Output:
	   steps = the final discretization along a specific axis (list of steps)
	'''

	steps = [0]*len(d_seg)
	for i in np.argsort(d_seg):
		s=[]
		da = d_cp[i]                               # left step          
		db = d_cp[i+1]                             # right step
		dmax = d_seg[i]                            # maximum step  		
		ra = r                                     # left grading ratio
		rb = r		                           # right grading ratio 
		Na = int(floor(1.0+log(dmax/da)/log(r)))   # dmax/r<= da*r^(Na-1) <dmax  
		Nb = int(floor(1.0+log(dmax/db)/log(r)))   # dmax/r<= db*r^(Nb-1) <dmax 
		L  = cp[i+1]-cp[i]                         # segment length 		
		La = da*(1.0-r**Na)/(1.0-r)                # total length of the Na steps
		Lb = db*(1.0-r**Nb)/(1.0-r) 		   # total length of the Nb steps
		Lm = L-La-Lb    	                   # m=middle	
		if tol<Lm:
			Nm = int(floor(Lm/dmax))
			if Na==Nb==1 and Nm>0:
				dmax = Lm/Nm 		   # maximum step violation 	
			else:         
				if Na>Nb:
					Na = Na+1
					La = La+Lm-Nm*dmax
					ra, w = newton(r, Na, La/da)
					if w:
						Na = Na-1
						ra, _ = newton(r, Na, La/da, higher=True)     # grading ratio violation 
				else:
					Nb = Nb+1
					Lb = Lb+Lm-Nm*dmax
					rb, w = newton(r, Nb, Lb/db)             
					if w:
						Nb = Nb-1					
						rb, _ = newton(r, Nb, Lb/db, higher=True)     # grading ratio violation
           		s.extend([da*ra**n for n in range(Na)])
			s.extend([dmax]*Nm)
			s.extend([db*rb**n for n in range(Nb-1,-1,-1)])
		else:
			x = (L+(db-da)/(r-1.0))/2.0        # theoretical intersection
			if x<da:                                       
				b = [db]
				while sum(b)<L: b.append(r*b[-1])
				Nb = len(b)
				rb, _ = newton(r, Nb, L/db)                 
				s.extend([db*rb**n for n in range(Nb-1,-1,-1)])
			elif x>L-db:                    
				a = [da]
				while sum(a)<L: a.append(r*a[-1])
				Na = len(a)
				ra, _ = newton(r, Na, L/da)
				s.extend([da*ra**n for n in range(Na)])					
			else:
				a = [da]
				while sum(a)<x: a.append(r*a[-1])
				Na = len(a)
				ra, _ = newton(r, Na, x/da)
				s.extend([da*ra**n for n in range(Na)])
				b = [db]
				while sum(b)<(L-x): b.append(r*b[-1])
				Nb = len(b)
				rb, _ = newton(r, Nb, (L-x)/db)
				s.extend([db*rb**n for n in range(Nb-1,-1,-1)])

		steps[i] = [e for e in s if e>tol]	
	
		if not all([1.0/r<steps[i][j+1]/steps[i][j]<r for j in range(len(steps[i])-1)]): 
			print('WARNING: grading ratio violation detected.')	

		if abs(sum(steps[i])-L)>tol: 
			print('WARNING: segment is not properly discretized.')

	steps = [w for v in steps for w in v] 

	if abs(sum(steps)-(cp[-1]-cp[0]))>tol: 
		print('WARNING: axis is not properly discretized.')

	return steps


#--------------------------------------------------------------------------------------------------------------------------------------------------


def rayCasting(mesh, points):
	'''
	Check if a mesh contains a set of points using ray tests, 
	also called ray casting or parity of intersection count.

	Input:
	   mesh   = a single Trimesh object
	   points = (n,3) points in space

	Output:
	   contains = (n) boolean array, whether point is inside mesh or not


	Remarks:

	1. Depending on whether or not you installed Anaconda, the ray-triangle intersections are done by:
	    - slow rtree/numpy queries (standard)
	    - fast pyembree queries (Anaconda)
 	
	   The r-tree algorithm groups nearby triangles inside rectangular (r) minimum bounding boxes organized in a tree data 
	   structure, also called a bounding volume hierarchy. Note that determining whether or not a point is inside a rectangular 
           volume is a very cheap computation and does not require any advanced intersection algorithm. Hence, the effective number 
           of ray-triangle intersection operations is strongly reduced by the r-tree.    

	   Pyembree, which uses Intel's embree code, is more than 50 times faster, but is imcompatible with Windows. 

	2. Some measures are taken to increase the robustness:
	    - rtree/numpy uses one set of rays with a random direction
	    - pyembree uses the majority vote of the outcome of three sets of axis-aligned rays 

	3. Further acceleration could be achieved by exploiting the fact that axis-aligned rays contain information about the parity of
	   intersection for mulitple points in the tensor-product grid. The number of ray queries could be reduced from a 3D problem to 
	   a 2D problem. This only poses a gain in overall performance if the intersection count along the third dimension is done with 
	   low-level loops (and not in Python).  
	'''	      

	contains = np.zeros(len(points), dtype=np.bool)
	ray_origins = np.asanyarray(points)

	try:
		# Embree is more robust for axis-aligned rays
		from pyembree import rtcore_scene  
		direction = [1,0,0]                                                    
		n = 3                                                                                                             
	except ImportError: 
		# r-tree/numpy is more robust for random directions 
		direction = [0.483215846215448752, 0.912215483121, 0.154831486512]    
		n = 1

	c = [None]*n

	for i in range(n):
		ray_directions = np.tile(np.roll(direction,i), (len(points), 1))  
		locations, index_ray = mesh.ray.intersects_location(ray_origins, ray_directions)

		if len(locations)==0: return contains 

		hit_count = np.bincount(index_ray, minlength=len(points))
		c[i] = np.mod(hit_count, 2)==1
		del locations, index_ray, hit_count

	if n==1:
		contains = c[i]
	elif n==3: 
		# majority vote	
		contains = np.logical_or(np.logical_and(c[0],c[1]), np.logical_and(c[2], np.logical_or(c[0],c[1])))

	return contains


#--------------------------------------------------------------------------------------------------------------------------------------------------


def materialMapping(g):

	'''
	Assign integer labels to the different materials (matlut) and give every Yee cell a material label (matmap).	

	Input:
	   g    = Grid object
	
	Output:
           Following attributes are added to g:
	      matmap = column-major list containing unsigned integer labels identifying the different materials assigned to the centers of the Yee cells
	      matlut = list containing material tuples (mur,epsr,sigma) of each mesh object, serving as a look-up table for matmap 	
	'''	
         
	# set up the matlut    
	background = (1.0,1.0,0.0)                                                    # the background medium is chosen to be vacuum	
	mesh = filter(lambda m: (m.mur,m.epsr,m.sigma)!=background, g.mesh)           # remove the mesh objects which are filled with the background medium 
	matlut = [(m.mur,m.epsr,m.sigma) for m in mesh]	  
	matlut.append(background)	                                                                       
	matlut = sorted(set(matlut))                                                  # material tuples should occur only once in the look-up table
	
	# initialize matmap (uniform background medium)
	dual = g.dual_edges()
	x, y, z = np.meshgrid(dual[0], dual[1], dual[2], indexing='ij')	              # Yee-cell centers
	N = (dual[0].size, dual[1].size, dual[2].size)		
	points = np.vstack((np.reshape(x,-1), np.reshape(y,-1), np.reshape(z,-1))).T  # column-major ordering (z-adjacency)	
	matmap = matlut.index(background)*np.ones(N)						              
	
	# ray casting
	for m in mesh:
		cont = np.zeros((points.shape[0],1))	                                 
		mask = np.logical_and(np.all(m.bounds[0]<points,axis=1),np.all(points<m.bounds[1],axis=1))           	   	
		cont[mask,0] = rayCasting(m, points[mask,:]).astype(np.int)	      # AABB test: only consider points inside the AABB 		
		matmap += matlut.index((m.mur,m.epsr,m.sigma))*np.reshape(cont,N)                             
	
	# add PML layers to matmap (materials touching the PML are extended inside the PML) 
	if hasattr(g,'Npml'):	
		matmap = np.concatenate( ( np.tile(matmap[[0],:,:],(g.Npml[0],1,1)), matmap, np.tile(matmap[[-1],:,:],(g.Npml[1],1,1)) ), axis=0)
		matmap = np.concatenate( ( np.tile(matmap[:,[0],:],(1,g.Npml[2],1)), matmap, np.tile(matmap[:,[-1],:],(1,g.Npml[3],1)) ), axis=1)	
		matmap = np.concatenate( ( np.tile(matmap[:,:,[0]],(1,1,g.Npml[4])), matmap, np.tile(matmap[:,:,[-1]],(1,1,g.Npml[5])) ), axis=2)	

	# convert matmap to an integer list with column-major ordering		
	matmap = map(int,(matmap.ravel(order='C')).tolist())

	# pass-by-object 
	g.matmap = matmap
	g.matlut = matlut	


#--------------------------------------------------------------------------------------------------------------------------------------------------


def gridGeneration(mg, sg, mesh, gr, dist):
	
	'''
	Input:
	  mg   = MainGrid object with attributes Npml, Kmax and Amax
	  sg   = list of SubGrids objects with attributes pt1, pt2 and maxstep		
	  mesh = list of Trimesh objects containing material attributes
	  dist = list of distances between the bounding box of the CAD mesh and the FDTD simulation box for each of the 6 faces [mm]  
	  gr   = grading ratio (the maximum size ratio of two adjacent cells)

	Output:
	  Following attributes are added to mg and sg: 
	     edges  = 3x1 list containing 1d numpy arrays with locations of the primary and secondary-grid edges		  
	     matmap = 3d numpy array containing unsigned integer labels identifying the different materials assigned to the centers of the Yee cells
	     matlut = list containing material tuples (mur,epsr,sigma) of each mesh object, serving as a look-up table for matmap		
	'''

	# assign the mesh objects to the correct grid (keep in mind that meshes could be in multiple grids)
	for m in mesh:
		assigned=False
		for s in sg:
			if s.contains(m.bounds[0]) and s.contains(m.bounds[1]):    # only full overlap is supported at the moment
				s.mesh.append(m)
				assigned=True	 
		if not assigned:
			mg.mesh.append(m)


	print('Meshing of main grid ...')

	# discretization for each dimension (i={0,1,2} corresponds to {x,y,z})
	for i in range(3):

		# constraint points determined by the axis-aligned bounding box (AABB) of each object  
		cp = []                                          
		for m in mg.mesh:	
			cp.append(m.bounds[0][i])                          
			cp.append(m.bounds[1][i])
		for s in sg:
			cp.append(s.bounds[0][i])                          
			cp.append(s.bounds[1][i]) 
			if not s.impl[i]:			
				for m in s.mesh:
					cp.append(m.bounds[0][i])
	   				cp.append(m.bounds[1][i])
		cp = np.sort(cp)                                                                          
		cp = np.insert(cp,0,cp[0]-dist[2*i])                     # add exterior boundaries  
		cp = np.append(cp,cp[-1]+dist[2*i+1])    		
		cp = cp[np.insert(np.abs(np.diff(cp))>tol, 0, 1)]        # remove degenerate points                 
                                                          	
		# assign a maximum step to each segment in between two constraint points		
		d_seg = np.empty(len(cp)-1)	                         			
		for j in range(len(cp)-1):
			L = cp[j+1]-cp[j]                                # segment length
			a = (cp[j]+cp[j+1])/2.0	                         # segment center
			d_seg[j] = min(mg.maxstep[i], L)			                
			for m in mg.mesh:                                                                   
				if m.bounds[0][i]<a<m.bounds[1][i]:				
					d_seg[j] = min(d_seg[j], m.maxstep[i])	
					if not m.pec:
						d_seg[j] = min(d_seg[j], mg.maxstep[i]/sqrt(m.epsr*m.mur))														
			for s in sg: 
				if not s.impl[i]:
					for m in s.mesh:
						if m.bounds[0][i]<a<m.bounds[1][i]:					
							d_seg[j] = min([d_seg[j], m.maxstep[i], s.maxstep[i]])
							if not m.pec:
								d_seg[j] = min(d_seg[j], mg.maxstep[i]/sqrt(m.epsr*m.mur)) 
			# segments with few steps have less flexibility and are likely to cause grading ratio violations	                  
			if L/d_seg[j]<5: d_seg[j] = L/ceil(L/d_seg[j])         


		# assign a step to the constraint points 
		d_cp = np.empty(len(cp))
		d_cp[0]  = d_seg[0]
		d_cp[-1] = d_seg[-1]	
		for j in range(len(cp)-2):
			d_cp[j+1] = min(d_seg[j], d_seg[j+1])			
			for m in mg.mesh:
                           	if m.sigma and (abs(cp[j+1]-m.bounds[0][i])<tol  or abs(cp[j+1]-m.bounds[1][i])<tol):
					# approximately 20 samples per skin depth, mm units				
					d_cp[j+1] = min(d_cp[j+1], max(1.0/sqrt(mg.fmax*m.sigma), tol))       

		# reduce spurious reflections by matching the normal cell size at both sides of the maingrid-subgrid interface 
		

		# nonuniform discretization of all segments along one dimension                                                           	
		steps = discr1D(cp, d_seg, d_cp, gr)

		# store locations of primary edges 	
		mg.edges[i] = np.insert(cp[0]+np.cumsum(np.asarray(steps)), 0, cp[0])                             
	                               										                                                                                                 
	# ray casting
	print('Ray casting of main grid ...')
	materialMapping(mg)
	


	for k,s in enumerate(sg):
 
		print('Meshing of subgrid ' + str(k) + ' ...')

		# discretization for each dimension (i={0,1,2} corresponds to {x,y,z})
		for i in range(3):
			
			if not s.impl[i]:

				# copy the discretization of the main grid 
				s.edges[i] = mg.edges[i][np.logical_and(s.bounds[0][i]-tol<mg.edges[i], mg.edges[i]<s.bounds[1][i]+tol)] 
 
			else:
				# constraint points determined by the axis-aligned bounding box (AABB) of each object  
				cp = []       
				cp.append(s.bounds[0][i])
				cp.append(s.bounds[1][i])                                     
				for m in s.mesh:	  
					cp.append(m.bounds[0][i])                          
					cp.append(m.bounds[1][i])  					
				cp = np.sort(cp)                                                                          
				cp = cp[np.insert(np.abs(np.diff(cp))>tol, 0, 1)]                                                            
						                          	
				# assign a maximum step to each segment in between two constraint points		
				d_seg = np.empty(len(cp)-1)		
				for j in range(len(cp)-1):
					L = cp[j+1]-cp[j]                                
					a = (cp[j]+cp[j+1])/2.0					
					d_seg[j] = min(mg.maxstep[i], L)			                
					for m in s.mesh:
						if m.bounds[0][i]<a<m.bounds[1][i]:					
							d_seg[j] = min([d_seg[j], m.maxstep[i], s.maxstep[i]])
							if not m.pec:
								d_seg[j] = min(d_seg[j], mg.maxstep[i]/sqrt(m.epsr*m.mur))     	                  
					if L/d_seg[j]<5: d_seg[j] = L/ceil(L/d_seg[j])                                   
				
				# assign a step to the constraint points 
				d_cp[0]  = d_seg[0]
				d_cp[-1] = d_seg[-1]	
				for j in range(len(cp)-2):
					d_cp[j+1] = min(d_seg[j], d_seg[j+1])	
					for m in s.mesh:
				           	if m.sigma and (abs(cp[j+1]-m.bounds[0][i])<tol  or abs(cp[j+1]-m.bounds[1][i])<tol):				
							d_cp[j+1] = min(d_cp[j+1], max(1.0/sqrt(mg.fmax*m.sigma), tol))             


				# nonuniform discretization of all segments along one dimension                                                      	
				steps = discr1D(cp, d_seg, d_cp, gr)

				# store locations of primary edges 		
				s.edges[i] = np.insert(cp[0]+np.cumsum(np.asarray(steps)), 0, cp[0])                     

		# ray casting
		print('Ray casting of subgrid ' + str(k) + ' ...')
		materialMapping(s)

