# REMARKS: 
# (1) Change the extension to '.pyw' to ensure compatibility with Mac and Windows
# (2) The layout manager takes care of the parent-child relation of the widgets: 
#       if the main window is closed, all widgets are closed as well without error.
# (3) Mayavi is notorious for memory leaks: 
#       - 'stop()', 'remove()', 'mlab.clf()', 'gc.collect()', ... do not correctly free the occupied memory
#       - 'mlab.close()' seems to correclty free the memory, but requires to open a new window for the next plot 
#       - 'mlab.options.offscreen=True' seems to be free of memory leaks 
# (4) ADHIE: 
#       - the step thresholds are solely based on the primary edges
#       - the two dual steps overlapping with an implicitized primary step are implicitized as well 
# (5) UNITS:
#	- All distances used by the gui are in millimeters (mm). The FDTD C-module converts this to meters.
#       - All time quantities used by gui are in nanoseconds (ns).
#       - All frequencies used by the gui are in gigaHertz (GHz).
#       - The gui uses the physical electric field, whereas the FDTD C-module uses the rescaled electric field E/Z0. 
#       - The global constants and the material constants have SI units.     
# (6) You can define only one source (in the main grid and/or one or more of the subgrids)  
# (7) One sensor can only produce one scalar quantity for each time step (so no cuboid field sensors are possible, 
#     unless you insert them one by one). 
# (8) Everything that happens in the subgrids, also happens in the underlying main grid, but is shielded. 
# (9) Be careful with list concatination because it sometimes copies the reference to the value and not the value itself. 
#      For example, try 'a=[[0]*3]*2; a[0][0]=1' 
#      Similarly, extend() cannot act on seperate list elements but instead extends every element in the list.    
# (10) Python function arguments are "passed by object" => giving the mayavi scene to a function allows to change the scene 
# (11) Maybe, it is better to create a class Trimat that extends the Trimesh class by a set of material constants.
#	This could be done by inheritance or by composition. Even a third option is to modify the original Trimesh class.
#       Now, material constants are dynamically added to the Trimesh objects, which is an ugly work-around.   
# (12) Subgrid material traverse is not yet implemented.        


from math import pi, sqrt, exp, cos, sin, floor, log10
import numpy as np
import matplotlib.pyplot as plt
import h5py          
from time import time
from random import random

from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor  

import trimesh
from mesher import gridGeneration
from cython_wrapper import run_fdtd_python 



###################################################################################################################################################
####### GLOBAL CONSTANTS AND FUNCTIONS ############################################################################################################
###################################################################################################################################################


global Z0, c0, tol, inf 
Z0 = 376.730313461    # (ohms)
c0 = 299792458.0      # (m/s)
tol = 1e-9            # tolerance of the spatial discretization (mm)


def deleteAttr(obj,attr):
	''' Free the memory of a list of attributes attr belonging to object obj. '''
	for a in attr:
		if hasattr(obj,a):
			setattr(obj,a,None) 
			delattr(obj,a)   	   



def le2float(le,a=0):
	''' Convert the text from a Qt line editor to a floating point number. 
	    If the text is not a valid number, a default value (a) is returned. '''
	try:    return float(le.text())
	except: return a



def closestPow10(x):
	''' Computes the closest power of 10. '''
	if abs(x)>1e-100: return 10**floor(log10(x))	
	else:             return 1.0                   # no rescaling to avoid division by zero



def my_fft(s, dt):
	''' Compute the amplitude of the fast Fourier transform of the time signal s which is uniformly sampled with step dt.
            The time signal s is padded with zeros to bring the number of samples to the closest power of 2.  '''
	n  = 2**(s.shape[0]-1).bit_length() 		
	S  = dt*np.abs(np.fft.fft(s,n)[:n/2])				
	f  = np.linspace(0, 0.5/dt, num=n/2)	
	return S, f 



def my_plot(t, s, f, S, t_min=0, t_max=100, f_min=0, f_max=100, t_label='amplitude', f_label='amplitude'): 
	''' Plot time-domain data s(t) and frequency-domain data S(f) above each other. '''					
	bt = (t_min<=t)*(t<=t_max)   
	bf = (f_min< f)*(f<=f_max)	 # b=boolean		
	plt.subplot(211)		
	plt.plot(t[bt],s[bt])
	plt.xlabel('time (ns)')
	plt.ylabel(t_label)
	plt.grid(True)
	plt.subplot(212)
	plt.plot(f[bf],S[bf])	
	plt.xlabel('frequency (GHz)')
	plt.ylabel(f_label)
	plt.grid(True)
	plt.show()



def showVoxels(s, g):
	''' Plot irregular (cell-conform) cuboids to visualize the voxelized objects in the grid g.'''	
	vox  = [] 	                                                       # voxel object consisting of Mayavi points3d plots	
	prim = g.steps                                                         # primary steps including PMLs
	dual = g.dual_edges()                                                  # dual edges excluding PMLs
	Npml = [0]*6 if not hasattr(g,'Npml') else g.Npml                      # number of PML layers     
	N = tuple([len(prim[i]) for i in range(3)])                            # number of cells in each dimension including PMLs
	m = np.reshape(g.matmap, N)                                            
	m = m[Npml[0]:N[0]-Npml[1],                                            # material map 3d numpy array excluding PMLs 
	      Npml[2]:N[1]-Npml[3],
              Npml[4]:N[2]-Npml[5]]
	x, y, z    = np.meshgrid(dual[0], dual[1], dual[2], indexing='ij')     # cell positions excluding PMLs	
	dx, dy, dz = np.meshgrid(prim[0][Npml[0]:N[0]-Npml[1]],
                                 prim[1][Npml[2]:N[1]-Npml[3]],                # cell dimensions excluding PMLs 
                                 prim[2][Npml[4]:N[2]-Npml[5]], indexing='ij') 	             
	for u in np.unique(m)[1:]:	     				
		v = s.mlab.quiver3d(x[m==u]-dx[m==u]/2, y[m==u], z[m==u], dx[m==u], dy[m==u], dz[m==u], 
                                    scale_factor=1, mode='cube', color=(random(),random(),random()))
		v.glyph.glyph.orient = False 
		v.glyph.glyph.scale_mode = 'scale_by_vector_components'        		
		vox.append(v)                   
	return vox



def rectGrid(g):
	''' Create a tvtk rectilinear grid structure. The input is a Grid object. '''
	from tvtk.api import tvtk
	rg = tvtk.RectilinearGrid()	
	rg.dimensions = [a.shape[0] for a in g.edges]
	rg.point_data.scalars = np.zeros(rg.dimensions).ravel()
	rg.point_data.scalars.name = 'scalars'	
	rg.x_coordinates, rg.y_coordinates, rg.z_coordinates = g.edges	
	return rg	



def showGridPlanes(s,g):
	''' Display the edges of the Grid g only in each of the three coordinate planes 
	    of the Cartesian system (s is the Mayavi scene). '''
	gp = [s.mlab.pipeline.grid_plane(rectGrid(g)) for i in range(3)] 
	gp[0].grid_plane.axis = 'x'
	gp[1].grid_plane.axis = 'y'
	gp[2].grid_plane.axis = 'z'
	return gp



###################################################################################################################################################
####### BACKEND CLASSES ###########################################################################################################################
###################################################################################################################################################



class Grid(object):
	
	def __init__(self, bounds=None, maxstep=None):		
		self.bounds  = bounds              # two cartesian points [[x1,y1,z1],[x2,y2,z2]] define a cuboid grid (mm)
		self.maxstep = maxstep             # maximum allowed step inside the grid (user-specified or depending on spw)  
		self.edges   = [0,0,0]             # primary-grid edges (PML layers excluded, numpy arrays, in mm)
		self.steps   = [0,0,0]             # primary-grid steps (PML layers included, lists, in mm)		
		self.matmap  = []                  # material map (PML layers included, integers)
		self.matlut  = []                  # material look-up table (mur,epsr,sigma)
		self.mesh    = []                  # Trimesh objects 		
		self.source  = None                # source  
		self.sensor  = []                  # sensors 
	
	def add_bounds(self):	
		''' Adds the bounds. Only useful for a grid that is already discretized.'''
		self.bounds = [[self.edges[0][0] , self.edges[1][0] , self.edges[2][0] ],
			       [self.edges[0][-1], self.edges[1][-1], self.edges[2][-1]]]	
	
	def contains(self, pt):
		''' Returns True if the point pt [x,y,z] is inside the grid.'''
		return all([self.bounds[0][i]<=pt[i]<=self.bounds[1][i] for i in range(3)])
	         
	def dual_edges(self):
		''' Returns the dual-grid edges, which are exactly in the middle between two 
                    adjacent primary-grid edges in order to yield supraconvergence.'''
		return [(self.edges[i][:-1]+self.edges[i][1:])/2.0 for i in range(3)] 
	
	def add_steps(self):
		''' Adds the primary-grid steps to the Grid object.'''									
		for i in range(3): 
			self.steps[i] = (self.edges[i][1:]-self.edges[i][:-1]).tolist()



class MainGrid(Grid):
	
	def __init__(self, Npml, Kmax, Amax, fmax, bounds=None, maxstep=None):
		super(MainGrid, self).__init__(bounds, maxstep) 
		self.Npml  = Npml                 # number of PML layers along each of the six faces (-x,+x,-y,+y,-z,+z) 
		self.Kmax  = Kmax                 # maximum value of the real stretch along each of the six faces
		self.Amax  = Amax                 # maximum value of the complex-frequency shift along each of the six faces
		self.fmax  = fmax	          # maximum frequency of interest (GHz)	
		self.adhie = Adhie()              # ADHIE implicitization parameters 
		self.dtau  = None
	
	def add_steps(self):
		''' Adds the primary-grid steps (including PML layers) to the Grid object.'''							
		super(MainGrid, self).add_steps()		
		for i in range(3): 
			self.steps[i] = [self.steps[i][0]]*self.Npml[2*i] + self.steps[i] + [self.steps[i][-1]]*self.Npml[2*i+1]					

	def minstep(self):
		''' Returns for each dimension the smallest primary step that is not eliminated by implicitization.'''	
		minstep = []		
		for i in range(3):
			try:    minstep.append(min([s for s in self.steps[i] if s>self.adhie.thr[i]]))	
			except:	minstep.append(None)
		return minstep
	
	def add_dtau_expl(self):
		''' Computes the maximum time step (in mm) without implicitization and adds it to the MainGrid object.'''
		self.dtau = 1.0/sqrt(sum([(cos(pi/2.0/len(self.steps[i]))/min(self.steps[i]))**2 for i in range(3)])) 
	
	def dtau_impl(self):
		''' Computes the maximum time step (in mm) with implicitization.
		    For a fully implicit grid, the Nyquist limit is returned. '''
		minstep = self.minstep()
		try:    return (1.0-self.adhie.alpha)/sqrt(sum([(cos(pi/2.0/len(self.steps[i]))/minstep[i])**2 for i in range(3) if minstep[i]]))
		except: return 1e-6*c0/2/self.fmax    # (mm)           

	def dtau_curlnorm(self):
		''' Computes a more exact time step upper bound (in mm), with implicitization, based on the curl norm.''' 
		pass

	def dtau2dt(self):
		''' Returns the time step in ns.'''
		return self.dtau/c0*1e6



class SubGrid(Grid):
	
	def __init__(self, bounds=None, maxstep=None, impl=None, bi=None):
		super(SubGrid, self).__init__(bounds, maxstep)
		self.impl = impl                  # list of zeros and ones specifying which dimensions are resolved implicitly
		self.bi   = bi                    # main grid's dual indices [i1,j1,k1,i2,j2,k2] of the bounds of the subgrid (PMLs included)
	
	def overlaps(self, bounds):
		''' Returns True if the cuboid delimited by bounds=[[x1,y1,z1],[x2,y2,z2]] overlaps with the subgrid.'''
		if any([bounds[0][i]>bounds[1][i] for i in range(3)]): raise ValueError('Incorrect bounds.')					
		return not (any([(self.bounds[1][i]<bounds[0][i] or self.bounds[0][i]>bounds[1][i]) for i in range(3)]))

	def add_bounds_indices(self, maingrid):
		''' Returns the bounds of the subgrid in terms of the spatial indices [i1,j1,k1,i2,j2,k2] of the 
		    outermost Yee cells inside the main grid (including PMLs) that are covered by the subgrid. '''
		ind = [[np.argmin(np.abs(self.bounds[j][i]-maingrid.edges[i])) for i in range(3)] for j in range(2)] 
		ind = [ind[j][i]+maingrid.Npml[i]-j for i in range(3) for j in range(2)]
		self.ind = ind


class WaveForm(object):                                        
	
	def __init__(self, t, amp=1.0, f0=None, tw=None, t0=None, ts=None):                  				
		self.type = t                     # type		
		self.amp  = amp                   # amplitude
		self.f0   = f0                    # center frequency (GHz)
		self.tw   = tw                    # temporal width (ns)
		self.t0   = t0                    # time delay (ns)
		self.timesignal = ts              # precomputed time signal (list) --> requires dt



class NearField(object):                                        
	
	def __init__(self, t, ind):
		self.type = t                     # source/sensor type
		self.ind  = ind                   # two sets of three grid indices [i1,j1,k1,i2,j2,k2] define a cuboid source/sensor 



class PlaneWave(object):                                        	
	
	def __init__(self, polarization, azimuth, elevation):
		self.pola = polarization          # 0=TE, 1=TM
		self.azim = azimuth               # phi (degrees)         
		self.elev = elevation             # theta (degrees)



class Adhie(object):                                        	
	
	def __init__(self, threshold=[0]*3, alpha=0):      
		self.thr   = threshold            # thresholds for each of the three spatial steps (mm)
		self.alpha = alpha                # ADHIE parameter (not needed for Crank-Nicolson)



class FDTDvisual(object):                                        	
	
	def __init__(self, t, start, stride, scale):      
		self.type   = t                   # field component that should be visualized
		self.start  = start               # (ns)
		self.stride = stride              # (ns)
		self.scale  = scale               # rescaling factor as to ensure a good contrast in the plotted data  



###################################################################################################################################################
####### FRONTEND CLASSES (WIDGETS) ################################################################################################################
###################################################################################################################################################



class Visualization(HasTraits):
	
	scene = Instance(MlabSceneModel, ())                                                                                    

	@on_trait_change('scene.activated')             
	def update_plot(self): 
		self.scene.render()                                                  

	view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
		     width=600, height=600, show_label=False), resizable=True)


class MayaviWindow(QtGui.QWidget):

	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		layout = QtGui.QVBoxLayout(self)
		layout.setContentsMargins(0,0,0,0)    
		layout.setSpacing(0)
		self.visualization = Visualization()
		self.ui = self.visualization.edit_traits(parent=self,kind='subpanel').control
		layout.addWidget(self.ui)
		self.ui.setParent(self)



###################################################################################################################################################



class MaterialInputWidget(QtGui.QWidget):                                        

	def __init__(self):
		super(MaterialInputWidget, self).__init__()                      				
      
		self.layout = QtGui.QGridLayout(self)		
		self.layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
                
		self.layout.addWidget(QtGui.QLabel('material:', self), 0, 0)	

		self.le_mur   = self.createLineEdit(1, u'\u03BC' + '<sub>r</sub> (-): '      , '1')
		self.le_epsr  = self.createLineEdit(2, u'\u03B5' + '<sub>r</sub> (-): '      , '1')
		self.le_sigma = self.createLineEdit(3, u'\u03C3' + ' (S/m): ', '0')        
          
		self.cb_pec = QtGui.QComboBox(self)
		self.cb_pec.addItems(('--','PEC'))
		self.layout.addWidget(self.cb_pec, 4, 1)
		self.cb_pec.currentIndexChanged.connect(self.slot_pec)

		self.layout.addWidget(QtGui.QLabel('', self), 5, 0)
		self.layout.addWidget(QtGui.QLabel('max step (mm):', self), 6, 0)

		self.le_dx = self.createLineEdit(7, '     x: ')
		self.le_dy = self.createLineEdit(8, '     y: ')
		self.le_dz = self.createLineEdit(9, '     z: ')        

		self.layout.addWidget(QtGui.QLabel('', self), 10, 0)

		self.btn = QtGui.QPushButton('assign',self)                              
		self.btn.setFixedSize(80,30)	
		self.layout.addWidget(self.btn,11,0)

		self.btn_all = QtGui.QPushButton('assign to all',self)                              
		self.btn_all.setFixedSize(120,30)	
		self.layout.addWidget(self.btn_all,12,0)	

	def createLineEdit(self, row, text, init=None):
		self.layout.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)  		
		le.setValidator(QtGui.QDoubleValidator(bottom=0))      		
		self.layout.addWidget(le, row, 1)		
		return le

	def slot_pec(self):
		if self.cb_pec.currentIndex() == 0:
			self.le_mur.setEnabled(True)
			self.le_epsr.setEnabled(True)
			self.le_sigma.setEnabled(True)
		else:
			self.le_mur.setEnabled(False);   self.le_mur.setText('1')  
			self.le_epsr.setEnabled(False);  self.le_epsr.setText('inf')
			self.le_sigma.setEnabled(False); self.le_sigma.setText('0')



##################################################################################################################################################



class SubgridInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, mesh, scene):
		super(SubgridInputWidget, self).__init__()                      				

		# main layout is a vertical box
		layout = QtGui.QVBoxLayout(self)		
		layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)

		# divide vertical box into two grids 
		grid1 = QtGui.QWidget()
		grid2 = QtGui.QWidget()	
		self.layout1 = QtGui.QGridLayout(grid1)                  
		self.layout2 = QtGui.QGridLayout(grid2)  
		self.layout1.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
		self.layout2.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
		
		# set up the first grid
		for i,t in enumerate(['x', 'y', 'z']):
			self.layout1.addWidget(QtGui.QLabel(t, self), 0, i+1, alignment=QtCore.Qt.AlignCenter)	
		self.le_pt1 = self.createLineEdits(1, 'point 1 (mm):')
		self.le_pt2 = self.createLineEdits(2, 'point 2 (mm):')
		self.le_max = self.createLineEdits(3, 'max step (mm):', enabled=False)
		self.layout1.addWidget(QtGui.QLabel('implicit:', self), 4, 0)		
		self.chk_impl = [0,0,0]		
		for i in range(3):
			self.chk_impl[i] = QtGui.QCheckBox(self)		
      			self.chk_impl[i].setCheckState(QtCore.Qt.Unchecked)
			self.layout1.addWidget(self.chk_impl[i], 4, i+1, alignment=QtCore.Qt.AlignCenter)	
			self.chk_impl[i].stateChanged.connect(lambda state, j=i: self.slot_enableMaxstep(j))			

		# set up the second grid
		self.btn_select = QtGui.QPushButton('select object', self)                              
		self.btn_select.setFixedSize(160,30) 	
		self.layout2.addWidget(self.btn_select, 5, 0)
		self.connect(self.btn_select, QtCore.SIGNAL('clicked()'), lambda m=mesh: self.slot_selectObject(m))	

		self.layout2.addWidget(QtGui.QLabel('', self), 6, 0)

		self.btn_save = QtGui.QPushButton('save', self)                              
		self.btn_save.setFixedSize(80,30) 	
		self.layout2.addWidget(self.btn_save, 7, 0)
			
		self.btn_quit = QtGui.QPushButton('quit', self)                              
		self.btn_quit.setFixedSize(80,30) 	
		self.layout2.addWidget(self.btn_quit, 8, 0)

		# add everything to the main layout		
		layout.addWidget(grid1)
		layout.addWidget(grid2)		
	
		# visualize the subgrid
		self.box = scene.mlab.outline(line_width=3)
		self.box.bounds = np.zeros(6)                       
		for i in range(3): 
			self.connect(self.le_pt1[i], QtCore.SIGNAL('textChanged(const QString&)'), self.slot_visualize)
			self.connect(self.le_pt2[i], QtCore.SIGNAL('textChanged(const QString&)'), self.slot_visualize)

	def createLineEdit(self, row, column, init, enabled):	
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator())
		le.setEnabled(enabled)    
		le.setMaxLength(20)  		
		self.layout1.addWidget(le, row, column)		
		return le

	def createLineEdits(self, row, text, init=None, enabled=True):
		self.layout1.addWidget(QtGui.QLabel(text, self), row, 0)	
		le_x = self.createLineEdit(row, 1, init, enabled)
		le_y = self.createLineEdit(row, 2, init, enabled)
		le_z = self.createLineEdit(row, 3, init, enabled)		
		return [le_x, le_y, le_z]
	
	def slot_enableMaxstep(self, j):		
		if self.chk_impl[j].isChecked(): self.le_max[j].setEnabled(True)
		else:				 self.le_max[j].setEnabled(False)

	def slot_selectObject(self, mesh):
		m, ok = QtGui.QInputDialog.getInt(self,'Object','Select object number:', 0, 0, len(mesh)-1)
		if ok:
			for i in range(3):			
				self.le_pt1[i].setText(format(mesh[m].bounds[0][i],'.64f'))
				self.le_pt2[i].setText(format(mesh[m].bounds[1][i],'.64f'))
				self.le_pt1[i].setCursorPosition(0)
				self.le_pt2[i].setCursorPosition(0)

	def slot_visualize(self):		
		x1 = le2float(self.le_pt1[0])
		y1 = le2float(self.le_pt1[1])
		z1 = le2float(self.le_pt1[2])
		x2 = le2float(self.le_pt2[0])
		y2 = le2float(self.le_pt2[1])
		z2 = le2float(self.le_pt2[2])
		if x1<x2 and y1<y2 and z1<z2:		
			self.box.bounds = np.asarray([x1,x2,y1,y2,z1,z2])


##################################################################################################################################################



class MainGridInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, mesh, subgrid, scene):
		super(MainGridInputWidget, self).__init__()                      				
      
		# main layout is a vertical box
		layout = QtGui.QVBoxLayout(self)		
		layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)

		# divide vertical box into three grids 
		grid1 = QtGui.QWidget()
		grid2 = QtGui.QWidget()
		grid3 = QtGui.QWidget()		
		self.layout1 = QtGui.QGridLayout(grid1)                  
		self.layout2 = QtGui.QGridLayout(grid2)  
		self.layout3 = QtGui.QGridLayout(grid3)	
		self.layout3.setAlignment(QtCore.Qt.AlignLeft)
		
		# set up first grid     
		self.le_spw  = self.createLineEdit1(0, 'samples per wavelength: ' , '15')				
		self.le_fmax = self.createLineEdit1(1, 'highest frequency (GHz): ', '10')
		self.le_gr   = self.createLineEdit1(2, 'grading ratio: '          , '1.5') 	

		# set up second grid		
		for i,t in enumerate(['d (mm)', 'BC', 'Npml', 'Kmax', 'Amax']):
			self.layout2.addWidget(QtGui.QLabel(t, self), 0, i+1, alignment=QtCore.Qt.AlignCenter)			
		self.le_dist, self.cb, self.le_Npml, self.le_Kmax, self.le_Amax = [], [], [], [], []		
		for i,t in enumerate(['-x', '+x', '-y', '+y', '-z', '+z']):					
			le_dist, cb, le_Npml, le_Kmax, le_Amax = self.createLineEdit2(i+1, t, '0', '0', '5','0.05')  
			self.le_dist.append(le_dist)
			self.cb.append(cb)
			self.le_Npml.append(le_Npml)
			self.le_Kmax.append(le_Kmax)	
			self.le_Amax.append(le_Amax)						

		# set up third grid 	
		self.chk_vox = QtGui.QCheckBox('voxelize', self)
      		self.chk_vox.setCheckState(QtCore.Qt.Unchecked)
		self.layout3.addWidget(self.chk_vox, 0, 0, 1, 3)

		self.btn_mesh = QtGui.QPushButton('start', self)                              
		self.btn_mesh.setFixedSize(80,30) 	
		self.layout3.addWidget(self.btn_mesh, 1, 0)

		self.btn_save = QtGui.QPushButton('save', self)                              
		self.btn_save.setFixedSize(80,30) 
		self.btn_save.setEnabled(False)	
		self.layout3.addWidget(self.btn_save, 1, 1)
		
		self.btn_quit = QtGui.QPushButton('quit', self)                              
		self.btn_quit.setFixedSize(80,30) 
		self.btn_quit.setEnabled(False)	
		self.layout3.addWidget(self.btn_quit, 1, 2)
		
		# add everything to the main layout		
		layout.addWidget(grid1)
		layout.addWidget(QtGui.QLabel('exterior boundaries:', self))
		layout.addWidget(grid2)
		layout.addWidget(grid3)	
	
		# plot axis-aligned bounding box of the full mesh (stb = standard bounds)		
		self.stb = np.asarray([[min([m.bounds[0][i] for m in mesh]+[sg.bounds[0][i] for sg in subgrid]),
					max([m.bounds[1][i] for m in mesh]+[sg.bounds[1][i] for sg in subgrid])] for i in range(3)]).ravel()	
		self.box = scene.mlab.outline(line_width=3, color=(0,0,0))                      		 		
		self.box.bounds = self.stb

		for j in range(6):
			self.connect(self.le_dist[j], QtCore.SIGNAL('textChanged(const QString&)'), lambda state, i=j: self.slot_plotbox(i))

	def createLineEdit1(self, row, text, init):
		self.layout1.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator(bottom=0))      		
		self.layout1.addWidget(le, row, 1)		
		return le

	def createLineEdit2(self, row, text, init_dist, init_Npml, init_Kmax, init_Amax):	
		self.layout2.addWidget(QtGui.QLabel(text, self), row, 0)	
		
		le_dist = QtGui.QLineEdit(init_dist, self)
		le_dist.setValidator(QtGui.QDoubleValidator(bottom=0))      		
		self.layout2.addWidget(le_dist, row, 1)		

		cb = QtGui.QComboBox(self)
		cb.addItems(('PEC','PML'))
		self.layout2.addWidget(cb, row, 2)
	
		le_Npml = QtGui.QLineEdit(init_Npml, self)
		le_Npml.setValidator(QtGui.QIntValidator(bottom=0))      		
		le_Npml.setEnabled(False)		
		self.layout2.addWidget(le_Npml, row, 3)
	
		le_Kmax = QtGui.QLineEdit(init_Kmax, self)
		le_Kmax.setValidator(QtGui.QDoubleValidator(bottom=0))      		
		le_Kmax.setEnabled(False)		
		self.layout2.addWidget(le_Kmax, row, 4)

		le_Amax = QtGui.QLineEdit(init_Amax, self)
		le_Amax.setValidator(QtGui.QDoubleValidator(bottom=0))      		
		le_Amax.setEnabled(False)		
		self.layout2.addWidget(le_Amax, row, 5)

		def slot_EnablePML():
			if cb.currentIndex() ==0: 
				le_Npml.setEnabled(False)
				le_Npml.setText('0')
				le_Kmax.setEnabled(False)
				le_Amax.setEnabled(False)
			else:
				le_Npml.setEnabled(True)
				le_Npml.setText('5')
				le_Kmax.setEnabled(True)
				le_Amax.setEnabled(True)				

		self.connect(cb, QtCore.SIGNAL('currentIndexChanged(const QString&)'), slot_EnablePML)	

		return le_dist, cb, le_Npml, le_Kmax, le_Amax
		
	def slot_plotbox(self, i): 
		b = np.copy(self.box.bounds) 
		b[i] = self.stb[i] + (-1)**(i+1) * le2float(self.le_dist[i]) 
		self.box.bounds = b



##################################################################################################################################################



class WaveFormInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, fmax):
		super(WaveFormInputWidget, self).__init__()                      				
      		
		self.fmax = fmax

		self.layout = QtGui.QGridLayout(self)		
		self.layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
				
		self.layout.addWidget(QtGui.QLabel('time dependence:'), 0, 0)
		self.cb = QtGui.QComboBox(self)
		self.cb.addItems(('monochromatic sine wave',
                                  'gaussian pulse',
                                  'differentiated gaussian pulse',
                                  'gaussian-modulated sinusoidal pulse',
                                  'dirac pulse'))
		self.layout.addWidget(self.cb, 0, 1)

		self.le_amp = self.createLineEdit(1, 'amplitude: ', '1')
		self.le_f0  = self.createLineEdit(2, 'central freq (GHz): ', str(fmax))
		self.le_bw  = self.createLineEdit(3, 'bandwidth (GHz): ' , enabled=False)
		self.le_t0  = self.createLineEdit(4, 'time delay (ns): ' , enabled=False)				    

		self.btn_plot = QtGui.QPushButton('plot', self)                              
		self.btn_plot.setFixedSize(40,30) 
		self.connect(self.btn_plot, QtCore.SIGNAL('clicked()'), self.slot_plot)
		self.layout.addWidget(self.btn_plot, 5, 0)

		self.btn_save = QtGui.QPushButton('save and quit', self)                              
		self.btn_save.setFixedSize(120,30) 
		self.layout.addWidget(self.btn_save, 6, 0)

		self.connect(self.cb, QtCore.SIGNAL('currentIndexChanged(const QString&)'), self.slot_changed)	

	def createLineEdit(self, row, text, init=None, enabled=True):		
		self.layout.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator(bottom=0))  	
		le.setEnabled(enabled)  		
		self.layout.addWidget(le, row, 1)		
		return le
	
	def slot_changed(self):
		if self.cb.currentIndex() == 0:
			self.le_amp.setEnabled(True); self.le_amp.setText('1')			
			self.le_f0.setEnabled(True);  self.le_f0.setText(str(self.fmax))
			self.le_bw.setEnabled(False); self.le_bw.setText(None)
			self.le_t0.setEnabled(False); self.le_t0.setText(None)
		elif self.cb.currentIndex() == 1 or self.cb.currentIndex() == 2:
			self.le_amp.setEnabled(True); self.le_amp.setText('1')			
			self.le_f0.setEnabled(False); self.le_f0.setText(None)
			self.le_bw.setEnabled(True);  self.le_bw.setText(str(self.fmax))
			self.le_t0.setEnabled(True);  self.le_t0.setText(str(4.0/self.fmax))
		elif self.cb.currentIndex() == 3:					
			self.le_amp.setEnabled(True); self.le_amp.setText('1')			
			self.le_f0.setEnabled(True);  self.le_f0.setText(str(self.fmax/7.0))
			self.le_bw.setEnabled(True);  self.le_bw.setText(str(6.0*self.fmax/7.0))
			self.le_t0.setEnabled(True);  self.le_t0.setText(str(7.0/4.0/self.fmax))	
		elif self.cb.currentIndex() == 4:
			self.le_amp.setEnabled(True); self.le_amp.setText('1')
			self.le_f0.setEnabled(False); self.le_f0.setText(None)
			self.le_bw.setEnabled(False); self.le_bw.setText(None)
			self.le_t0.setEnabled(False); self.le_t0.setText(None)	
	
	def slot_plot(self):
		amp = le2float(self.le_amp,1e-100)                    
		f0  = le2float(self.le_f0 ,1e-100)           # units: t is expressed in ns, f in GHz. 
		bw  = le2float(self.le_bw ,1e-100)
		t0  = le2float(self.le_t0 ,1e-100)
		tw  = 2/pi/bw		
		n   = 2**10
		q   = 10                                     # determines the quality of the Fourier transform		
		if self.cb.currentIndex() == 0: 			
			tmax = 5.0/f0			
			t = np.linspace(0, q*tmax, num=n)				
			s = amp * np.piecewise(t, [t<1.0/f0, t>=1.0/f0], 
						  [lambda t: np.sin(2*pi*f0*t)*np.exp(-(4*f0*(t-1/f0))**2), 
						   lambda t: np.sin(2*pi*f0*t)                             ])		
		elif self.cb.currentIndex() == 1:			
			tmax = 2.0*(t0+tw)			
			t = np.linspace(0, q*tmax, num=n)
			s = amp * np.exp(-((t-t0)/tw)**2)
		elif self.cb.currentIndex() == 2:	
			tmax = 2.0*(t0+tw)			
			t = np.linspace(0, q*tmax, num=n)
			s = -amp*sqrt(2)*exp(0.5)*(t-t0)/tw*np.exp(-((t-t0)/tw)**2)	
		elif self.cb.currentIndex() == 3:
			tmax = 2.0*(t0+tw)		
			t = np.linspace(0, q*tmax, num=n)
			s = amp * np.sin(2*pi*f0*t)*np.exp(-((t-t0)/tw)**2)
		elif self.cb.currentIndex() == 4:
			tmax = 1e-9		
			t = np.linspace(0, q*tmax, num=n)
			s = np.zeros(n)
			s[0] = 1  				
		S, f = my_fft(s, t[1]-t[0])		
		my_plot(t, s, f, S, t_max=tmax, f_max=5*self.fmax)



##################################################################################################################################################


class NearFieldInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, nf , maingrid, subgrid, scene):
		super(NearFieldInputWidget, self).__init__()           				

		if nf=='source':
			types = ('Hx','Hy','Hz','Ex','Ey','Ez','Ix','Iy','Iz')     # I = current [A]
			color = (0,0,0.3) 
		elif nf=='sensor':
			types = ('Hx','Hy','Hz','Ex','Ey','Ez','Vx_qs','Vy_qs','Vz_qs','Ix','Iy','Iz','Ix_qs','Iy_qs','Iz_qs')		
			color = (0,0.3,0)			
	
		self.layout = QtGui.QGridLayout(self)		
		self.layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
		
		self.layout.addWidget(QtGui.QLabel('type:'), 0, 0)
		self.cb_type = QtGui.QComboBox(self)
		self.cb_type.addItems(types)        
		self.cb_type.setFixedSize(70,30) 		
		self.layout.addWidget(self.cb_type, 0, 1)

		self.slider = QtGui.QSlider(QtCore.Qt.Horizontal,self)
		self.slider.setMinimum(-30.0)
		self.slider.setMaximum(30.0)
		self.slider.setValue(0.0)
		self.slider.setTickPosition(QtGui.QSlider.TicksBelow)
      		self.slider.setTickInterval(10.0)		
		self.slider.setFixedSize(160,30)		
		self.layout.addWidget(self.slider, 0, 2, 1, 3, alignment=QtCore.Qt.AlignCenter)
		
		self.layout.addWidget(QtGui.QLabel(' ', self), 1, 0)	
		self.layout.addWidget(QtGui.QLabel('desired', self), 2, 1, 1, 2, alignment=QtCore.Qt.AlignCenter)
		self.layout.addWidget(QtGui.QLabel('closest', self), 2, 3, 1, 2, alignment=QtCore.Qt.AlignCenter)

		self.le_position = [0,0,0]
		self.le_position[0] = self.createLineEdits(3, 'x (mm):')
		self.le_position[1] = self.createLineEdits(4, 'y (mm):')
		self.le_position[2] = self.createLineEdits(5, 'z (mm):')		

		self.layout.addWidget(QtGui.QLabel('reference:'), 6, 0)
		self.cb_ref = QtGui.QComboBox(self)
		self.cb_ref.addItems(('stl-file','grid corner'))		
		self.layout.addWidget(self.cb_ref, 6, 1, 1, 2)

		self.layout.addWidget(QtGui.QLabel(' ', self), 7, 0)	

		self.btn_save = QtGui.QPushButton('save and quit', self)                              
		self.btn_save.setFixedSize(120,30) 
		self.layout.addWidget(self.btn_save, 8, 0)	

		self.nf_type = [None]*(len(subgrid)+1)
		self.nf_ind  = [None]*(len(subgrid)+1)
		self.points  = []           		
		
		if nf=='sensor': 
			self.slot_enablePositions()			
			self.connect(self.cb_type, QtCore.SIGNAL('currentIndexChanged(const QString&)'), self.slot_enablePositions)
			self.connect(self.cb_ref , QtCore.SIGNAL('currentIndexChanged(const QString&)'), self.slot_enablePositions)
			for i in range(3): 
				self.connect(self.le_position[i][0], QtCore.SIGNAL('textChanged(const QString&)'), self.slot_enablePositions)
				self.connect(self.le_position[i][1], QtCore.SIGNAL('textChanged(const QString&)'), self.slot_enablePositions)

		self.slot_closest(maingrid, subgrid, color, scene)			
		lambdaslot = lambda state, mg=maingrid, sg=subgrid, co=color, sc=scene: self.slot_closest(mg,sg,co,sc)
		self.connect(self.slider , QtCore.SIGNAL('valueChanged(int)')                  , lambdaslot)			
		self.connect(self.cb_type, QtCore.SIGNAL('currentIndexChanged(const QString&)'), lambdaslot)
		self.connect(self.cb_ref , QtCore.SIGNAL('currentIndexChanged(const QString&)'), lambdaslot)		
		for i in range(3): 
			self.connect(self.le_position[i][0], QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)
			self.connect(self.le_position[i][1], QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)

	def createLineEdit(self, row, column, init, enabled=True):			
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator())  	
		le.setEnabled(enabled) 
		le.setFixedSize(65,30)  	
		self.layout.addWidget(le, row, column)	
		return le

	def createLineEdits(self, row, text, init=['0']*4):		
		self.layout.addWidget(QtGui.QLabel(text, self), row, 0)	
		le_min_desired = self.createLineEdit(row, 1, init[0])
		le_max_desired = self.createLineEdit(row, 2, init[1])
		le_min_closest = self.createLineEdit(row, 3, init[2], enabled=False)
		le_max_closest = self.createLineEdit(row, 4, init[3], enabled=False)
		return [le_min_desired, le_max_desired, le_min_closest, le_max_closest]

	def slot_enablePositions(self):
		t = int(self.cb_type.currentIndex()) 		
		if t<6: 
			for i in range(3): 
				self.le_position[i][1].setEnabled(False)
				self.le_position[i][1].setText(self.le_position[i][0].text()) 
		elif t<9:
			self.le_position[t%3][1].setEnabled(True)
			self.le_position[(t+1)%3][1].setEnabled(False)
			self.le_position[(t+2)%3][1].setEnabled(False)	
			self.le_position[(t+1)%3][1].setText(self.le_position[(t+1)%3][0].text())
			self.le_position[(t+2)%3][1].setText(self.le_position[(t+2)%3][0].text())	
		else:
			self.le_position[t%3][1].setEnabled(False) 
			self.le_position[(t+1)%3][1].setEnabled(True)
			self.le_position[(t+2)%3][1].setEnabled(True)
			self.le_position[t%3][1].setText(self.le_position[t%3][0].text())

	def slot_closest(self, mg, sg, co, sc):
		self.nf_type = [None]*(len(sg)+1)
		self.nf_ind  = [None]*(len(sg)+1)
		for p in self.points: p.stop()
		t = int(self.cb_type.currentIndex()) 						
		r = int(self.cb_ref.currentIndex())			
		d=[[le2float(self.le_position[i][j])+r*mg.edges[i][0] for i in range(3)] for j in range(2)]					
		if any(d[0][i]>d[1][i] for i in range(3)): return		
		for k,g in enumerate([mg]+sg): 
			if k==0 or g.overlaps(d):		
				prim = g.edges                                                 # prim = primary grid edges
				dual = g.dual_edges()                                          # dual = dual grid edges				 
				if t<3:		                                               # e = nf component locations
					if   t%3==0: e = [prim[0][1:-1],dual[1],dual[2]]
					elif t%3==1: e = [dual[0],prim[1][1:-1],dual[2]]
					elif t%3==2: e = [dual[0],dual[1],prim[2][1:-1]]
				else:
					if   t%3==0: e = [dual[0],prim[1][1:-1],prim[2][1:-1]]
					elif t%3==1: e = [prim[0][1:-1],dual[1],prim[2][1:-1]]
					elif t%3==2: e = [prim[0][1:-1],prim[1][1:-1],dual[2]]      							
				ind = [[np.argmin(np.abs(d[j][i]-e[i])) for i in range(3)] for j in range(2)] 				
				for i in range(3): 
					self.le_position[i][2].setText(str(e[i][ind[0][i]]-r*mg.edges[i][0])) 
					self.le_position[i][3].setText(str(e[i][ind[1][i]]-r*mg.edges[i][0]))
					self.le_position[i][2].setCursorPosition(0)
					self.le_position[i][3].setCursorPosition(0)
				x,y,z = [a.ravel() for a in np.meshgrid(e[0][ind[0][0]:ind[1][0]+1],
                                   		    		        e[1][ind[0][1]:ind[1][1]+1],
                                    		                        e[2][ind[0][2]:ind[1][2]+1])]
				if k==0:
					for i in range(len(x)):
						for s in sg:
							if s.contains([x[i],y[i],z[i]]):
								np.delete(x,i)
								np.delete(y,i)
								np.delete(z,i) 			
				self.points.append(sc.mlab.points3d(x,y,z,mode='sphere',color=co,scale_factor=10**(0.1*self.slider.value())))				
				self.nf_type[k] = t				
				self.nf_ind[k]  = [ind[int(i/3)][i%3]+mg.Npml[0::2][i%3] for i in range(6)]



##################################################################################################################################################


 
class PlaneWaveInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, scene):
		super(PlaneWaveInputWidget, self).__init__()                      				

		self.layout = QtGui.QGridLayout(self)		
		self.layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
				
		self.layout.addWidget(QtGui.QLabel('polarization: '), 0, 0)
		self.cb = QtGui.QComboBox(self)
		self.cb.addItems(('TE','TM'))
		self.layout.addWidget(self.cb, 0, 1)

		self.le_azimuth   = self.createLineEdit(1, 'azimuth (phi,degr): ', '0')
		self.le_elevation = self.createLineEdit(2, 'elevation (theta,degr): ', '0')

		self.btn_save = QtGui.QPushButton('save and quit', self)                              
		self.btn_save.setFixedSize(120,30) 
		self.layout.addWidget(self.btn_save, 3, 0)	

		self.arrow = scene.mlab.quiver3d(0,0,0,0,0,1,mode='arrow',color=(0,0,0.3),scale_factor=3)
               	
		lambdaslot = lambda state, s=scene: self.slot_plotarrow(s)
		self.connect(wid_pw.le_azimuth  , QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)
		self.connect(wid_pw.le_elevation, QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)

	def createLineEdit(self, row, text, init=None, enabled=True):		
		self.layout.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator(bottom=0))  	
		le.setEnabled(enabled)  		
		self.layout.addWidget(le, row, 1)		
		return le

	def slot_plotarrow(self, scene): 
		self.arrow.stop()					
		phi   = le2float(wid_pw.le_azimuth)/180.0*pi  				
		theta = le2float(wid_pw.le_elevation)/180.0*pi 
		self.arrow = scene.mlab.quiver3d(0,0,0,sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta),mode='arrow',color=(0,0,0.3),scale_factor=3)



##################################################################################################################################################



class AdhieInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, maingrid):
		super(AdhieInputWidget, self).__init__()                      				

		_layout = QtGui.QVBoxLayout(self)		
		_layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)

		_grid1 = QtGui.QWidget()
		_grid2 = QtGui.QWidget()	
		self.layout1 = QtGui.QGridLayout(_grid1)                  
		self.layout2 = QtGui.QGridLayout(_grid2)  
		
		self.layout1.addWidget(QtGui.QLabel('threshold', self), 0, 1)
		self.layout1.addWidget(QtGui.QLabel('min explicit step', self), 0, 2)

		self.le_thr = [0,0,0]
		self.le_min = [0,0,0]
		self.le_thr[0], self.le_min[0] = self.createLineEdit1(1, 'x (mm):', init='0')
		self.le_thr[1], self.le_min[1] = self.createLineEdit1(2, 'y (mm):', init='0')
		self.le_thr[2], self.le_min[2] = self.createLineEdit1(3, 'z (mm):', init='0') 
		
		self.le_alpha     = self.createLineEdit2(1, 'ADHIE factor:'        , '1'                       ) 
		self.le_dtau      = self.createLineEdit2(2, 'time step (mm):'      , str(maingrid.dtau) ) 	
		self.le_dtaubound = self.createLineEdit2(3, 'time step limit (mm):', enabled=False             ) 		

		self.btn_save = QtGui.QPushButton('save and quit', self)                              
		self.btn_save.setFixedSize(120,30) 
		self.layout2.addWidget(self.btn_save, 4, 0)
			
		self.btn_choose = QtGui.QPushButton('choose', self)                              
		self.btn_choose.setFixedSize(60,30) 
		self.layout2.addWidget(self.btn_choose, 4, 1)
		self.connect(self.btn_choose, QtCore.SIGNAL('clicked()'), self.slot_choose)

		_layout.addWidget(_grid1)
		_layout.addWidget(_grid2)	

		self.slot_changed(maingrid)

		lambdaslot = lambda state, mg=maingrid: self.slot_changed(mg)
		self.connect(self.le_alpha, QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)
		for i in range(3):
			self.connect(self.le_thr[i], QtCore.SIGNAL('textChanged(const QString&)'), lambdaslot)

	def createLineEdit1(self, row, text, init=None):		
		self.layout1.addWidget(QtGui.QLabel(text, self), row, 0)	
		le1 = QtGui.QLineEdit(init, self)
		le2 = QtGui.QLineEdit(self) 
		le1.setValidator(QtGui.QDoubleValidator(bottom=0))  	
		le1.setEnabled(True)  
		le2.setEnabled(False)  				
		self.layout1.addWidget(le1, row, 1)
		self.layout1.addWidget(le2, row, 2)	 			
		return le1, le2

	def createLineEdit2(self, row, text, init=None, enabled=True):			
		self.layout2.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)
		le.setValidator(QtGui.QDoubleValidator(bottom=0))  	
		le.setEnabled(enabled)    
		le.setCursorPosition(0) 				
		self.layout2.addWidget(le, row, 1)	
		return le
	
	def slot_changed(self, mg):
		mg.adhie = Adhie([le2float(self.le_thr[i]) for i in range(3)], le2float(self.le_alpha))	   # the pass-by-object changes mg also outside this function  	
		minstep = mg.minstep()		
		for i in range(3): self.le_min[i].setText(str(minstep[i]))	 		                              
		self.le_dtaubound.setText(str(max(mg.dtau_impl(), mg.dtau)))             
		self.le_dtaubound.setCursorPosition(0) 

	def slot_choose(self):
		self.le_dtau.setText(self.le_dtaubound.text())

		

##################################################################################################################################################



class RunInputWidget(QtGui.QWidget):                                        
	
	def __init__(self, dt, fmax):
		super(RunInputWidget, self).__init__()                      				
	
		self.layout = QtGui.QGridLayout(self)		
		self.layout.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)
				
		self.le_Nt = self.createLineEdit(0, 'number of iterations:','1000', integer=True)		
		self.le_st = self.createLineEdit(1, 'simulated time (ns):', str(le2float(self.le_Nt)*dt))	

		self.layout.addWidget(QtGui.QLabel('visualize field:'), 2, 0)
		self.cb = QtGui.QComboBox(self)
		self.cb.addItems(('--','Ex','Ey','Ez','E'))
		self.layout.addWidget(self.cb, 2, 1)		

		self.le_start   = self.createLineEdit(3, 'visualization start:' , '0' , integer=True)
		self.le_stride  = self.createLineEdit(4, 'visualization stride:', '10', integer=True)
		self.le_scale   = self.createLineEdit(5, 'colormap max (V/m):', '0.001')

		self.btn_start = QtGui.QPushButton('start', self)                              
		self.btn_start.setFixedSize(60,30) 
		self.layout.addWidget(self.btn_start, 6, 0)	
		
		self.btn_plot = QtGui.QPushButton('plot', self)                              
		self.btn_plot.setFixedSize(60,30) 
		self.btn_plot.setEnabled(False)
		self.layout.addWidget(self.btn_plot, 7, 0)

		self.btn_save = QtGui.QPushButton('save', self)                              
		self.btn_save.setFixedSize(60,30) 
		self.btn_plot.setEnabled(False)
		self.layout.addWidget(self.btn_save, 8, 0)

		self.connect(self.le_Nt, QtCore.SIGNAL('editingFinished()'), lambda d=dt: self.slot_change_st(d))
		self.connect(self.le_st, QtCore.SIGNAL('editingFinished()'), lambda d=dt: self.slot_change_Nt(d))	

	def createLineEdit(self, row, text, init=None, enabled=True, integer=False):		
		self.layout.addWidget(QtGui.QLabel(text, self), row, 0)	
		le = QtGui.QLineEdit(init, self)
		if integer: le.setValidator(QtGui.QIntValidator(bottom=0))  	
		else:       le.setValidator(QtGui.QDoubleValidator(bottom=0))   
		le.setEnabled(enabled) 
		le.setMaxLength(19)
		le.setCursorPosition(0) 		
		self.layout.addWidget(le, row, 1)		
		return le

	def slot_change_st(self, dt):
		self.le_st.setText(format(le2float(self.le_Nt)*dt),)
		self.le_st.setCursorPosition(0) 

	def slot_change_Nt(self, dt):
		self.le_Nt.setText(str(int(round(le2float(self.le_st)/dt,0))))




##################################################################################################################################################



class DialogWithList(QtGui.QDialog):

	def __init__(self):
		super(DialogWithList, self).__init__()  
		layout = QtGui.QHBoxLayout(self)		
		layout.addWidget(QtGui.QListWidget(self))
		layout.addWidget(QtGui.QLabel('test',self))
		self.setLayout(layout)




###################################################################################################################################################
####### MAIN WINDOW CLASS #########################################################################################################################
###################################################################################################################################################



class MainWindow(QtGui.QMainWindow):


	def __init__(self):
		super(MainWindow, self).__init__()                     					     
		
		# set up back-end attributes
		self.surfplot = []                         # list containing the plotted Mayavi triangular mesh objects
		self.mesh     = []                         # list containing the Trimesh triangular mesh objects

		# set up central widget with grid layout
		_container = QtGui.QWidget()                          
		_layout = QtGui.QGridLayout(_container)	
		self.setCentralWidget(_container)		
		
		# set up Mayavi widget (layout manager automatically assigns ownership to container)
		win_mayavi = MayaviWindow()		
		_layout.addWidget(win_mayavi,0,1)	 
		self.scene = win_mayavi.visualization.scene 				
		self.scene.mlab.figure(self.scene.mlab.gcf(),bgcolor=(1,1,1),fgcolor=(0,0,0))

		# set up input panel 
		_input = QtGui.QWidget()
		_input.setFixedSize(400,600)	
		_layout.addWidget(_input,0,0)
		self.layout_input = QtGui.QGridLayout(_input)	
		self.layout_input.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)		
		self.layout_input.addWidget(QtGui.QLabel('<font color=gray size=4><b> Input panel </b></font>'),0,0)		
			
		# set up workspace panel
		_workspace = QtGui.QWidget()
		_workspace.setFixedSize(150,600)
		_layout.addWidget(_workspace,0,2)
		self.layout_workspace = QtGui.QVBoxLayout(_workspace)	
		self.layout_workspace.setAlignment(QtCore.Qt.AlignTop|QtCore.Qt.AlignLeft)		
		self.layout_workspace.addWidget(QtGui.QLabel('<font color=gray size=4><b> Workspace </b></font>'),0)	
		self.layout_workspace.addWidget(QtGui.QLabel(' '))		
		self.createWorkspaceEntry('materials', self.slot_ws_materials)
		self.createWorkspaceEntry('subgrids' , self.slot_ws_subgrids )		
		self.createWorkspaceEntry('sources'  , self.slot_ws_sources  )		
		self.createWorkspaceEntry('sensors'  , self.slot_ws_sensors  )

		# create actions
		self.act_openstl      = self.createAction('stl-file'         , self.slot_openstl,  enabled=True)			
		self.act_openhdf5     = self.createAction('hdf5-file'        , self.slot_openhdf5, enabled=True)
		self.act_subgrid      = self.createAction('subgrid'          , self.slot_subgrid)		
		self.act_maingrid     = self.createAction('main grid'        , self.slot_maingrid)
		self.act_waveform     = self.createAction('waveform'         , self.slot_waveform)		
		self.act_NFsource     = self.createAction('near-field source', self.slot_NFsource)
		self.act_planewave    = self.createAction('plane wave'       , self.slot_planewave)		
		self.act_NFsensor     = self.createAction('near-field sensor', self.slot_NFsensor)	
		self.act_FFsensor     = self.createAction('far-field sensor' , self.slot_FFsensor)	
		self.act_energysensor = self.createAction('total energy'     , self.slot_energysensor)		
		self.act_adhie        = self.createAction('adhie'            , self.slot_adhie)
		self.act_parallel     = self.createAction('parallel'         , self.slot_parallel)
		self.act_amr          = self.createAction('amr'              , self.slot_amr)
		self.act_run          = self.createAction('run'              , self.slot_run)

		# add actions to the menu bar 		
		openMenu   = self.menuBar().addMenu('Open');        openMenu.addActions((self.act_openstl, self.act_openhdf5)) 		
		meshMenu   = self.menuBar().addMenu('Grid');        meshMenu.addActions((self.act_subgrid, self.act_maingrid,)) 			
		sourceMenu = self.menuBar().addMenu('Excitation');  sourceMenu.addActions((self.act_waveform, self.act_NFsource, self.act_planewave)) 
		sensorMenu = self.menuBar().addMenu('Sensor');      sensorMenu.addActions((self.act_NFsensor, self.act_FFsensor, self.act_energysensor))
		runMenu    = self.menuBar().addMenu('Run');         runMenu.addActions((self.act_adhie, self.act_parallel, self.act_run)) 	 

		# decorate main window
		self.status = self.statusBar()		
		self.setWindowTitle('SimuWave')		
	

	def createAction(self, text, slot, enabled=False):
		action = QtGui.QAction(text, self)
		self.connect(action, QtCore.SIGNAL('triggered()'), slot)
		action.setEnabled(enabled)	
		return action

	def createWorkspaceEntry(self, text, slot):
		self.btn = QtGui.QPushButton(text, self)                               
		self.btn.setFixedSize(140,30) 
		self.layout_workspace.addWidget(self.btn)
		self.connect(self.btn, QtCore.SIGNAL('clicked()'), slot)

	def message(self, string, time=2000):
		self.status.clearMessage ()	
		self.status.showMessage(string, time)
		print(string)


#------- define workspace slots ------------------------------------------------------------------------------------------------------------------
	
	def slot_ws_materials(self):		
		dia = DialogWithList()		
		dia.exec_()

	def slot_ws_subgrids(self): 
		dia = DialogWithList()		
		dia.exec_()		

	def slot_ws_sources(self): 
		dia = DialogWithList()		
		dia.exec_()

	def slot_ws_sensors(self): 
		dia = DialogWithList()		
		dia.exec_()	


#------- define input panel slots ----------------------------------------------------------------------------------------------------------------

	def slot_openstl(self):	

		# disable actions
		self.act_openstl.setEnabled(False)
		self.act_openhdf5.setEnabled(False)
	
		# load stl-file
		filepath = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', '*.stl')			
		if hasattr(filepath, '__iter__'): filepath = filepath[0]	
		if filepath: 
			fullmesh = trimesh.load_mesh(filepath)		
		else: 
			self.act_openstl.setEnabled(True)
			self.act_openhdf5.setEnabled(True)
			return		

		# split the mesh in submeshes which can all have their own material properties
		if not fullmesh.is_watertight: self.message('The triangular mesh (stl) is not watertight.')
		splitmesh = fullmesh.split(only_watertight=False)	                                         		
			
		# visualize the surface mesh                                 	               			
		for m in splitmesh:		
			self.surfplot += [self.scene.mlab.triangular_mesh( 
					   m.vertices[:,0], m.vertices[:,1], m.vertices[:,2], m.faces, color=(0.5,0.5,0.5), opacity=0.5 )]                          	

		# outline the selected body		
		self.counter = 0		
		outline = self.scene.mlab.outline(line_width=3)                       
		outline.outline_mode = 'cornered'	
		outline.bounds = splitmesh[self.counter].bounds.T.ravel() 	

	  	# material assignment  		     			
	        def slot_singleMaterialAssignment():	
			mur   = le2float(wid_mat.le_mur  ,-1)			
			epsr  = le2float(wid_mat.le_epsr ,-1)
			sigma = le2float(wid_mat.le_sigma,-1)
			pec   = bool(wid_mat.cb_pec.currentIndex())
			dx    = max(le2float(wid_mat.le_dx,1e100),tol)
			dy    = max(le2float(wid_mat.le_dy,1e100),tol)
			dz    = max(le2float(wid_mat.le_dz,1e100),tol)
			co    = (random(),random(),random())
			op    = 0.5 if not (mur==1 and epsr==1 and sigma==0) else 0 
			if mur>=1 and epsr>=1 and sigma>=0:		
				m = splitmesh[self.counter]				
				m.mur     = mur				
				m.epsr    = epsr
				m.sigma   = sigma
				m.pec     = pec				
				m.maxstep = [dx,dy,dz]
				self.surfplot[self.counter].stop()
				self.surfplot[self.counter] = self.scene.mlab.triangular_mesh( 
                                                               m.vertices[:,0], m.vertices[:,1], m.vertices[:,2], m.faces, color=co, opacity=op)	
				self.counter+=1		
				if self.counter==splitmesh.shape[0]: 			
					self.mesh = np.append(self.mesh, splitmesh) 
					self.message('Mesh object(s) added.')
					outline.stop()
					self.act_openstl.setEnabled(True)
					self.act_subgrid.setEnabled(True)
					self.act_maingrid.setEnabled(True)
					deleteAttr(self,['counter'])			
					wid_mat.deleteLater()				
				else:
					outline.bounds = splitmesh[self.counter].bounds.T.ravel() 										
			else: 
				QtGui.QMessageBox.critical(self, 'Warning', 'Insert physical values for mu_r, eps_r and sigma.')
        
		def slot_fullMaterialAssignment():	
			mur   = le2float(wid_mat.le_mur  ,-1)			
			epsr  = le2float(wid_mat.le_epsr ,-1)
			sigma = le2float(wid_mat.le_sigma,-1)
			pec   = bool(wid_mat.cb_pec.currentIndex())
			dx    = max(le2float(wid_mat.le_dx,1e100),tol)
			dy    = max(le2float(wid_mat.le_dy,1e100),tol)
			dz    = max(le2float(wid_mat.le_dz,1e100),tol)
			co    = (random(),random(),random())
			op    = 0.5 if not (mur==1 and epsr==1 and sigma==0) else 0
			if mur>=1 and epsr>=1 and sigma>=0:		
				for i,m in enumerate(splitmesh):				
					m.mur     = mur				
					m.epsr    = epsr
					m.sigma   = sigma
					m.pec     = pec		
					m.maxstep = [dx,dy,dz]
					self.surfplot[i].stop()
					self.surfplot[i] = self.scene.mlab.triangular_mesh( 
							    m.vertices[:,0], m.vertices[:,1], m.vertices[:,2], m.faces, color=co, opacity=op)    
				self.mesh = np.append(self.mesh, splitmesh)
				self.message('Mesh object(s) added.')  
				outline.stop()                                        						
				self.act_openstl.setEnabled(True)
				self.act_subgrid.setEnabled(True)	
				self.act_maingrid.setEnabled(True)	
				wid_mat.deleteLater()											
			else: 
				QtGui.QMessageBox.critical(self, 'Warning', 'Insert physical values for mu_r, eps_r and sigma.')

		wid_mat = MaterialInputWidget()			
		self.layout_input.addWidget(wid_mat)		 					
		self.connect(wid_mat.btn    , QtCore.SIGNAL('clicked()'), slot_singleMaterialAssignment)	
		self.connect(wid_mat.btn_all, QtCore.SIGNAL('clicked()'), slot_fullMaterialAssignment  )

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_openhdf5(self):		

		# disable actions
		self.act_openstl.setEnabled(False)
		self.act_openhdf5.setEnabled(False)

		# read hdf5-file
		filepath = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', '*.hdf5')			
		if hasattr(filepath, '__iter__'): filepath = filepath[0]	
		if filepath: 		
			with h5py.File(filepath, 'r') as h5f:
				self.maingrid = MainGrid( list(h5f['Npml'][:]), 
						          list(h5f['Kmax'][:]), 
						          list(h5f['Amax'][:]), 
							       h5f['fmax'][0]   )				
				Ng=0 
				while 'edge_x_'+str(Ng) in h5f: Ng=Ng+1 				
				self.subgrid = [SubGrid()]*(Ng-1)	
				self.message(str(Ng-1) + ' subgrids.')			
				for i,g in enumerate([self.maingrid]+self.subgrid):					
					g.edges[0] = h5f['edge_x_'+str(i)][:]						
					g.edges[1] = h5f['edge_y_'+str(i)][:]		
					g.edges[2] = h5f['edge_z_'+str(i)][:]	
					g.matmap = h5f['matmap_'+str(i)][:]	 		
					g.matlut = [tuple(l) for l in h5f['matlut_'+str(i)][:]]
				        g.add_bounds()
					g.add_steps()
				self.maingrid.add_dtau_expl()
				for s in self.subgrid: s.add_bounds_indices(self.maingrid) 
		else: 
			self.act_openstl.setEnabled(True)
			self.act_openhdf5.setEnabled(True)
			return	

		# gridplanes
		self.gridplanes=[]		
		for g in [self.maingrid]+self.subgrid:
			self.gridplanes.extend(showGridPlanes(self.scene, g))							

		# voxelization
		self.voxels=[]								   																					       			
		for g in [self.maingrid]+self.subgrid: 
			self.voxels.extend(showVoxels(self.scene, g))

		# enable actions
		self.act_waveform.setEnabled(True)	
		self.act_NFsource.setEnabled(True)	

#--------------------------------------------------------------------------------------------------------------------------------------------------
	
	def slot_subgrid(self): 	
		
		# disable actions
		self.act_openstl.setEnabled(False)
		self.act_openhdf5.setEnabled(False)
		self.act_subgrid.setEnabled(False)
		self.act_maingrid.setEnabled(False)

		# subgrid input
		wid_sg = SubgridInputWidget(self.mesh, self.scene)
		self.layout_input.addWidget(wid_sg)
		
		# save
		def slot_save():
			_s = SubGrid( bounds = [[le2float(wid_sg.le_pt1[i]) for i in range(3)],
						[le2float(wid_sg.le_pt2[i]) for i in range(3)]], 				
				      maxstep = [max(le2float(wid_sg.le_max[i],1e100),tol) for i in range(3)], 
                                      impl = [int(wid_sg.chk_impl[i].isChecked()) for i in range(3)] )
			if  any([_s.bounds[0][i]>=_s.bounds[1][i] for i in range(3)]):
				QtGui.QMessageBox.critical(self,'Bad bounds',
				'Every coordinate of point 2 should be higher than that of point 1.') 
				return
			if not any([chk.isChecked() for chk in wid_sg.chk_impl]): 
				QtGui.QMessageBox.critical(self,'No implicitization',
				'At least one dimension should be implicit.') 
				return
			wid_sg.btn_save.setEnabled(False)
			partial_mesh_overlap = False
			for m in self.mesh:
				partial_mesh_overlap = _s.overlaps(m.bounds) and not (_s.contains(m.bounds[0]) and _s.contains(m.bounds[1])) 							       	 
				if partial_mesh_overlap:	
					ans = QtGui.QMessageBox.question(self,'Partial mesh overlap',
					'Subgrids are not yet allowed to partially overlap with ' +
		                        'the axis-aligned bounding box of a mesh object. ' +
					'Only the case for full overlap is correctly implemented. ' +
					'Are you sure you want to add the subgrid anyway?' , 
					 QtGui.QMessageBox.Yes| QtGui.QMessageBox.No)
					break
			if not partial_mesh_overlap or ans==QtGui.QMessageBox.Yes:
				if hasattr(self,'subgrid'): 				
					subgrid_overlap = False				
					for i,s in enumerate(self.subgrid):
						subgrid_overlap = _s.overlaps(s.bounds)
						if subgrid_overlap:	
							QtGui.QMessageBox.critical(self,'Subgrid overlap',
							'Your subgrid overlaps with subgrid ' + str(i) + '.' +
					                'Redefine your current subgrid or delete this already ' +
							'existing subgrid from the workspace.')
							break
					if not subgrid_overlap:				
						self.subgrid.append(_s)
						self.message('Subgrid added.') 					                                    
				else:		         				
					self.subgrid = [_s]
					self.message('Subgrid added.')	
			wid_sg.btn_save.setEnabled(True)			
		self.connect(wid_sg.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

		# enable actions and quit
		def slot_quit():
			wid_sg.box.stop()		
			self.act_subgrid.setEnabled(True)
			self.act_maingrid.setEnabled(True)				
			wid_sg.deleteLater()			
		self.connect(wid_sg.btn_quit, QtCore.SIGNAL('clicked()'), slot_quit)		

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_maingrid(self):

		# disable actions
		self.act_openstl.setEnabled(False)
		self.act_openhdf5.setEnabled(False)
		self.act_subgrid.setEnabled(False)
		self.act_maingrid.setEnabled(False)
		
		# if non-existing, declare a subgrid list
		if not hasattr(self, 'subgrid'): self.subgrid = []

		# print useful information
		self.message(str(len(self.mesh)) + ' mesh objects, ' + str(len(self.subgrid)) + ' subgrids')

		# ugly work-around to activate an object in the scene (needed for outline())
		self.scene.mlab.points3d(0,0,0,opacity=0)

		# main grid input	     			
		wid_mg = MainGridInputWidget(self.mesh, self.subgrid, self.scene)	
		self.layout_input.addWidget(wid_mg)	
		
		# perform the actual grid generation and visualization
		def slot_startGridGeneration():									
			
			# disable buttons while meshing
			wid_mg.btn_mesh.setEnabled(False)
			wid_mg.btn_save.setEnabled(False)
			wid_mg.btn_quit.setEnabled(False)	

			# import data from the line edit forms
			spw  = le2float(wid_mg.le_spw,-1)  
			fmax = le2float(wid_mg.le_fmax,-1) 
			gr   = le2float(wid_mg.le_gr,-1) 			  				
			dist = [     le2float(e)  for e in wid_mg.le_dist ]			
			Npml = [ int(le2float(e)) for e in wid_mg.le_Npml ]
			Kmax = [     le2float(e)  for e in wid_mg.le_Kmax ]			
			Amax = [     le2float(e)  for e in wid_mg.le_Amax ]					     			
			if spw<10: 
				QtGui.QMessageBox.critical(self,'Warning',
					'Numerical dispersion is expected to pollute your computations. ' +
					'The generally accepted lower bound is 10 samples per wavelength, ' +
					'but 20 or more are preferred.')
			elif gr<1.2: 
				QtGui.QMessageBox.critical(self,'Warning',
					'The grading ratio must be larger or equal to 1.2 in order to avoid ' + 
					'excessive mesh generation times.')
			elif any(e>20 for e in Npml):
				QtGui.QMessageBox.critical(self,'Warning',
					'Generally 5-10 PML layers are used. Beyond 20 layers, the reflection ' +
					'error is hardly affected and the computational burden is huge.')
			elif any((e<1 or e>20) for e in Kmax):
				QtGui.QMessageBox.critical(self,'Warning',
					'Kappa is the real part of the stretching factor and serves to absorb ' +
					'evanescent waves. It is polynomially graded inside the PML up to the ' +
					'specified value Kmax. For example, for long waveguides exclusively supporting ' +
                                        'propagating modes, the best absorption is achieved for Kmax=1. However, ' + 
                                        'for short waveguides with insufficiently damped evanescent modes, for ' +  
                                        'good conductors inside the PML, for small dipoles close to the PML, and ' +  
					'so on ... kappa=5-15 is a good choice.')
			elif any((e<0 or e>1) for e in Amax):
				QtGui.QMessageBox.critical(self,'Warning',
					'The CFS-PML parameter "a" downgrades the theoretical performance of the PML at ' +
                                        'low frequencies, but enables the practical absorption of strongly evanescent waves ' + 
                                        'on a grid consisting of finite-sized cells. A good choice is Amax = 0.05-0.30.')
			else:				
				self.maingrid = MainGrid(Npml, Kmax, Amax, fmax, maxstep=[c0/spw/fmax*1e-6]*3)          				

				# actual grid generation 				
				gridGeneration(self.maingrid, self.subgrid, self.mesh, gr, dist)			
				
				Nx = self.maingrid.edges[0].size - 1	
				Ny = self.maingrid.edges[1].size - 1
				Nz = self.maingrid.edges[2].size - 1

				if (Nx-1)*(Ny-1)*(Nz-1)>0: 			
					
					# add step lists, bounds and time step 
					for g in [self.maingrid]+self.subgrid: g.add_steps()
					self.maingrid.add_bounds()
					self.maingrid.add_dtau_expl()
					for s in self.subgrid: s.add_bounds_indices(self.maingrid)	

					# gridplanes
					if hasattr(self,'gridplanes'): 
						for gp in self.gridplanes: gp.stop()
					else:
						self.gridplanes=[]				
					for g in [self.maingrid]+self.subgrid:
						self.gridplanes.extend(showGridPlanes(self.scene, g))					
					
					# voxelization
					if hasattr(self,'voxels'):
						for v in self.voxels: v.stop() 	
					else:
						self.voxels=[]					
					if wid_mg.chk_vox.isChecked():
						for sp in self.surfplot: sp.stop()			   																					       			
						for g in [self.maingrid]+self.subgrid: self.voxels.extend(showVoxels(self.scene, g))		

					# let the user know that the mesh is finished and show some useful information 					
					self.message('Meshing is finished. Number of cells (PML excluded): ' + str((Nx,Ny,Nz)), time=10000)					

					# activate save and quit buttons			
					wid_mg.btn_save.setEnabled(True)
					wid_mg.btn_quit.setEnabled(True)

				else:
					QtGui.QMessageBox.critical(self,'Warning','Invalid discretization. Increase the sampling density.');  

			# enable mesh button
			wid_mg.btn_mesh.setEnabled(True)					

		self.connect(wid_mg.btn_mesh, QtCore.SIGNAL('clicked()'), slot_startGridGeneration)
		
		# save the FDTD mesh to a hdf5-file                                                        
		def slot_save2hdf5():
			filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save file', '', '*.hdf5')			
			if hasattr(filepath, '__iter__'): filepath = filepath[0]	
			if filepath: 
				with  h5py.File(filepath, 'w') as h5f:
					for i,g in enumerate([self.maingrid]+self.subgrid):					
						h5f.create_dataset('edge_x_'+str(i), data=g.edges[0], dtype='float64' )
						h5f.create_dataset('edge_y_'+str(i), data=g.edges[1], dtype='float64' )
						h5f.create_dataset('edge_z_'+str(i), data=g.edges[2], dtype='float64' )
						h5f.create_dataset('matmap_'+str(i), data=g.matmap     , dtype='uint32'  ) 	
						h5f.create_dataset('matlut_'+str(i), data=np.asarray(g.matlut), dtype='float64' ) 					
					h5f.create_dataset('Npml', data=self.maingrid.Npml  , dtype='uint32'  )
					h5f.create_dataset('Kmax', data=self.maingrid.Kmax  , dtype='float64' ) 
					h5f.create_dataset('Amax', data=self.maingrid.Amax  , dtype='float64' )
					h5f.create_dataset('fmax', data=[self.maingrid.fmax], dtype='float64' ) 
		self.connect(wid_mg.btn_save, QtCore.SIGNAL('clicked()'), slot_save2hdf5)

		# quit mesher 
		def slot_quit():				
			# enable actions			
			self.act_waveform.setEnabled(True)		
			self.act_NFsource.setEnabled(True) 			
			# delete attributes that are not needed anymore
			deleteAttr(self,['gridplanes','voxels','mesh'])						
			wid_mg.deleteLater()
		
		self.connect(wid_mg.btn_quit, QtCore.SIGNAL('clicked()'), slot_quit)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_waveform(self):
		wid_wf = WaveFormInputWidget(self.maingrid.fmax)
		self.layout_input.addWidget(wid_wf)
		
		def slot_save():		
			self.wf = WaveForm( int(wid_wf.cb.currentIndex()), 
				            le2float(wid_wf.le_amp), 
				            le2float(wid_wf.le_f0), 
				            2/pi/max(le2float(wid_wf.le_bw),1e-100),
					    le2float(wid_wf.le_t0) ) 		
			wid_wf.deleteLater()	
		self.connect(wid_wf.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_planewave(self):
		wid_pw = PlaneWaveInputWidget(self.scene)
		self.layout_input.addWidget(wid_pw)

		def slot_save():		
			self.pw = PlaneWave( int(wid_pw.cb.currentIndex()),
			                     le2float(wid_pw.le_azimuth),
                                             le2float(wid_pw.le_elevation) )	
			self.arrow.stop()					
			# enable actions			
			self.act_NFsensor.setEnabled(True)
			self.act_adhie.setEnabled(True)	
			self.act_run.setEnabled(True)				
			wid_pw.deleteLater()	
		self.connect(wid_pw.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_NFsource(self):		
		wid_src = NearFieldInputWidget('source', self.maingrid, self.subgrid, self.scene)   
		self.layout_input.addWidget(wid_src)			

		def slot_save():	  			
			added = False		
			for k,g in enumerate([self.maingrid]+self.subgrid):			
				if wid_src.nf_type[k] is not None:
					g.source = NearField(wid_src.nf_type[k], wid_src.nf_ind[k])
					added = True
			if added: 
				self.message('Near-field source modified.')							 							
				for p in wid_src.points: p.stop()			
				self.act_NFsensor.setEnabled(True)
				self.act_adhie.setEnabled(True)	
				self.act_run.setEnabled(True)				
				wid_src.deleteLater()	

		self.connect(wid_src.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_NFsensor(self):		
		wid_nfs = NearFieldInputWidget('sensor', self.maingrid, self.subgrid, self.scene)   
		self.layout_input.addWidget(wid_nfs)					
	
		def slot_save():
			added = False	 			
		        for k,g in enumerate([self.maingrid]+self.subgrid):
				if wid_nfs.nf_type[k] is not None:			
					g.sensor.append(NearField(wid_nfs.nf_type[k], wid_nfs.nf_ind[k])) 
					added = True
			if added: 
				self.message('Near-field sensor added.')			 							
				for p in wid_nfs.points: p.stop() 							
				self.act_adhie.setEnabled(True)	
				self.act_run.setEnabled(True)				
				wid_nfs.deleteLater()
	
		self.connect(wid_nfs.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_FFsensor(self): pass    # NTFF-transformation

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_energysensor(self): pass

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_adhie(self):
		self.maingrid.add_dtau_expl()
		wid_adhie = AdhieInputWidget(self.maingrid)
		self.layout_input.addWidget(wid_adhie)

		def slot_save():
			alpha = le2float(wid_adhie.le_alpha)
			if not 0<alpha<=1:									
				QtGui.QMessageBox.critical(self,'Warning','Choose alpha in the interval ]0,1].')
			else:
				self.maingrid.adhie = Adhie([le2float(wid_adhie.le_thr[i]) for i in range(3)], alpha)			
				self.maingrid.dtau = le2float(wid_adhie.le_dtau)     # (mm)								
				wid_adhie.deleteLater()	
		self.connect(wid_adhie.btn_save, QtCore.SIGNAL('clicked()'), slot_save)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_parallel(self):					
		# OpenCL
		pass

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_amr(self):					
		# adaptive mesh refinement (finite-integration subgridding with refinement ratio 2 and local time stepping)
		pass

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def slot_run(self):
		wid_run = RunInputWidget(self.maingrid.dtau2dt(), self.maingrid.fmax)
		self.layout_input.addWidget(wid_run)

		def slot_FDTDsimulation():
			
			self.Nt  = int(le2float(wid_run.le_Nt))        	
			self.vis = FDTDvisual( int(wid_run.cb.currentIndex()),			
			                       int(le2float(wid_run.le_start,0)), 
		                               int(le2float(wid_run.le_stride,1)), 			
			                       le2float(wid_run.le_scale,1) ) 			

			# add a standard waveform if the user did not specify one
			if not hasattr(self,'wf'):
				self.wf = WaveForm(2, tw=2.0/pi/self.maingrid.fmax, t0=4.0/self.maingrid.fmax)

			if self.Nt>0 and self.vis.start>=0 and self.vis.stride>0:			  			
									
				# disable buttons
				wid_run.btn_start.setEnabled(False)
				wid_run.btn_plot.setEnabled(False)
				wid_run.btn_save.setEnabled(False)
										               
				# precompute source (avoid if-statements in the timestepping loop)
				t = np.arange(self.Nt)*self.maingrid.dtau2dt()                             															
				if self.wf.type == 0: 									
					s = np.piecewise(t, [t<1.0/self.wf.f0, t>=1.0/self.wf.f0], 
                                                            [lambda t: np.sin(2*pi*self.wf.f0*t)*np.exp(-(4*self.wf.f0*(t-1.0/self.wf.f0))**2), 
                                                             lambda t: np.sin(2*pi*self.wf.f0*t)])		
				elif self.wf.type == 1:			
					s =  np.exp(-((t-self.wf.t0)/self.wf.tw)**2)
				elif self.wf.type == 2:	
					s = -sqrt(2)*exp(0.5)*(t-self.wf.t0)/self.wf.tw*np.exp(-((t-self.wf.t0)/self.wf.tw)**2)	
				elif self.wf.type == 3:
					s =  np.sin(2*pi*self.wf.f0*t)*np.exp(-((t-self.wf.t0)/self.wf.tw)**2)
				elif self.wf.type == 4:
					s = np.zeros(self.Nt) 
                                        s[0] = 1.0 
				s *= self.wf.amp  
				self.wf.timesignal = s.tolist()							
	
				# actual FDTD simulation
				print('FDTD simulation is started.') 				
				time_start = time()				
				self.output = run_fdtd_python(self)      # FDTD C-module				
				time_end = time()  
				self.CPUtime = time_end-time_start
				self.message('FDTD simulation is finished. Elapsed time: ' + str(self.CPUtime) + 's.', time=10000)

				# enable buttons
				wid_run.btn_start.setEnabled(True)
				wid_run.btn_plot.setEnabled(True)
				wid_run.btn_save.setEnabled(True)		
				
		self.connect(wid_run.btn_start, QtCore.SIGNAL('clicked()'), slot_FDTDsimulation)
 
		def slot_plotOutput():
			if len(self.maingrid.sensor)>0: 
				s, ok = QtGui.QInputDialog.getInt(self,'Sensor','Select sensor number:', 0, 0, len(self.maingrid.sensor)-1)
				if ok:	
					dt = self.maingrid.dtau2dt()		
					t = np.arange(self.Nt)*dt
					out = np.asarray(self.output).reshape((self.Nt,len(self.maingrid.sensor)))   				
					t_sensor = out[:,s]								
					f_sensor, f = my_fft(t_sensor,dt)				
					t_scale = closestPow10(np.absolute(t_sensor).max()) 
					f_scale = closestPow10(f_sensor.max())              			
					if   self.maingrid.sensor[s].type<3: units = 'A/m' 
					elif self.maingrid.sensor[s].type<6: units = 'V/m'
					elif self.maingrid.sensor[s].type<9: units = 'V'
					else:     		             units = 'A'    
					my_plot(t, t_sensor/t_scale, f, f_sensor/f_scale, f_max=self.maingrid.fmax, 
						  t_label='amplitude ['+str(t_scale)+units+']', 
						  f_label='amplitude ['+str(f_scale)+units+']')
			else:
				QtGui.QMessageBox.critical(self,'Information','There are no sensors.')

		self.connect(wid_run.btn_plot, QtCore.SIGNAL('clicked()'), slot_plotOutput)

		def slot_save2hdf5():
			if len(self.maingrid.sensor)>0: 
				out = np.asarray(self.output).reshape((self.Nt,len(self.maingrid.sensor)))
				filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save file', '', '*.hdf5')			
				if hasattr(filepath, '__iter__'): filepath = filepath[0]	
				if filepath: 			
					with  h5py.File(filepath, 'w') as h5f:
						h5f.create_dataset('data'    , data=out                       , dtype='float64')
						h5f.create_dataset('dt'      , data=[self.maingrid.dtau2dt()] , dtype='float64')  # (ns)
						h5f.create_dataset('CPU time', data=[self.CPUtime]            , dtype='float64')
			else:
				QtGui.QMessageBox.critical(self,'Information','There are no sensors.')

		self.connect(wid_run.btn_save, QtCore.SIGNAL('clicked()'), slot_save2hdf5)

#--------------------------------------------------------------------------------------------------------------------------------------------------

	def closeEvent(self, event):
	    	r = QtGui.QMessageBox.question(self,
		              'Exit',
		              'Are you sure you want to exit ?',
		              QtGui.QMessageBox.Yes| QtGui.QMessageBox.No)
		if r==QtGui.QMessageBox.Yes: event.accept()
		else:                        event.ignore() 






