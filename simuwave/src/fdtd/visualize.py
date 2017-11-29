
#Issues: 
#  How to remove the incorrect axes ticks? 
#  How to preserve the ctf and otf efficiently?



from numpy import zeros, reshape, asarray, diff, amax, arange             
from scipy.interpolate import interp1d
from mayavi import mlab                    
from tvtk.util.ctf import ColorTransferFunction, PiecewiseFunction	
from pyface.qt import QtGui

from os import environ
environ['ETS_TOOLKIT'] = 'qt4'




def NewColorMap(h):
	
	ctf = ColorTransferFunction()                                      
	otf = PiecewiseFunction()      	
	ctf.remove_all_points()	
	otf.remove_all_points()		
	ctf.add_rgb_point(0,1,1,1) 	                   
	otf.add_point(0,0) 	
	for i in range(128):			
		q=(i+1.0)/128.0		
		ctf.add_rgb_point(-q, 1-q*q, 1-q*q, 1    )               # -1 --> blue,  0 --> white,  1-->red            	                    	
		ctf.add_rgb_point( q, 1    , 1-q*q, 1-q*q)    
		otf.add_point(-q, 0.5*q)                                   	
		otf.add_point( q, 0.5*q)                  				
	ctf.range = [-1,1]	
	h._volume_property.set_color(ctf)		
	h._ctf = ctf		
	h.update_ctf = True
	h._otf = otf
        h._volume_property.set_scalar_opacity(otf)






def initVisualization(mg):		

	global h, nc, xc, yc, zc, xv, yv, zv 


	### INTERPOLATION ###

	# cell centers 	
	[xc,yc,zc] = mg.dual_edges()

	# number of cells in each dimension
	nc = (xc.size, yc.size, zc.size)

	# voxel size
	d = min(amax(diff(xc)),amax(diff(yc)),amax(diff(zc)))

	# voxel centers
	xv = arange(start=xc[0], stop=xc[-1], step=d)
	yv = arange(start=yc[0], stop=yc[-1], step=d)
	zv = arange(start=zc[0], stop=zc[-1], step=d)

	# number of voxels in each dimension
	nv = (xv.size, yv.size, zv.size)


	### VISUALIZATION ###
	
	# zero-initialize the 3D scalar data visualization
	fig = mlab.figure(1,bgcolor=(1,1,1),fgcolor=(0,0,0),size=(800,800))						
	mlab.clf(fig)			
	h = mlab.pipeline.volume(mlab.pipeline.scalar_field(zeros(nv))) 		

	# change the color and opacity transfer functions 
	NewColorMap(h)			                        
	
	# preserve the resolution of the input data 	
	h.volume_mapper.lock_sample_distance_to_input_spacing = True

	# interpolation (nearest or linear) 
	h.volume_property.interpolation_type = 'linear'
	
	# disable shade
	h.volume_property.shade = False

	# decorations			 	  
	mlab.axes(xlabel='x',ylabel='y',zlabel='z')





def visualize(l):

	# construct 3d array containing the Yee-cell data
	a = reshape(asarray(l),nc)

	# interpolate Yee-cell data to voxel data
	f = interp1d(xc, a, axis=0, kind='linear');  a = f(xv)       # linear, nearest, zero, slinear, quadratic, cubic 
	f = interp1d(yc, a, axis=1, kind='linear');  a = f(yv)
	f = interp1d(zc, a, axis=2, kind='linear');  a = f(zv)

	# actual volume rendering 	
	h.mlab_source.scalars =	a               

	# change the color and opacity transfer functions  			
	NewColorMap(h)		                                       			

	# real-time on-screen rendering
	QtGui.QApplication.processEvents()	
  



'''
Important python files:
/usr/local/lib/python2.7/dist-packages/mayavi/tools/modules.py
/usr/local/lib/python2.7/dist-packages/mayavi/tools/show.py
'''


# speed up tricks when saving images: 
#   fig.scene.disable_render=True   
#   fig.scene.off_screen_rendering=True  or  mlab.options.offscreen=True
#   fig.scene.anti_aliasing_frames=0


