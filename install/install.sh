#!/bin/bash

cd packages

sudo dpkg -i cython_0.23.4-0ubuntu5_amd64.deb

sudo dpkg -i python-pip_8.1.1-2ubuntu0.4_all.deb
sudo dpkg -i python-numpy_1%3a1.11.0-1ubuntu1_amd64.deb
sudo dpkg -i python-matplotlib_1.5.1-1ubuntu1_amd64.deb 
sudo dpkg -i python-h5py_2.6.0-1_amd64.deb

sudo dpkg -i python-vtk_5.10.1+dfsg-2.1build1_amd64.deb
sudo dpkg -i python-setuptools_20.7.0-1_all.deb
sudo dpkg -i python-traits_4.5.0-1ubuntu2_amd64.deb  
sudo dpkg -i python-qt4_4.11.4+dfsg-1build4_amd64.deb  
sudo dpkg -i python-qt4-gl_4.11.4+dfsg-1build4_amd64.deb  
sudo dpkg -i python-pyface_4.5.2-1_all.deb  
sudo dpkg -i python-configobj_5.0.6-2_all.deb
sudo -H pip2 install mayavi-4.5.0.tar.gz
  
sudo dpkg -i libspatialindex-dev_1.8.5-3_amd64.deb 
sudo dpkg -i libgeos-dev_3.5.0-1ubuntu2_amd64.deb 
sudo dpkg -i python-rtree_0.8.2+ds-2_all.deb 
sudo dpkg -i python-scipy_0.17.0-1_amd64.deb  
sudo dpkg -i python-networkx_1.11-1ubuntu1_all.deb  
sudo -H pip2 install trimesh-2.10.18.tar.gz

