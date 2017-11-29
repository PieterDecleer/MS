
Requires Python 2.7



==== STANDARD INSTALLATION =================================================================================


Installation of SimuWave with Ubuntu16.04 terminal (amd64):
  chmod +x install.sh
  sudo bash install.sh 

Installation on other operating systems and other hardware (via pip): 
  for example:
    pip install trimesh==2.10.18



===== MAKE CHANGES TO SIMUWAVE =============================================================================


Info:
  Mayavi is embedded in C via the Python/C API.
  C is embedded in Python with a Cython wrapper. 

Make changes to the C-code:
  add the c-files to the ext-list in 'cython_setup.py'
  and enter 'python cython_setup.py build_ext --inplace' in the terminal



