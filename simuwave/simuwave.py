
import os, sys 
sys.path.insert(0, os.path.abspath('./src'))
sys.path.insert(0, os.path.abspath('./src/fdtd'))
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui
import gui

__version__ = '1.0.0'


# Main program 	
def main():
	try: app = QtGui.QApplication.instance()
	except: app = QtGui.QApplication(sys.argv)	
	win_main = gui.MainWindow()
	win_main.show()		
	app.setActiveWindow(win_main)
	app.setApplicationName('SimuWave')
	app.setOrganizationName('Ghent University')
	app.setOrganizationDomain('ugent.be')
	app.setWindowIcon(QtGui.QIcon('./src/logo.svg'))		
	app.exec_()
	win_main.deleteLater()	
	app.deleteLater()
	sys.exit()	


# execute main program if this file is not imported
if __name__ == "__main__": main()
