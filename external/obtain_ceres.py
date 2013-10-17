import obtain
import os
import sys
import shutil


# Download GLOG
obtain.obtain('GLOG v0.3.3', 'glog-0.3.3', 'http://google-glog.googlecode.com/files/glog-0.3.3.tar.gz', 'glog-0.3.3.tar.gz')

if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'glog'):
	if os.name == 'nt':
		print '<<WARNING>> no automatic installation for this package'
	else:
		print '<<INSTALL>> configure and build GLOG v0.3.3'
		obtain.configure_build('glog-0.3.3')
	#end
else:
	print '<<INSTALL>> GLOG already installed'
#end

# Download Eigen3
execfile('obtain_eigen.py')
 
# Download CERES
obtain.obtain('CERES v1.7.0', 'ceres-solver-1.7.0', 'http://ceres-solver.googlecode.com/files/ceres-solver-1.7.0.tar.gz', 'ceres-solver-1.7.0.tar.gz')
print '<<WARNING>> CERES installation requires CMake. You need to run it yourself'

if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'ceres'):
	if os.name == 'nt':
		print '<<WARNING>> no automatic installation for this package'
	else:
		print '<<INSTALL>> configure and build CERES'
		os.chdir('.' + os.sep + 'ceres-solver-1.7.0')
		build_dir = os.pardir + os.sep + 'build' + os.sep
		libname = ''
		if os.name == 'posix':
			libname = 'libglog.a'
		else:
			libname = 'glog.lib'
		ret = os.system('cmake -DGLOG_LIB=' + build_dir + 'lib' + os.sep + libname + ' -DGLOG_INCLUDE=' + build_dir + 'include -DGFLAGS=OFF -DEIGEN_INCLUDE=' + build_dir + 'include -DCMAKE_INSTALL_PREFIX=' + build_dir + ' .')
		ret = os.system('make install')
		os.chdir(os.pardir)
	#end
else:
	print '<<INSTALL>> CERES already installed'
#end