import obtain
import os
import sys
import shutil
import subprocess
import SCons.Warnings as W
import SCons.SConf as C

# Download GLOG
version   = '0.3.3'
base      = 'glog'
name      = 'GLOG v' + version
directory = base + '-' + version
url       = 'http://google-glog.googlecode.com/files/glog-' + version + '.tar.gz'
filename  = 'glog-' + version + '.tar.gz'
sha256    = 'fbf90c2285ba0561db7a40f8a4eefb9aa963e7d399bd450363e959929fe849d0'

glog_obtained = False
glog_compiled = False
if os.name != 'nt':
   glog_obtained = obtain.obtain(name, directory, url, filename, sha256)
   glog_compiled = os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'glog')

if glog_obtained and not glog_compiled:

   if sys.platform == 'darwin':
      obtain.patch('glog-0.3.3/src/glog/stl_logging.h.in', 'glog.patch')

   C.progress_display('configuring and building GLOG for CERES')
   glog_compiled = obtain.configure_build(directory,
                                          ['--disable-shared',
                                           '--enable-static', '--with-pic'])


# Download CERES.  Assume Eigen is already available.
version   = '1.9.0'
base      = 'ceres-solver'
name      = 'CERES v' + version
directory = base + '-' + version
url       = 'http://ceres-solver.org/ceres-solver-' + version + '.tar.gz'
filename  = 'ceres-solver-' + version + '.tar.gz'
sha256    = '30ac0729249f908afe80cb6fd06ae6d037f25a60d9fac54f61344389adab9c1a'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

## Test for the presence of already compiled ceres version in
## the $ALTA/external/build directory. Then test for the
## presence of cmake.
compile_test = not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'ceres')

if obtained and compile_test:
   C.progress_display('configuring and building ' + name)

   # Build cmake script
   build_dir = os.pardir + os.sep + 'build' + os.sep
   options   = ''
   if glog_compiled:
      options = '-DGLOG_LIBRARY=' + build_dir + 'lib' + os.sep + 'libglog.a' + ' -DGLOG_INCLUDE_DIR=' + build_dir + 'include' + ' -DMINIGLOG=OFF'
   else:
      options = '-DMINIGLOG=ON'
   options = options + ' include -DGFLAGS=OFF -DEIGEN_INCLUDE_DIR=' + build_dir + os.sep + 'include' + ' -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF'

   obtain.cmake_build(directory, options)

else:
   W.warn(obtain.AltaDependencyWarning,
             'CERES already installed or cannot be installed automatically')
