import obtain
import os
import sys
import shutil
import SCons.Errors as Errors

# Download Eigen3
version   = '3.3.9'
base      = 'Eigen'
name      = 'Eigen v' + version
directory = 'eigen-' + version
url       = 'https://gitlab.com/libeigen/eigen/-/archive/' + version + '/eigen-' + version + '.tar.gz'
filename  = 'eigen-' + version + '.tar.gz'
sha256    = '7985975b787340124786f092b3a07d594b2e9cd53bbfe5f3d9b1daee7d55f56f'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

if obtained:
   rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
   if not os.path.exists(rep):
      shutil.copytree(directory + os.sep + 'Eigen', rep)

   unsup_rep = 'build' + os.sep + 'include' + os.sep + 'unsupported'
   if not os.path.exists(unsup_rep):
      shutil.copytree(directory + os.sep + 'unsupported', unsup_rep)

else:
   msg = "download from '" + url + "' failed"
   raise Errors.BuildError(errstr = msg,
                           filename = filename,
                           exitstatus = 1)
