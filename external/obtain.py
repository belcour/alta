import urllib
import os
import sys
import shutil
import tarfile
import hashlib
import SCons.Errors
import SCons.Warnings as W
import SCons.SConf as C
import mimetypes
from subprocess import Popen, PIPE
from urllib.request import urlretrieve

# Warning class to report stuff about our dependencies.
class AltaDependencyWarning(W.WarningOnByDefault):
   pass

# Check if the build dir exists
if not os.path.exists('build'):
   os.mkdir('build')

# Make sure that FILENAME has the given SHA256, and raise an error
# otherwise.
def check_integrity(filename, sha256):
  hasher = hashlib.sha256()
  BLOCK_SIZE = 65536
  with open(filename, 'rb') as input:
          buf = input.read(BLOCK_SIZE);
          while len(buf) > 0:
                  hasher.update(buf)
                  buf = input.read(BLOCK_SIZE)
  obtained = hasher.hexdigest()
  if obtained != sha256:
          # Keep the faulty file around but move it out of the
          # way so we don't end up using it by mistake.
          os.rename(filename, filename + ".wrong-hash")

          print >> sys.stderr, "error: downloaded file '" + \
                  filename + "' is inauthentic"
          print >> sys.stderr, \
                  "error: got sha256 hash {0} but expected {1}".format(obtained, sha256)
          raise SCons.Errors.BuildError(errstr = "downloaded file is inauthentic",
                                        filename = filename,
                                        exitstatus = 1)

# Download from URL to FILENAME.
# Return True if the package was successfully downloaded, else return False.
def download(url, filename):
   if os.path.exists(filename):
      C.progress_display('Package already downloaded')
      return True

   try:
      urllib.request.urlretrieve(url, filename)
      return True
   except IOError as e:
      # In Python 2.x, 'urlretrieve' raises IOError upon HTTP errors,
      # except for some errors such as 404.  Oh well.
      W.warn(AltaDependencyWarning,
             "{0}: download failed: {1}".format(url, e.strerror))
      return False

# Uncompress the archive
def uncompress(filename):
   if mimetypes.guess_type(filename)[1] == 'gzip':
      tfile = tarfile.open(filename, 'r:gz')
      tfile.extractall()

# Apply a patch to some file
def patch(filename, patch):
   ret = os.system('patch ' + filename + ' ' + patch)

# Obtain a package, check its integrity, and uncompress it.
def obtain(name, rep, url, filename, sha256):
   if not os.path.exists(rep):
      C.progress_display('obtaining ' + name)

      # Try to download the file if it exist. If an error happens, return
      # False.
      if not download(url, filename):
         return False

      check_integrity(filename, sha256)
      uncompress(filename)
      #if os.path(filename).splitext in ['.tar', '.gz', '.tgz', '.zip']:
      #    uncompress(filename)
      return True

   else:
      C.progress_display(name + ' source code is already available')
      return True

# Launch './configure' and 'make install' for UNIX like archives.
# The command is runned silent. TODO: Add the output of the command
# to the configuration file.
def configure_build(rep, options = []):
   os.chdir(rep)

   args = ['./configure', '-q', '--prefix=' + os.getcwd() + os.sep + os.pardir + os.sep + 'build']
   #call = Popen(args + options.split(), stdout=PIPE, stderr=PIPE, stdin=PIPE)
   call = Popen(args + options)
   ret  = call.wait()
   if ret != 0:
      C.progress_display('Warning: unable to configure package' + rep)
      os.chdir(os.pardir)
      return False

   #call = Popen(['make', 'install'], stdout=PIPE, stderr=PIPE, stdin=PIPE)
   call = Popen(['make', 'install'])
   ret  = call.wait()
   if ret != 0:
      C.progress_display('Warning: unable to build & install package' + rep)
      os.chdir(os.pardir)
      return False

   os.chdir(os.pardir)
   return True

# Launch 'cmake' and 'make && make install' (or equivalent) for
# package 'rep'. Default options are to build static libs and to
# install then in the #external/build/ directory.
def cmake_build(rep, options = ''):

   os.chdir(rep)

   # Construct the CMake command
   build_dir = os.pardir + os.sep + 'build' + os.sep
   cmake_cmd = ['cmake', '-DBUILD_SHARED_LIBS=OFF',
                '-DCMAKE_INSTALL_PREFIX=' + build_dir,
                '-DCMAKE_BUILD_TYPE=Release']
   if os.name == 'nt' :
      cmake_cmd = cmake_cmd + ['-G', 'NMake Makefiles']

   try:
      ret = Popen(cmake_cmd + options.split() + ['.']).wait()
      if ret != 0:
         C.progress_display('Warning: unable to configure package' + rep)
         os.chdir(os.pardir)
         return False
   except OSError as e:
      W.warn(AltaDependencyWarning, "execution of '{0}' failed: {1}".format(cmake_cmd, e.strerror))
      os.chdir(os.pardir)
      return False

   # Construct the build command and run it
   build_cmd = []
   if os.name == 'nt':
      build_cmd = build_cmd + ['nmake']
   else:
      build_cmd = build_cmd + ['make']

   ret = Popen(build_cmd + ['install']).wait()
   if ret != 0:
      C.progress_display('Warning: unable to build & install package ' + rep)
      os.chdir(os.pardir)
      return False

   os.chdir(os.pardir)
   return True
