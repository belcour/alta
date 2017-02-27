# ALTA --- Analysis of Bidirectional Reflectance Distribution Functions
#
# Copyright (C) 2014, 2015 CNRS
# Copyright (C) 2013, 2014, 2015, 2016 Inria
# Copyright (C) 2015 Universite de Montreal
#
# This file is part of ALTA.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0.  If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import SCons.SConf as C

# Arrange so that the 'obtain' module can always be imported.
top_srcdir = Dir('.').srcnode().abspath
sys.path += [ top_srcdir + '/external' ]

## Add ALTA custom cmd configurations
##
AddOption('--cfg', help='Specify a configuration file')


## Import configuration from a config file
##
configFile = GetOption('cfg')
if configFile == None:

	if sys.platform == 'win32':
		configFile = "./configs/scons/config-windows-cl.py"
	elif sys.platform == 'darwin':
		configFile = "./configs/scons/config-macos-clang.py"
	elif sys.platform == 'linux2':
		configFile = "./configs/scons/config-linux-gcc.py"
	else:
		print '<<ERROR>> You need to specify a configuration file using:'
		print '<<ERROR>>    scons --cfg=[filename]'
		print '<<ERROR>> Please find example of configuration files in ${ALTA}/configs/scons/'
		Exit(1)
#end

if not os.path.exists(configFile):
	print '<<ERROR>> the config file you specified \"' + configFile + '\" does not exists'
	Exit(1)
else:
	print '<<INFO>> Using config file \"' + configFile + '\"'
#end

# List of C++ compilers we look for by default.
cxx_compilers = ['g++', 'c++', 'clang++', 'cl' ]
if(os.environ.has_key('CXX')):
   cxx_compilers = [os.environ['CXX'], cxx_compilers]

# List of Python 2.x interpreter names to look for.
python_programs = [ 'python2.7', 'python2', 'python', sys.executable ]

def program_file_name(choices):
  """
  Look for the programs listed in CHOICES.  Return the absolute file
  name of one that matches, or the first element of CHOICES.

  """
  for program in choices:
    full = WhereIs(program)
    if full != None:
      return full
  return choices[0]

vars = Variables(configFile)

vars.Add('INSTALL_PREFIX',    'Parent installation directory',
         default = '/usr/local')
vars.Add('CXX',               'C++ compiler',
         default = program_file_name(cxx_compilers))
vars.Add('CCFLAGS',           'Compiler\'s flags',
         default = '-std=c++11 -g -Wall')
vars.Add('LINKFLAGS',         'Linker\'s flags',
         default = '')
vars.Add('EIGEN_INC',         'Eigen include directory (mandatory)')
vars.Add('PYTHON_INC',        'Python and boost-python include directory')
vars.Add('PYTHON_DIR',        'Python and boost-python libraries directory')
vars.Add('PYTHON_LIB',        'Python and boost-python libraries', default = [])
vars.Add('OPENEXR_INC',       'OpenEXR include directory', default=[])
vars.Add('OPENEXR_LIB',       'OpenEXR libraries', default = [])
vars.Add('OPENEXR_DIR',       'OpenEXR libraries directory', default=[])
vars.Add('FLANN_INC',         'FLANN include directory')
vars.Add('FLANN_DIR',         'FLANN libraries directory')
vars.Add('FLANN_LIB',         'FLANN libraries')
vars.Add('CGAL_INC',          'CGAL include directory', default = [])
vars.Add('CGAL_DIR',          'CGAL libraries directory', default = [])
vars.Add('CGAL_LIB',          'CGAL libraries', default = [])
vars.Add('OPENMP_CFLAGS',     'OpenMP compiler flags', default = None)
vars.Add('OPENMP_LINKFLAGS',  'OpenMP linker flags', default = None)
vars.Add('QUADPROG_INC',      'QUADPROG include directory')
vars.Add('QUADPROG_DIR',      'QUADPROG libraries directory')
vars.Add('QUADPROG_LIB',      'QUADPROG libraries')
vars.Add('CERES_INC',         'CERES include directory')
vars.Add('CERES_DIR',         'CERES libraries directory')
vars.Add('CERES_LIB',         'CERES libraries')
vars.Add('CERES_OPT_LIB',     'CERES optional libraries')
vars.Add('NLOPT_INC',         'NLOPT include directory')
vars.Add('NLOPT_DIR',         'NLOPT libraries directory')
vars.Add('NLOPT_LIB',         'NLOPT libraries', default = [])
vars.Add('NLOPT_OPT_LIB',     'NLOPT optional libraries', default = [])
vars.Add('IPOPT_INC',         'IPOPT include directory')
vars.Add('IPOPT_DIR',         'IPOPT libraries directory')
vars.Add('IPOPT_LIB',         'IPOPT libraries', default = [])
vars.Add('IPOPT_OPT_LIB',     'IPOPT optional libraries')
vars.Add('MATLAB_INC',        'MATLAB include directory')
vars.Add('MATLAB_DIR',        'MATLAB directory')
vars.Add('MATLAB_LIB',        'MATLAB libraries')
vars.Add('PYTHON',            'Python interpreter',
         default = program_file_name(python_programs))


##
# Copy the system environment.
#
# Update the PKG_CONFIG_PATH variable to add the package configuration
# files contained is #external/build/lib/pkgconfig
##
envVars = {}
for var in [ 'PATH', 'CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH',
             'PKG_CONFIG_PATH', 'TMP', 'TMPDIR' ]:
	if var in os.environ:
		envVars[var] = os.environ[var]
	else:
		envVars[var] = '';

if len(envVars['PKG_CONFIG_PATH']) > 0:
	envVars['PKG_CONFIG_PATH'] += ':'
envVars['PKG_CONFIG_PATH'] += os.path.abspath('external' + os.sep + 'build' + os.sep + 'lib' + os.sep + 'pkgconfig')
env = Environment(variables = vars, ENV = envVars)


# Generate help text for the build variables.
Help(vars.GenerateHelpText(env))

# Rule to build the documentation.
if 'doc' in COMMAND_LINE_TARGETS:
  env.Alias('doc', env.Command(Dir('#documents/doxygen/html'),
                               '#documents/doxygen.conf',
                               'doxygen doxygen.conf',
                               chdir = Dir('#documents').srcnode().abspath))

C.progress_display('the current platform is: ' + env['PLATFORM'])


def CheckPKG(context, name):
		"""Return True if package NAME can be found with 'pkg-config'."""
		context.Message('Checking for %s using pkg-config... ' % name)
		ret = context.TryAction('pkg-config --exists \'%s\'' % name)[0]
		context.Result(ret)
		return ret


def library_available(env, pkgspec='', lib='', header='',
                      language='c++', inc_var='', lib_var=''):
  """Return True if the given library is available.

  First look for the LIB_VAR and INC_VAR construction variables,
  honoring them if they are defined.  Then look for PKGSPEC using
  pkg-config.  Last, try to build LIB with HEADER.  Configure ENV
  accordingly.

  """
  conf = Configure(env, custom_tests = { 'CheckPKG' : CheckPKG })

  # If a XXX_LIB is specified in the environment, add the various path
  # and link flags. Check if the library is correctly compiling and
  # linking with the header.
  if (lib in env) and (len(env[lib]) > 0):
    env.AppendUnique(LIBPATH = env[lib_var])
    env.AppendUnique(CPPPATH = env[inc_var])
    env.AppendUnique(LIBS = env[lib])

    # Check whether the library is usable.
    result = conf.CheckLibWithHeader(env[lib], header, language)
  elif conf.CheckPKG(pkgspec):
    env.ParseConfig('pkg-config --cflags --libs "' + pkgspec + '"')
    # XXX: We can't use 'CheckLibWithHeader' to verify that
    # everything's alright because we don't know the library name.  So
    # assume that pkg-config got it right.
    result = True
  else:
		result = conf.CheckLibWithHeader(pkgspec, header, language)

  conf.Finish()
  return result

def catch_available(env):
  """Return True if Catch is available."""
  env.AppendUnique(CPPPATH = '#external')
  conf = Configure(env)
  has_catch = conf.CheckCXXHeader('Catch/single_include/catch.hpp')
  conf.Finish()
  return has_catch

def openexr_available(env):
   """Return True if OpenEXR is available."""
   env.AppendUnique(CPPPATH = '#external')
   conf = Configure(env)
   has_tinyexr = conf.CheckCXXHeader('tinyexr/tinyexr.h')
   conf.Finish()
   return has_tinyexr
   #return library_available(env, pkgspec='OpenEXR',
   #						 inc_var='OPENEXR_INC',
   #						 lib_var='OPENEXR_DIR',
   #						 lib='OPENEXR_LIB',
   #						 header='ImfRgbaFile.h')

def CheckOpenMP(context):
  """
  Check whether the compiler supports '-fopenmp'.  If it does, add
  '-fopenmp' to the compilation and link flags.

  """
  context.Message("Checking whether '-fopenmp' is supported... ")
  env = context.env
  save_CCFLAGS = env['CCFLAGS']
  save_LINKFLAGS = env['LINKFLAGS']

  env['CCFLAGS'] = ' '.join([save_CCFLAGS, '-fopenmp'])
  env['LINKFLAGS'] = ' '.join([env['LINKFLAGS'], '-fopenmp'])

  # Users of Clang++ sometimes lack libgomp (which is the OpenMP
  # runtime that Clang uses), so better check that things link.
  has_fopenmp = conf.TryLink("""
#include <stdlib.h>
int frob (int x)
{
  return x * x;
}
int main (int argc, char *argv[])
{
  int z;
#pragma omp parallel for
  for (z = 0; z < atoi (argv[1]); z++)
    frob (z);
  return 0;
}
""", '.cpp')
  context.Result('yes' if has_fopenmp else 'no')

  env['CCFLAGS'] = save_CCFLAGS
  env['LINKFLAGS'] = save_LINKFLAGS
  if has_fopenmp:
    env['OPENMP_CFLAGS'] = ' -fopenmp'
    env['OPENMP_LINKFLAGS'] = ' -fopenmp'
  else:
    # Don't leave None in there.
    env['OPENMP_CFLAGS'] = ''
    env['OPENMP_LINKFLAGS'] = ''

  return has_fopenmp

def CheckCXX11(context):
  """
  Check whether the compiler supports a small part of C++11 extensions.
  """
  context.Message("Checking whether 'c++11' is supported... ")
  env = context.env
  has_cxx11 = conf.TryLink("""
#include <vector>
#include <memory>

struct A {
   A(int a) : _a(a) {}
   int _a;
};

int main(int argc, char* argv[]) {
   std::shared_ptr<A> myA(new A(1));
   std::vector<int> vec;
   vec.assign(10, myA->_a);
   for(auto& x : vec) {
      x += 10;
   }
   return 0;
}""", '.cpp')
  context.Result('yes' if has_cxx11 else 'no')

  return has_cxx11

# Export these for use in SConscripts.
Export('CheckPKG', 'library_available', 'openexr_available', 'catch_available')

conf = Configure(env, custom_tests = { 'CheckOpenMP': CheckOpenMP,
                                       'CheckCXX11' : CheckCXX11 })

# Determine the extra libraries that libcore (and thus everything
# else) depends on.  Plugins need to specify it in addition to -lcore
# because libcore is not a shared library.
ALTA_LIBS = []

# Libcore's uses 'dlopen', which is in libdl in GNU libc.
if conf.CheckLibWithHeader('dl', 'dlfcn.h', 'c++'):
		ALTA_LIBS = ['dl']

# Libcore uses 'clock_gettime', which is in librt in GNU libc.
if conf.CheckLibWithHeader('rt', 'sched.h', 'c++'):
		ALTA_LIBS = ALTA_LIBS + ['rt']

Export('ALTA_LIBS')

if not 'OPENMP_CFLAGS' in env:
  conf.CheckOpenMP()

if not conf.CheckCXX11():
   C.progress_display('Error: your compiler is not compatible with c++11 using the \'-std=c++11\' flag.')
   Exit(1);

env = conf.Finish()

## Load the configuration file if it exists. The configuration file
## is a python script that updates the env variable with different
## paths.
##
env.AppendUnique(LIBPATH = ['#external/build/lib'])
env.AppendUnique(LIBPATH = ['#build/core'])
env.AppendUnique(CPPPATH = ['#external/build/include'])
env.AppendUnique(CPPPATH = ['#sources'])

# Consider files changed as soon as their modification time changes.
env.Decider('timestamp-newer')

## Launch the compilations
##
Export('env')

def alta_sconscript(script):
  return env.SConscript(script,
                        variant_dir='build/' + os.path.dirname(script),
                        duplicate=0)

external = env.SConscript('external/SConscript')
core     = alta_sconscript('sources/core/SConscript')
plugins  = alta_sconscript('sources/plugins/SConscript')
softs    = alta_sconscript('sources/softs/SConscript')

if 'python' in COMMAND_LINE_TARGETS:
  python = alta_sconscript('sources/python/SConscript')
  env.Depends(python, core)

env.Depends(core, external)
env.Depends(plugins, core)
env.Depends(softs, core)

if 'tests' in COMMAND_LINE_TARGETS:
  tests = alta_sconscript('sources/tests/SConscript')
  env.Depends(tests, core)
  env.Depends(tests, plugins)
  if 'python' in COMMAND_LINE_TARGETS:
    env.Depends(tests, python)

env.Alias('install', env['INSTALL_PREFIX'])
