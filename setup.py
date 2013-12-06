from numpy.distutils.core import setup, Extension
from numpy import get_include
import numpy.version
from distutils.version import StrictVersion

import generatePyCamb
import os.path
import sys
if '--nonstop' in sys.argv:
  sys.argv.remove('--nonstop')
  from nonstopf2py import f2py
else:
  from numpy import f2py


cambsources = ['camb/%s' % f for f in [
    'constants.f90',
    'utils.F90',
    'subroutines.f90',
    'inifile.f90',
    'power_tilt.f90',
    'recfast.f90',
    'reionization.f90',
    'modules.f90',
    'bessels.f90',
    'equations.f90',
    'halofit.f90',
    'lensing.f90',
    'SeparableBispectrum.F90',
    'cmbmain.f90',
    'camb.f90',
]]

for f in cambsources:
  if not os.path.exists(f):
    raise Exception("At least one of CAMB code file: '%s' is not found. Download and extract to camb/.  You can use teh extract_camb.sh script to help" % f)

try: os.mkdir('src')
except: pass
generatePyCamb.main()

f2py.run_main(['-m', '_pycamb', '-h', '--overwrite-signature', 'src/py_camb_wrap.pyf', 
         'src/py_camb_wrap.f90', 'skip:', 'makeparameters', ':'])

# Newer versions of f2py (from numpy >= 1.6.2) use specific f90 compile args
if StrictVersion(numpy.version.version) > StrictVersion('1.6.1'):
    pycamb_ext = Extension("pycamb._pycamb",
                           ['src/py_camb_wrap.pyf'] + cambsources + ['src/py_camb_wrap.f90'],
                           extra_f90_compile_args=['-O0', '-g', '-Dintp=npy_intp', '-fopenmp'],
                           libraries=['gomp'],
                           include_dirs=[get_include()],
                           )
else:
    Extension("pycamb._pycamb",
             ['src/py_camb_wrap.pyf'] + cambsources + ['src/py_camb_wrap.f90'],
             extra_compile_args=['-O0', '-g', '-Dintp=npy_intp'],
             include_dirs=[get_include()],
             )

# Perform setup
setup(name="pycamb", version="0.3",
      author="Joe Zuntz",
      author_email="joezuntz@googlemail.com",
      description="python binding of camb, you need sign agreement and obtain camb source code to build this. Thus we can not GPL this code.",
      url="http://github.com/joezuntz/pycamb",
      zip_safe=False,
      install_requires=['numpy'],
      requires=['numpy'],
      packages=[ 'pycamb' ],
      package_dir={'pycamb': 'src'},
      data_files=[('pycamb/camb', ['camb/HighLExtrapTemplate_lenspotentialCls.dat'])],
      scripts=[],
      ext_modules=[pycamb_ext]
    )

