from numpy.distutils.core import setup, Extension
from numpy import get_include
import generatePyCamb
import os.path
from nonstopf2py import f2py

# Get CAMB from http://camb.info, untar and copy *.[fF]90 to src/
# this is done by the script extract_camb.sh

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
  if not os.path.exists('camb/Makefile'):
    raise Exception("At least one of CAMB code file: '%s' is not found. Download and extract to camb/" % f)

try: os.mkdir('src')
except: pass
generatePyCamb.main()

f2py.run_main(['-m', '_pycamb', '-h', '--overwrite-signature', 'src/py_camb_wrap.pyf', 
         'src/py_camb_wrap.f90', 'skip:', 'makeparameters', ':'])

setup(name="pycamb", version="0.1",
      author="Joe Zuntz",
      author_email="jaz@astro.ox.ac.uk",
      description="python binding of camb, you need sign agreement and obtain camb source code to build this. Thus we can not GPL this code.",
      url="http://github.com/rainwoodman/#",
      download_url="http://web.phys.cmu.edu/~yfeng1/#",
      zip_safe=False,
      install_requires=['numpy'],
      requires=['numpy'],
      packages = [ 'pycamb' ],
      package_dir = {'pycamb': 'src'},
      data_files = [('pycamb/camb', ['camb/HighLExtrapTemplate_lenspotentialCls.dat'])],
      scripts = [],
      ext_modules = [
        Extension("pycamb._pycamb", 
             ['src/py_camb_wrap.pyf'] + cambsources +['src/py_camb_wrap.f90'],
             extra_compile_args=['-O0', '-g', '-Dintp=npy_intp'],
#             libraries = {'noexit': ['src/noexit.c']},
             include_dirs=[get_include()],
        )]
    )

