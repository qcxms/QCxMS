# This file is part of QCxMS.

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
static = get_option('default_library') == 'static' and get_option('fortran_link_args').contains('-static')

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
  add_project_arguments(
    '-DLINUX',
    language: 'c',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
  error('This project does not support compilation with NVHPC yet')
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  exe_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel'
    mkl_dep += fc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += fc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += fc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += fc.find_library('mkl_core')
  exe_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  exe_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas_dep')
  endif
  exe_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
    endif
    exe_deps += lapack_dep
  endif

else
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack')
  endif
  exe_deps += lapack_dep
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas')
  endif
  exe_deps += blas_dep
endif

# Look for tool-chain library as external dependency
mctc_dep = dependency('mctc-lib', required: false, static: static)
if not mctc_dep.found()
  # Create the tool chain library as subproject
  mctc_prj = subproject(
    'mctc-lib',
    version: '>=0.2',
    default_options: [
      'default_library=static',
    ],
  )
  mctc_dep = mctc_prj.get_variable('mctc_dep')

  if install
    install_data(
      mctc_prj.get_variable('mctc_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-lib'
    )
  endif
endif
exe_deps += mctc_dep

# Look tight-binding framework as external dependency
tblite_dep = dependency('tblite', required: false, static: static)
if not tblite_dep.found()
  # Create the tight-binding framework as subproject
  tblite_prj = subproject(
    'tblite',
    version: '>=0.1',
    default_options: [
      'default_library=static',
    ],
  )
  tblite_dep = tblite_prj.get_variable('tblite_dep')

  if install
    install_data(
      tblite_prj.get_variable('tblite_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'tblite'
    )
  endif
endif
exe_deps += tblite_dep
