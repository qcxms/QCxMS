# This file is part of QCxMS.

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
  if fc.version().version_compare('<8')
    error('GCC version 8 or higher is required.')
  endif
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
  add_project_arguments(
    '-DLINUX',
    language: 'c',
  )
  if fc.version().version_compare('<19.1')
    error('Intel version 19.1 or higher is required.')
  endif
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
  cc = fc
  if add_languages('c', required: false, native: false)
    cc = meson.get_compiler('c')
  endif
  if fc_id == 'intel'
    mkl_dep += cc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
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

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.3',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static', 'json=disabled'],
)
exe_deps += mctc_dep

# Create the tight-binding framework as subproject
tblite_dep = dependency(
  'tblite',
  version: '>=0.2',
  fallback: ['tblite', 'tblite_dep'],
  default_options: ['default_library=static', 'api=false'],
)
exe_deps += tblite_dep

if tblite_dep.version().version_compare('>=0.3')
  error('tblite version 0.3 or higher is not supported yet.')
endif
 
# Create the rmsd tool as subproject
rmsd_dep = dependency(
  'rmsd-tool',
  version: '>=0.1',
  fallback: ['rmsd-tool', 'rmsd_dep'],
  default_options: ['default_library=static'],
)
exe_deps += rmsd_dep
