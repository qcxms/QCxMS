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
  error('This project does not support compilation with GCC yet')
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

if get_option('tblite')
  # Create the tight-binding framework as subproject
  tblite_prj = subproject(
    'tblite',
    version: '>=0.1',
    default_options: [
      'default_library=static',
    ],
  )
  tblite_dep = tblite_prj.get_variable('tblite_dep')
  exe_deps += tblite_dep
  exe_deps += declare_dependency(
    compile_args: '-DWITH_TBLITE=1',
  )

  if install
    install_data(
      tblite_prj.get_variable('tblite_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'tblite'
    )
  endif
endif
