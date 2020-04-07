project('timeseries_process','fortran')
add_global_arguments('-cpp',language:'fortran')
add_global_arguments('-ffixed-line-length-none',language:'fortran')

################################################################################
# host system detection

host_os = host_machine.system().to_lower()
message('Host os: ' + host_os)

platform_linux = (
  host_os.contains('linux')
)

platform_windows = (
  host_os.contains('mingw') or
  host_os.contains('cygwin') or
  host_os.contains('windows')
)

platform_osx = (
  host_os.contains('machten') or
  host_os.contains('rhapsody') or
  host_os.contains('darwin')
)

#
################################################################################

if platform_windows
add_global_arguments('-DPLATFORM_IS_WINDOWS',language:'fortran')
add_global_link_arguments('-static',language:'fortran')
endif

executable('bern2time','src/bern2time.for','src/bubblesort.for',install : true)
executable('getpl','src/getpl.for','src/bubblesort.for',install : true)
executable('comfilt','src/comfilt.for','src/bubblesort.for','src/int2char.for','src/lsqr.for',install : true)
executable('raw','src/comfilt.for','src/bubblesort.for','src/int2char.for','src/lsqr.for',fortran_args : '-DRAW',install : true)
executable('comfilt_trend','src/comfilt_trend.for','src/bubblesort.for','src/int2char.for','src/lsqr.for',install : true)
executable('remove_trend','src/remove_trend.for',install : true)
executable('remove_period','src/remove_period.for',install : true)
executable('remove_antenna','src/remove_antenna.for',install : true)
executable('fit','src/fit.for','src/int2char.for','src/lsqr.for',install : true)
executable('plot_ts','src/plot_ts.for',install : true)
