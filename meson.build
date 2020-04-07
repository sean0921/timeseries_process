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

bern2time_src       = ['src/bern2time.for','src/bubblesort.for']
getpl_src           = ['src/getpl.for','src/bubblesort.for']
comfilt_src         = ['src/comfilt.for','src/bubblesort.for','src/int2char.for','src/lsqr.for']
raw_src             = ['src/comfilt.for','src/bubblesort.for','src/int2char.for','src/lsqr.for']
comfilt_trend_src   = ['src/comfilt_trend.for','src/bubblesort.for','src/int2char.for','src/lsqr.for']
remove_trend_src    = ['src/remove_trend.for']
remove_period_src   = ['src/remove_period.for']
remove_antenna_src  = ['src/remove_antenna.for']
fit_src             = ['src/fit.for','src/int2char.for','src/lsqr.for']
plot_ts_src         = ['src/plot_ts.for']

#### TODO: Add WINDRES module commands ####

#if platform_windows
#bern2time_src       = ['']
#getpl_src           = ['']
#comfilt_src         = ['']
#raw_src             = ['']
#comfilt_trend_src   = ['']
#remove_trend_src    = ['']
#remove_period_src   = ['']
#remove_antenna_src  = ['']
#fit_src             = ['']
#plot_ts_src         = ['']
#endif

############################################

executable('bern2time'      , bern2time_src     )
executable('getpl'          , getpl_src         )
executable('comfilt'        , comfilt_src       )
executable('raw'            , raw_src           , fortran_args : '-DRAW')
executable('comfilt_trend'  , comfilt_trend_src )
executable('remove_trend'   , remove_trend_src  )
executable('remove_period'  , remove_period_src )
executable('remove_antenna' , remove_antenna_src)
executable('fit'            , fit_src           )
executable('plot_ts'        , plot_ts_src       )
