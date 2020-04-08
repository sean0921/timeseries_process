project('timeseries_process',
        'fortran',
        version: '2020.408.0',
        default_options:['buildtype=debugoptimized']
)
add_global_arguments('-cpp',language:'fortran')
add_global_arguments('-ffixed-line-length-none',language:'fortran')

################################################################################
# host system detection

host_os = host_machine.system().to_lower()
message('Host os: ' + host_os)

platform_windows = (
  host_os.contains('mingw') or
  host_os.contains('cygwin') or
  host_os.contains('windows')
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

####### Add WINDRES module commands #######

if platform_windows

rcinfo_windows = import('windows')

bern2time_rcinfo        = rcinfo_windows.compile_resources('windres/bern2time.rc'     )
getpl_rcinfo            = rcinfo_windows.compile_resources('windres/getpl.rc'         )
comfilt_rcinfo          = rcinfo_windows.compile_resources('windres/comfilt.rc'       )
raw_rcinfo              = rcinfo_windows.compile_resources('windres/raw.rc'           )
comfilt_trend_rcinfo    = rcinfo_windows.compile_resources('windres/comfilt_trend.rc' )
remove_trend_rcinfo     = rcinfo_windows.compile_resources('windres/remove_trend.rc'  )
remove_period_rcinfo    = rcinfo_windows.compile_resources('windres/remove_period.rc' )
remove_antenna_rcinfo   = rcinfo_windows.compile_resources('windres/remove_antenna.rc')
fit_rcinfo              = rcinfo_windows.compile_resources('windres/fit.rc'           )
plot_ts_rcinfo          = rcinfo_windows.compile_resources('windres/plot_ts.rc'       )

bern2time_src        += bern2time_rcinfo
getpl_src            += getpl_rcinfo
comfilt_src          += comfilt_rcinfo
raw_src              += raw_rcinfo
comfilt_trend_src    += comfilt_trend_rcinfo
remove_trend_src     += remove_trend_rcinfo
remove_period_src    += remove_period_rcinfo
remove_antenna_src   += remove_antenna_rcinfo
fit_src              += fit_rcinfo
plot_ts_src          += plot_ts_rcinfo

endif

###########################################

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
