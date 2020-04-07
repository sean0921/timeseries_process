# 時間序列處理工具

![](https://travis-ci.org/sean0921/timeseries_process.svg?branch=master)

## 使用 GNU Fortran 與 Meson Build System 編譯

* Windows 平臺：

    - 在 Windows 7/10 底下安裝 [MSYS2 環境](https://www.msys2.org/) ([教學](https://magiclen.org/msys2/)) 編譯
    - 安裝 Meson Build System(`meson`), Ninja Build(`ninja`), 和 GNU Fortran(`gfortran`): `pacman -S meson ninja gfortran`
    - 在程式碼最上層目錄輸入 `meson builddir; cd builddir; ninja`，就會開始編譯了，你會在 `builddir/` 資料夾底下發現編譯好的 `*.exe` 程式

## 使用方式

* 執行順序
  1. `bern2time.exe`
  2. `raw.exe`
  3. `remove_trend.exe`
  4. `remove_period.exe`
  5. `remove_antenna.exe`
  6. `Fit_velocity.exe`

* 其他待補充

## 參考文獻

* Wdowinski, S., Bock, Y., Zhang, J., Fang, P., and Genrich, J. ( 1997), Southern California permanent GPS geodetic array: Spatial filtering of daily positions for estimating coseismic and postseismic displacements induced by the 1992 Landers earthquake, J. Geophys. Res., 102( B8), 18057– 18070, doi:10.1029/97JB01378.
* Nikolaidis, Rosanne. (2002). Observation of Geodetic and Seismic Deformation with the Global Positioning System.

---

# Time Series Processing Tools

![](https://travis-ci.org/sean0921/timeseries_process.svg?branch=master)

## Compiling by GNU Fortran and Meson Build Systerm

* Install Meson Build System, Ninja Build, GNU Fortran  for your platform (Linux/Windows MSYS2)  first.
* `meson builddir; cd builddir; ninja`

## Reference for order

0. `bern2time.exe`
1. `raw.exe`
2. `remove_trend`
3. `remove_period`
4. `remove_antenna`
5. `Fit_velocity`

## Input file

* `bern2time`:
    - `getpl.gout` (from `getpl`)
    - `sta-file` (from `getpl`)

* `comfilt` (`raw`)
    - `comfilt.inp` (`raw.inp`, adapted from `sta-file`, format example:`comfilt.inp.example`/`raw.inp.example`)

* `remove_{trend,period}`
    - `remove.inp` (copied from `comfilt.inp`/`raw.inp`)

* `remove_attenna`:
    - `remove.inp` (copied from `comfilt.inp`/`raw.inp`)

* `fit`:
    - `fit.inp` (copied from `comfilt.inp`/`raw.inp`)

## References/Related Articles

* Wdowinski, S., Bock, Y., Zhang, J., Fang, P., and Genrich, J. ( 1997), Southern California permanent GPS geodetic array: Spatial filtering of daily positions for estimating coseismic and postseismic displacements induced by the 1992 Landers earthquake, J. Geophys. Res., 102( B8), 18057– 18070, doi:10.1029/97JB01378.
* Nikolaidis, Rosanne. (2002). Observation of Geodetic and Seismic Deformation with the Global Positioning System.
