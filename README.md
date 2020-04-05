# 時間序列處理工具

![](https://travis-ci.org/sean0921/timeseries_process.svg?branch=master)
![](https://i.imgur.com/gwfKEH5.png)

## 使用 GNU Fortran 與 Meson Build System 編譯

* Windows 平臺：

    - 在 Windows 7/10 底下安裝 [MSYS2 環境](https://www.msys2.org/) ([教學](https://magiclen.org/msys2/)) 編譯
    - 安裝 Meson Build System(`meson`), Ninja Build(`ninja`), 和 GNU Fortran(`gfortran`): `pacman -S meson ninja gfortran`
    - 在程式碼最上層目錄輸入 `meson builddir; cd builddir; ninja`，就會開始編譯了，你會在 `builddir/` 資料夾底下發現編號的程式，但 DLL 檔請自己抓或把 `C:\msys64\mingw32\bin` 加到你的 Windows 環境變數

## 使用方式

* 執行順序
  1. `bern2time`
  2. `raw`
  3. `remove_trend`
  4. `remove_period`
  5. `remove_antenna`
  6. `Fit_velocity`

* 其他待補充

## 參考文獻

* Wdowinski, S., Bock, Y., Zhang, J., Fang, P., and Genrich, J. ( 1997), Southern California permanent GPS geodetic array: Spatial filtering of daily positions for estimating coseismic and postseismic displacements induced by the 1992 Landers earthquake, J. Geophys. Res., 102( B8), 18057– 18070, doi:10.1029/97JB01378.
* Nikolaidis, Rosanne. (2002). Observation of Geodetic and Seismic Deformation with the Global Positioning System.

---

# Time Series Processing Tools

![](https://travis-ci.org/sean0921/timeseries_process.svg?branch=master)
![](https://i.imgur.com/gwfKEH5.png)

## Compiling by GNU Fortran and Meson Build Systerm

* Install Meson Build System, Ninja Build, GNU Fortran  for your platform (Linux/Windows MSYS2)  first.
* `meson builddir; cd builddir; ninja`

## Reference for order

0. bern2time
1. raw
2. remove_trend
3. remove_period
4. remove_antenna
5. Fit_velocity

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
