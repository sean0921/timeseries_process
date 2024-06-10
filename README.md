- I won't update this repo anymore.  Actually I don't know how this program works either. (?

# 時間序列處理工具

- [點我直接下載 .exe 檔](https://github.com/sean0921/timeseries_process/releases)
- [點我直接查看 Fortran 程式](https://github.com/sean0921/timeseries_process/tree/master/src)

![](https://travis-ci.org/sean0921/timeseries_process.svg?branch=master)
![](https://i.imgur.com/gwfKEH5.png)

## 使用 GNU Fortran 與 Meson Build System 編譯

* Windows 平臺：

    - 在 Windows 7/10 底下安裝 [MSYS2 環境](https://www.msys2.org/) ([教學](https://magiclen.org/msys2/)) 編譯
    - 按壓 `開始鍵+R` 組合鍵，在「執行」視窗輸入 `C:\msys64\mingw32.exe` 後按確認，開啟 **MSYS2 MinGW 32-bit** 環境
        + *或是* 按壓 `開始鍵+R` 組合鍵，在「執行」視窗輸入 `C:\msys64\mingw64.exe` 後按確認，開啟 **MSYS2 MinGW 64-bit** 環境
        + 開始使用 MSYS2 後，建議先更新環境：`pacman -Syu`，更新後若有出現相關提示，請按照提示直接關閉視窗，再於打開後執行 `pacman -Su`
    - 使用 `pacman` 指令來安裝 Meson Build System (`mingw-w64-i686-meson`)、Ninja Build (`mingw-w64-i686-ninja`)、編譯 32 位元 Windows 執行檔用的 GNU Fortran (`mingw-w64-i686-gcc-fortran`) 以及 Git 版本控制系統(`git`)：
        + `pacman -S mingw-w64-i686-meson mingw-w64-i686-ninja mingw-w64-i686-gcc-fortran git`
        + *或是* 使用 `pacman` 指令來安裝 Meson Build System (`mingw-w64-x86_64-meson`)、Ninja Build (`mingw-w64-x86_64-ninja`)、編譯 64 位元 Windows 執行檔用的 GNU Fortran (`mingw-w64-x86_64-gcc-fortran`) 以及 Git 版本控制系統(`git`)： `pacman -S mingw-w64-x86_64-meson mingw-w64-x86_64-ninja mingw-w64-x86_64-gcc-fortran git`
    - 使用 `git` 指令抓取本專案所有程式碼與版本紀錄：
        + `git clone https://github.com/sean0921/timeseries_process.git timeseries_process`
        + `cd timeseries_process` 即可進入本專案程式碼最上層目錄
        + 若已經透過 `git` 抓取下來，在程式碼最上層目錄執行 `git pull` 即可抓取更新
    - 設定原生環境下編譯 Windows 程式時才會用到的環境變數：(也可以直接加進 MSYS2 的 `~/.bashrc` 檔末端)
        + `export WINDRES=windres`
    - 在程式碼最上層目錄輸入：
        + `meson builddir; cd builddir; ninja`
    - 之後系統就會開始編譯程式了，結束後輸入 `ls`，即可在當前目錄底下發現已經編譯好的 `*.exe` 程式
    - 再輸入 `start .`，使用 Windows 檔案總管開啟當前目錄，就可以把編譯好的程式複製到你要的地方

## 使用方式

* 執行順序
  1. `bern2time.exe` (Transform Bersese 5.x `FN*.OUT` file to E/N/U format plain text)
  2. `raw.exe`/`comfilt.exe`
  3. `remove_trend.exe`
  4. `remove_period.exe`
  5. `remove_antenna.exe`
  6. `fit.exe` (Fit Velocity)

* 其他待補充

## 參考文獻

* Wdowinski, S., Bock, Y., Zhang, J., Fang, P., and Genrich, J. ( 1997), Southern California permanent GPS geodetic array: Spatial filtering of daily positions for estimating coseismic and postseismic displacements induced by the 1992 Landers earthquake, J. Geophys. Res., 102( B8), 18057– 18070, doi:10.1029/97JB01378.
* Nikolaidis, Rosanne. (2002). Observation of Geodetic and Seismic Deformation with the Global Positioning System.

---

# Time Series Processing Tools

## Compiling by GNU Fortran and Meson Build Systerm

* Install Meson Build System, Ninja Build, GNU Fortran  for your platform (Linux/Windows MSYS2)  first.
* `meson builddir; cd builddir; ninja`

## Reference for order

1. `bern2time`
2. `raw`/`comfilt`
3. `remove_trend`
4. `remove_period`
5. `remove_antenna`
6. `Fit_velocity`

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
