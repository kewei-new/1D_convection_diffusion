param(
  [string]$MfemDir = $env:MFEM_DIR,
  [string]$Generator = "Ninja",
  [string]$CxxCompiler = $env:CXX
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$build = Join-Path $root "build/mfem_dd"
if ([string]::IsNullOrWhiteSpace($MfemDir)) {
  $MfemDir = "E:/Projects/codes/2D_convection_diffusion/DD_model/模拟/V6/build"
}
if ([string]::IsNullOrWhiteSpace($CxxCompiler)) {
  $CxxCompiler = "D:/mingw64/bin/g++.exe"
}
cmake -G "$Generator" -S (Join-Path $root "cpp") -B $build -DMFEM_DIR="$MfemDir" -DCMAKE_CXX_COMPILER="$CxxCompiler"
cmake --build $build
ctest --test-dir $build --output-on-failure
