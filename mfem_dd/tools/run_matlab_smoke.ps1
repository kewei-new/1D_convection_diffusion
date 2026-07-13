$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
matlab -batch "cd('$root/matlab'); startup_mfem_dd; test_smoke"
