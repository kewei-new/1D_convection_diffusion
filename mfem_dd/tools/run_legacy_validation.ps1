param(
    [string]$BuildDir = "",
    [switch]$SkipCpp
)

$ErrorActionPreference = "Stop"
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$repoName = Split-Path $repoRoot -Leaf
$mfemRoot = Join-Path $repoRoot "mfem_dd"
$artifactDir = Join-Path $mfemRoot "artifacts\legacy_validation"
New-Item -ItemType Directory -Force -Path $artifactDir | Out-Null

if ([string]::IsNullOrWhiteSpace($BuildDir)) {
    $BuildDir = Join-Path $mfemRoot "build\mfem_dd"
}

function Convert-ToMatlabPath([string]$PathValue) {
    return $PathValue.Replace("\", "/")
}

function Invoke-Checked([scriptblock]$Command, [string]$Label) {
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

if ($repoName -like "1D_convection_diffusion") {
    $compareFunction = "mfemdd.compare_legacy_dd1d"
    $refineSteps = 5
    $order = 3
    $regressionFirstN = 4.507504048373722e-06
    $cppApp = "dd1d_mms"
    $cppArgs = @("-n", "20", "-o", "3", "-b", "legacy_baseline")
    $cppFirstN = 4.507504e-06
    $deviceCompareFunction = "mfemdd.compare_legacy_pn1d"
    $deviceForwardCurrent = 17981.0
    $deviceCqs = -24.835
    $nativeCompareFunction = "mfemdd.compare_legacy_dd1d"
} elseif ($repoName -like "2D_convection_diffusion") {
    $compareFunction = "mfemdd.compare_legacy_dd2d"
    $refineSteps = 4
    $order = 1
    $regressionFirstN = 4.507062201109525e-01
    $cppApp = "dd2d_mms"
    $cppArgs = @("-n", "6", "-o", "1", "-b", "legacy_baseline")
    $cppFirstN = 4.507062e-01
    $deviceCompareFunction = ""
    $deviceForwardCurrent = $null
    $deviceCqs = $null
    $nativeCompareFunction = ""
} else {
    throw "Unsupported repository for legacy validation: $repoName"
}

$matlabDir = Convert-ToMatlabPath (Join-Path $mfemRoot "matlab")
$compareJson = Convert-ToMatlabPath (Join-Path $artifactDir "legacy_runtime_compare.json")
$regressionDir = Convert-ToMatlabPath (Join-Path $artifactDir "matlab_regression")
$matlabCommand = "cd('$matlabDir'); startup_mfem_dd; " +
    "r=$compareFunction('backend','legacy_runtime','refine_steps',$refineSteps,'output_json','$compareJson'); " +
    "assert(r.matches_full_table); " +
    "s=run_regression_suite('quick',true,'order',$order,'backend','legacy_runtime','include_device',false,'output_dir','$regressionDir'); " +
    "assert(numel(s.rows)==2); assert(strcmp(char(s.rows(1).metrics.status),'legacy_runtime'));"
if (-not [string]::IsNullOrWhiteSpace($deviceCompareFunction)) {
    $deviceJson = Convert-ToMatlabPath (Join-Path $artifactDir "legacy_device_compare.json")
    $matlabCommand += " d=$deviceCompareFunction('backend','matlab_mfem','output_json','$deviceJson'); assert(d.matches_baseline);"
}
if (-not [string]::IsNullOrWhiteSpace($nativeCompareFunction)) {
    $nativeJson = Convert-ToMatlabPath (Join-Path $artifactDir "native_matlab_compare.json")
    $matlabCommand += " n=$nativeCompareFunction('backend','matlab_mfem','refine_steps',$refineSteps,'output_json','$nativeJson'); assert(n.matches_full_table);"
}

Invoke-Checked { matlab -batch $matlabCommand } "MATLAB legacy validation"

$regressionCsv = Join-Path (Join-Path $artifactDir "matlab_regression") "metrics.csv"
$regressionRows = Import-Csv $regressionCsv
$observedRegressionFirstN = [double]$regressionRows[0].n_l2_error
if ([Math]::Abs($observedRegressionFirstN - $regressionFirstN) -gt 1.0e-12) {
    throw "MATLAB regression first n_L2 mismatch: got $observedRegressionFirstN expected $regressionFirstN"
}

$matlabDeviceSummary = $null
if (-not [string]::IsNullOrWhiteSpace($deviceCompareFunction)) {
    $matlabDeviceSummary = [ordered]@{
        compare_json = (Join-Path $artifactDir "legacy_device_compare.json")
        iv_forward_current_1v = $deviceForwardCurrent
        cv_zero_bias_cqs = $deviceCqs
        status = "passed"
    }
}

$matlabNativeSummary = $null
if (-not [string]::IsNullOrWhiteSpace($nativeCompareFunction)) {
    $matlabNativeSummary = [ordered]@{
        compare_json = (Join-Path $artifactDir "native_matlab_compare.json")
        backend = "matlab_mfem"
        status = "passed"
    }
}

$cppSummary = $null
$cppDeviceSummary = $null
if (-not $SkipCpp) {
    if (-not (Test-Path (Join-Path $BuildDir "CMakeCache.txt"))) {
        throw "C++ build directory is missing or unconfigured: $BuildDir"
    }
    Invoke-Checked { cmake --build $BuildDir } "C++ build"
    Invoke-Checked { ctest --test-dir $BuildDir --output-on-failure } "CTest"

    $exe = Join-Path $BuildDir "$cppApp.exe"
    if (-not (Test-Path $exe)) {
        $exe = Join-Path $BuildDir $cppApp
    }
    if (-not (Test-Path $exe)) {
        throw "C++ executable not found: $cppApp"
    }

    $cppRunDir = Join-Path $artifactDir "cpp_baseline"
    New-Item -ItemType Directory -Force -Path $cppRunDir | Out-Null
    Push-Location $cppRunDir
    try {
        Invoke-Checked { & $exe @cppArgs } "C++ legacy baseline app"
    } finally {
        Pop-Location
    }

    $cppCsv = Join-Path $cppRunDir "metrics.csv"
    $cppRows = Import-Csv $cppCsv
    $observedCppFirstN = [double]$cppRows[0].n_l2_error
    if ([Math]::Abs($observedCppFirstN - $cppFirstN) -gt 1.0e-12) {
        throw "C++ baseline first n_L2 mismatch: got $observedCppFirstN expected $cppFirstN"
    }
    $cppSummary = [ordered]@{
        app = $cppApp
        metrics_csv = $cppCsv
        first_n_l2_error = $observedCppFirstN
        status = "passed"
    }

    if (-not [string]::IsNullOrWhiteSpace($deviceCompareFunction)) {
        $deviceExe = Join-Path $BuildDir "dd_device.exe"
        if (-not (Test-Path $deviceExe)) {
            $deviceExe = Join-Path $BuildDir "dd_device"
        }
        if (-not (Test-Path $deviceExe)) {
            throw "C++ executable not found: dd_device"
        }

        $cppDeviceRunDir = Join-Path $artifactDir "cpp_device_baseline"
        New-Item -ItemType Directory -Force -Path $cppDeviceRunDir | Out-Null
        Push-Location $cppDeviceRunDir
        try {
            Invoke-Checked { & $deviceExe -n 16 -o 1 -b legacy_baseline } "C++ device legacy baseline app"
        } finally {
            Pop-Location
        }

        $cppDeviceCsv = Join-Path $cppDeviceRunDir "metrics.csv"
        $cppDeviceRows = Import-Csv $cppDeviceCsv
        $observedForwardCurrent = [double]$cppDeviceRows[0].iv_forward_current_1v
        $observedCqs = [double]$cppDeviceRows[0].cv_zero_bias_cqs
        if ([Math]::Abs($observedForwardCurrent - $deviceForwardCurrent) -gt 1.0e-12) {
            throw "C++ device forward-current mismatch: got $observedForwardCurrent expected $deviceForwardCurrent"
        }
        if ([Math]::Abs($observedCqs - $deviceCqs) -gt 1.0e-12) {
            throw "C++ device Cqs mismatch: got $observedCqs expected $deviceCqs"
        }
        $cppDeviceSummary = [ordered]@{
            app = "dd_device"
            metrics_csv = $cppDeviceCsv
            iv_forward_current_1v = $observedForwardCurrent
            cv_zero_bias_cqs = $observedCqs
            status = "passed"
        }
    }
}

$summary = [ordered]@{
    repo = $repoName
    generated_at = (Get-Date).ToString("s")
    matlab_legacy_runtime = [ordered]@{
        compare_json = (Join-Path $artifactDir "legacy_runtime_compare.json")
        regression_csv = $regressionCsv
        first_n_l2_error = $observedRegressionFirstN
        status = "passed"
    }
    matlab_device_baseline = $matlabDeviceSummary
    matlab_native_mfem = $matlabNativeSummary
    cpp_legacy_baseline = $cppSummary
    cpp_device_legacy_baseline = $cppDeviceSummary
    status = "passed"
}

$summaryPath = Join-Path $artifactDir "validation_summary.json"
$summary | ConvertTo-Json -Depth 6 | Set-Content -Encoding UTF8 $summaryPath
Write-Host "Legacy validation passed: $summaryPath"
