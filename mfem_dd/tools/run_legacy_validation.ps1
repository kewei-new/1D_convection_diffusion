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
} elseif ($repoName -like "2D_convection_diffusion") {
    $compareFunction = "mfemdd.compare_legacy_dd2d"
    $refineSteps = 4
    $order = 1
    $regressionFirstN = 4.507062201109525e-01
    $cppApp = "dd2d_mms"
    $cppArgs = @("-n", "6", "-o", "1", "-b", "legacy_baseline")
    $cppFirstN = 4.507062e-01
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

Invoke-Checked { matlab -batch $matlabCommand } "MATLAB legacy validation"

$regressionCsv = Join-Path (Join-Path $artifactDir "matlab_regression") "metrics.csv"
$regressionRows = Import-Csv $regressionCsv
$observedRegressionFirstN = [double]$regressionRows[0].n_l2_error
if ([Math]::Abs($observedRegressionFirstN - $regressionFirstN) -gt 1.0e-12) {
    throw "MATLAB regression first n_L2 mismatch: got $observedRegressionFirstN expected $regressionFirstN"
}

$cppSummary = $null
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
    cpp_legacy_baseline = $cppSummary
    status = "passed"
}

$summaryPath = Join-Path $artifactDir "validation_summary.json"
$summary | ConvertTo-Json -Depth 6 | Set-Content -Encoding UTF8 $summaryPath
Write-Host "Legacy validation passed: $summaryPath"
