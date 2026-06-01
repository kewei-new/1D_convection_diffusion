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
    $cppCompareFields = @(
        "n_l2_error", "n_linf_error",
        "phi_l2_error", "phi_linf_error",
        "E_l2_error", "E_linf_error")
    $cppExpectedRows = @(
        [ordered]@{ elements = 20; order = 3; n_l2_error = 4.507504e-06; n_linf_error = 3.897874e-06; phi_l2_error = 5.547920e-06; phi_linf_error = 9.655373e-06; E_l2_error = 5.658847e-06; E_linf_error = 9.737335e-06 },
        [ordered]@{ elements = 40; order = 3; n_l2_error = 4.071754e-07; n_linf_error = 3.531198e-07; phi_l2_error = 3.751691e-07; phi_linf_error = 6.082643e-07; E_l2_error = 3.928129e-07; E_linf_error = 6.268181e-07 },
        [ordered]@{ elements = 80; order = 3; n_l2_error = 4.561335e-08; n_linf_error = 3.706344e-08; phi_l2_error = 2.966955e-08; phi_linf_error = 3.820909e-08; E_l2_error = 3.311548e-08; E_linf_error = 4.220075e-08 },
        [ordered]@{ elements = 160; order = 3; n_l2_error = 5.531522e-09; n_linf_error = 4.186841e-09; phi_l2_error = 2.940053e-09; phi_linf_error = 2.935649e-09; E_l2_error = 3.472993e-09; E_linf_error = 3.100092e-09 },
        [ordered]@{ elements = 320; order = 3; n_l2_error = 6.906262e-10; n_linf_error = 4.991136e-10; phi_l2_error = 3.419493e-10; phi_linf_error = 3.259086e-10; E_l2_error = 4.140009e-10; E_linf_error = 3.365440e-10 }
    )
    $deviceCompareFunction = "mfemdd.compare_legacy_pn1d"
    $deviceForwardCurrent = 17981.0
    $deviceCqs = -24.835
    $nativeCompareFunction = "mfemdd.compare_legacy_dd1d"
    $cppNativeBackend = "native_mfem"
    $cppNativeAbsTol = 2.0e-12
    $cppNativeRelTol = 1.0e-5
    $cppNativeFieldMap = @(
        [ordered]@{ csv = "n_l2_error"; matlab = "n_L2" },
        [ordered]@{ csv = "n_linf_error"; matlab = "n_Linf" },
        [ordered]@{ csv = "phi_l2_error"; matlab = "phi_L2" },
        [ordered]@{ csv = "phi_linf_error"; matlab = "phi_Linf" },
        [ordered]@{ csv = "E_l2_error"; matlab = "E_L2" },
        [ordered]@{ csv = "E_linf_error"; matlab = "E_Linf" }
    )
} elseif ($repoName -like "2D_convection_diffusion") {
    $compareFunction = "mfemdd.compare_legacy_dd2d"
    $refineSteps = 4
    $order = 1
    $regressionFirstN = 4.507062201109525e-01
    $cppApp = "dd2d_mms"
    $cppCompareFields = @("n_l2_error", "p_l2_error", "phi_l2_error", "Ex_l2_error", "Ey_l2_error")
    $cppExpectedRows = @(
        [ordered]@{ elements = 6; order = 1; n_l2_error = 4.507062e-01; p_l2_error = 3.550932e-01; phi_l2_error = 9.804350e-02; Ex_l2_error = 6.716957e-02; Ey_l2_error = 6.716957e-02 },
        [ordered]@{ elements = 12; order = 1; n_l2_error = 2.629943e-02; p_l2_error = 2.298078e-02; phi_l2_error = 3.181550e-02; Ex_l2_error = 2.051898e-02; Ey_l2_error = 2.051898e-02 },
        [ordered]@{ elements = 24; order = 1; n_l2_error = 1.211073e-02; p_l2_error = 9.843311e-03; phi_l2_error = 8.128524e-03; Ex_l2_error = 5.428426e-03; Ey_l2_error = 5.428426e-03 },
        [ordered]@{ elements = 48; order = 1; n_l2_error = 1.313986e-03; p_l2_error = 1.243116e-03; phi_l2_error = 2.030280e-03; Ex_l2_error = 1.379696e-03; Ey_l2_error = 1.379696e-03 }
    )
    $deviceCompareFunction = ""
    $deviceForwardCurrent = $null
    $deviceCqs = $null
    $nativeCompareFunction = ""
    $cppNativeBackend = ""
    $cppNativeAbsTol = $null
    $cppNativeRelTol = $null
    $cppNativeFieldMap = @()
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
$cppNativeSummary = $null
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
    $cppValidatedRows = @()
    foreach ($expected in $cppExpectedRows) {
        $rowRunDir = Join-Path $cppRunDir ("n_{0}" -f $expected.elements)
        New-Item -ItemType Directory -Force -Path $rowRunDir | Out-Null
        $rowArgs = @("-n", [string]$expected.elements, "-o", [string]$expected.order, "-b", "legacy_baseline")
        Push-Location $rowRunDir
        try {
            Invoke-Checked { & $exe @rowArgs } "C++ legacy baseline app"
        } finally {
            Pop-Location
        }

        $cppCsv = Join-Path $rowRunDir "metrics.csv"
        $cppRows = @(Import-Csv $cppCsv)
        if ($cppRows.Count -ne 1) {
            throw "C++ baseline emitted $($cppRows.Count) rows for elements=$($expected.elements); expected 1"
        }

        $observed = $cppRows[0]
        $summaryRow = [ordered]@{
            elements = [int]$expected.elements
            order = [int]$expected.order
            metrics_csv = $cppCsv
        }
        foreach ($field in $cppCompareFields) {
            $observedValue = [double]$observed.PSObject.Properties[$field].Value
            $expectedValue = [double]$expected[$field]
            if ([Math]::Abs($observedValue - $expectedValue) -gt 1.0e-12) {
                throw "C++ baseline $field mismatch for elements=$($expected.elements): got $observedValue expected $expectedValue"
            }
            $summaryRow[$field] = $observedValue
        }
        $summaryRow["status"] = "passed"
        $cppValidatedRows += $summaryRow
    }
    $cppSummary = [ordered]@{
        app = $cppApp
        row_count = $cppValidatedRows.Count
        rows = $cppValidatedRows
        first_n_l2_error = [double]$cppValidatedRows[0]["n_l2_error"]
        status = "passed"
    }

    if (-not [string]::IsNullOrWhiteSpace($cppNativeBackend)) {
        $nativeReport = Get-Content (Join-Path $artifactDir "native_matlab_compare.json") -Raw | ConvertFrom-Json
        $columns = @($nativeReport.current_table.columns)
        $nativeData = @($nativeReport.current_table.data)
        if ($nativeData.Count -lt $cppExpectedRows.Count) {
            throw "MATLAB native table has $($nativeData.Count) rows; expected at least $($cppExpectedRows.Count)"
        }

        $cppNativeRunDir = Join-Path $artifactDir "cpp_native_mfem"
        New-Item -ItemType Directory -Force -Path $cppNativeRunDir | Out-Null
        $cppNativeRows = @()
        $cppNativeMaxAbs = 0.0
        $cppNativeMaxRel = 0.0

        for ($rowIndex = 0; $rowIndex -lt $cppExpectedRows.Count; $rowIndex++) {
            $expectedRow = $cppExpectedRows[$rowIndex]
            $matlabRow = @($nativeData[$rowIndex])
            $rowRunDir = Join-Path $cppNativeRunDir ("n_{0}" -f $expectedRow.elements)
            New-Item -ItemType Directory -Force -Path $rowRunDir | Out-Null
            $rowArgs = @("-n", [string]$expectedRow.elements, "-o", [string]$expectedRow.order, "-b", $cppNativeBackend)
            Push-Location $rowRunDir
            try {
                Invoke-Checked { & $exe @rowArgs } "C++ native MFEM app"
            } finally {
                Pop-Location
            }

            $cppNativeCsv = Join-Path $rowRunDir "metrics.csv"
            $cppNativeCsvRows = @(Import-Csv $cppNativeCsv)
            if ($cppNativeCsvRows.Count -ne 1) {
                throw "C++ native backend emitted $($cppNativeCsvRows.Count) rows for elements=$($expectedRow.elements); expected 1"
            }

            $observed = $cppNativeCsvRows[0]
            $summaryRow = [ordered]@{
                elements = [int]$expectedRow.elements
                order = [int]$expectedRow.order
                metrics_csv = $cppNativeCsv
            }

            foreach ($field in $cppNativeFieldMap) {
                $csvField = [string]$field["csv"]
                $matlabField = [string]$field["matlab"]
                $columnIndex = [Array]::IndexOf($columns, $matlabField)
                if ($columnIndex -lt 0) {
                    throw "MATLAB native table missing column $matlabField"
                }
                $observedValue = [double]$observed.PSObject.Properties[$csvField].Value
                $expectedValue = [double]$matlabRow[$columnIndex]
                $absDiff = [Math]::Abs($observedValue - $expectedValue)
                $relDiff = $absDiff / [Math]::Max([Math]::Abs($expectedValue), [double]::Epsilon)
                $cppNativeMaxAbs = [Math]::Max($cppNativeMaxAbs, $absDiff)
                $cppNativeMaxRel = [Math]::Max($cppNativeMaxRel, $relDiff)
                $allowed = [Math]::Max($cppNativeAbsTol, $cppNativeRelTol * [Math]::Abs($expectedValue))
                if ($absDiff -gt $allowed) {
                    throw "C++ native $csvField mismatch for elements=$($expectedRow.elements): got $observedValue expected $expectedValue abs_diff=$absDiff allowed=$allowed"
                }
                $summaryRow[$csvField] = $observedValue
                $summaryRow["${csvField}_abs_diff"] = $absDiff
            }
            $summaryRow["status"] = "passed"
            $cppNativeRows += $summaryRow
        }

        $cppNativeSummary = [ordered]@{
            app = $cppApp
            backend = $cppNativeBackend
            row_count = $cppNativeRows.Count
            rows = $cppNativeRows
            max_abs_diff_vs_matlab_native = $cppNativeMaxAbs
            max_rel_diff_vs_matlab_native = $cppNativeMaxRel
            abs_tolerance = $cppNativeAbsTol
            rel_tolerance = $cppNativeRelTol
            status = "passed"
        }
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
    cpp_native_mfem = $cppNativeSummary
    cpp_device_legacy_baseline = $cppDeviceSummary
    status = "passed"
}

$summaryPath = Join-Path $artifactDir "validation_summary.json"
$summary | ConvertTo-Json -Depth 6 | Set-Content -Encoding UTF8 $summaryPath
& (Join-Path $PSScriptRoot "write_alignment_audit.ps1") -ValidationSummary $summaryPath
Write-Host "Legacy validation passed: $summaryPath"
