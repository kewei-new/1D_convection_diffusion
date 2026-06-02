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

function Assert-CsvFieldsClose(
    [object]$Row,
    [hashtable]$Expected,
    [double]$Tolerance,
    [string]$Label
) {
    foreach ($field in $Expected.Keys) {
        $property = $Row.PSObject.Properties[$field]
        if ($null -eq $property) {
            throw "$Label missing field $field"
        }
        $observed = [double]$property.Value
        $expectedValue = [double]$Expected[$field]
        $absDiff = [Math]::Abs($observed - $expectedValue)
        if ($absDiff -gt $Tolerance) {
            throw "$Label $field mismatch: got $observed expected $expectedValue abs_diff=$absDiff tolerance=$Tolerance"
        }
    }
}

function Compare-NumericCsvTable(
    [string]$ExpectedPath,
    [string]$ObservedPath,
    [string]$Label,
    [double]$Tolerance
) {
    if (-not (Test-Path $ExpectedPath)) {
        throw "$Label expected CSV not found: $ExpectedPath"
    }
    if (-not (Test-Path $ObservedPath)) {
        throw "$Label observed CSV not found: $ObservedPath"
    }
    $expectedRows = @(Import-Csv $ExpectedPath)
    $observedRows = @(Import-Csv $ObservedPath)
    if ($expectedRows.Count -ne $observedRows.Count) {
        throw "$Label row-count mismatch: got $($observedRows.Count) expected $($expectedRows.Count)"
    }
    $columns = @($expectedRows[0].PSObject.Properties.Name)
    $maxAbs = 0.0
    $maxRel = 0.0
    for ($i = 0; $i -lt $expectedRows.Count; $i++) {
        foreach ($column in $columns) {
            $observedProperty = $observedRows[$i].PSObject.Properties[$column]
            if ($null -eq $observedProperty) {
                throw "$Label missing observed column $column"
            }
            $expectedValue = [double]$expectedRows[$i].PSObject.Properties[$column].Value
            $observedValue = [double]$observedProperty.Value
            if ([double]::IsNaN($expectedValue) -and [double]::IsNaN($observedValue)) {
                continue
            }
            $absDiff = [Math]::Abs($observedValue - $expectedValue)
            $relDiff = $absDiff / [Math]::Max([Math]::Abs($expectedValue), [double]::Epsilon)
            $maxAbs = [Math]::Max($maxAbs, $absDiff)
            $maxRel = [Math]::Max($maxRel, $relDiff)
            if ($absDiff -gt $Tolerance) {
                throw "$Label $column row $i mismatch: got $observedValue expected $expectedValue abs_diff=$absDiff tolerance=$Tolerance"
            }
        }
    }
    return [ordered]@{
        rows = $expectedRows.Count
        columns = $columns
        max_abs_diff = $maxAbs
        max_rel_diff = $maxRel
        pass = $true
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
    $deviceTerminalCurrent = -13113.0
    $deviceExpectedSummary = @{
        iv_rows = 9.0
        cv_rows = 7.0
        transient_rows = 108.0
        iv_reverse_current_minus1v = -18006.0
        iv_zero_bias_current = 0.019603
        iv_forward_current_1v = 17981.0
        iv_zero_bias_qmag = 2612.5
        cv_zero_bias_cqs = -24.835
        transient_first_finite_time = 0.001875
        transient_first_finite_current = -406000.0
        transient_terminal_time = 0.2
        transient_terminal_current = -13113.0
    }
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
    $deviceExpectedSummary = @{}
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
    $pnOperatorJson = Convert-ToMatlabPath (Join-Path $artifactDir "native_pn_operator_snapshot.json")
    $matlabCommand += " pnos=mfemdd.dd1d_native_pn_operator_snapshot('output_json','$pnOperatorJson');" +
        " assert(pnos.space.total_dofs==480);" +
        " assert(abs(pnos.current.initial_right_contact-937499.5687136018)<1e-6);" +
        " assert(abs(pnos.current.step_right_contact-443761.54836405866)<1e-6);"
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
$matlabPNOperatorSummary = $null
if (-not [string]::IsNullOrWhiteSpace($deviceCompareFunction)) {
    $deviceReport = Get-Content (Join-Path $artifactDir "legacy_device_compare.json") -Raw | ConvertFrom-Json
    $matlabDeviceSummary = [ordered]@{
        compare_json = (Join-Path $artifactDir "legacy_device_compare.json")
        backend = "matlab_mfem"
        mode = "native_recompute_vs_legacy_csv"
        current_status = $deviceReport.current_status
        iv_forward_current_1v = $deviceForwardCurrent
        cv_zero_bias_cqs = $deviceCqs
        status = "passed"
    }

    $pnOperatorSnapshotPath = Join-Path $artifactDir "native_pn_operator_snapshot.json"
    $pnOperatorSnapshot = Get-Content $pnOperatorSnapshotPath -Raw | ConvertFrom-Json
    $expectedPNOperator = [ordered]@{
        total_dofs = 480.0
        n0_norm2 = 3542566.89286449
        rhs0_norm2 = 403973.48589950806
        n_step_norm2 = 3517692.30739153
        initial_right_contact = 937499.5687136018
        step_right_contact = 443761.54836405866
    }
    $observedPNOperator = [ordered]@{
        total_dofs = [double]$pnOperatorSnapshot.space.total_dofs
        n0_norm2 = [double]$pnOperatorSnapshot.vectors.n0.norm2
        rhs0_norm2 = [double]$pnOperatorSnapshot.vectors.rhs0.norm2
        n_step_norm2 = [double]$pnOperatorSnapshot.vectors.n_step.norm2
        initial_right_contact = [double]$pnOperatorSnapshot.current.initial_right_contact
        step_right_contact = [double]$pnOperatorSnapshot.current.step_right_contact
    }
    foreach ($field in $expectedPNOperator.Keys) {
        $absDiff = [Math]::Abs([double]$observedPNOperator[$field] - [double]$expectedPNOperator[$field])
        if ($absDiff -gt 1.0e-6) {
            throw "MATLAB PN operator snapshot $field mismatch: got $($observedPNOperator[$field]) expected $($expectedPNOperator[$field]) abs_diff=$absDiff"
        }
    }
    $matlabPNOperatorSummary = [ordered]@{
        snapshot_json = $pnOperatorSnapshotPath
        backend = "native_matlab_pn_operator_snapshot"
        total_dofs = [int]$observedPNOperator["total_dofs"]
        n0_norm2 = [double]$observedPNOperator["n0_norm2"]
        rhs0_norm2 = [double]$observedPNOperator["rhs0_norm2"]
        n_step_norm2 = [double]$observedPNOperator["n_step_norm2"]
        initial_right_contact = [double]$observedPNOperator["initial_right_contact"]
        step_right_contact = [double]$observedPNOperator["step_right_contact"]
        checked_summary_fields = $expectedPNOperator.Count
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
$cppDeviceBridgeSummary = $null
$cppDeviceNativeSummary = $null
$cppDeviceNativeTablesSummary = $null
$cppDeviceOperatorSummary = $null
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
        Assert-CsvFieldsClose $cppDeviceRows[0] $deviceExpectedSummary 1.0e-12 `
            "C++ device legacy baseline"
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
            transient_terminal_current = [double]$cppDeviceRows[0].transient_terminal_current
            checked_summary_fields = $deviceExpectedSummary.Count
            status = "passed"
        }

        $cppDeviceBridgeRunDir = Join-Path $artifactDir "cpp_device_matlab_bridge"
        New-Item -ItemType Directory -Force -Path $cppDeviceBridgeRunDir | Out-Null
        Push-Location $cppDeviceBridgeRunDir
        try {
            Invoke-Checked { & $deviceExe -n 16 -o 1 -b matlab_mfem } "C++ device MATLAB-native bridge app"
        } finally {
            Pop-Location
        }

        $cppDeviceBridgeCsv = Join-Path $cppDeviceBridgeRunDir "metrics.csv"
        $cppDeviceBridgeRows = Import-Csv $cppDeviceBridgeCsv
        if ($cppDeviceBridgeRows[0].status -ne "native_device_mfem") {
            throw "C++ device MATLAB bridge status mismatch: got $($cppDeviceBridgeRows[0].status)"
        }
        Assert-CsvFieldsClose $cppDeviceBridgeRows[0] $deviceExpectedSummary 1.0e-12 `
            "C++ device MATLAB bridge"
        $bridgeForwardCurrent = [double]$cppDeviceBridgeRows[0].iv_forward_current_1v
        $bridgeCqs = [double]$cppDeviceBridgeRows[0].cv_zero_bias_cqs
        $bridgeTerminalCurrent = [double]$cppDeviceBridgeRows[0].transient_terminal_current
        if ([Math]::Abs($bridgeForwardCurrent - $deviceForwardCurrent) -gt 1.0e-12) {
            throw "C++ device bridge forward-current mismatch: got $bridgeForwardCurrent expected $deviceForwardCurrent"
        }
        if ([Math]::Abs($bridgeCqs - $deviceCqs) -gt 1.0e-12) {
            throw "C++ device bridge Cqs mismatch: got $bridgeCqs expected $deviceCqs"
        }
        if ([Math]::Abs($bridgeTerminalCurrent - $deviceTerminalCurrent) -gt 1.0e-12) {
            throw "C++ device bridge terminal-current mismatch: got $bridgeTerminalCurrent expected $deviceTerminalCurrent"
        }
        $cppDeviceBridgeSummary = [ordered]@{
            app = "dd_device"
            backend = "matlab_mfem"
            metrics_csv = $cppDeviceBridgeCsv
            iv_forward_current_1v = $bridgeForwardCurrent
            cv_zero_bias_cqs = $bridgeCqs
            transient_terminal_current = $bridgeTerminalCurrent
            checked_summary_fields = $deviceExpectedSummary.Count
            status = "passed"
        }

        $cppDeviceNativeRunDir = Join-Path $artifactDir "cpp_device_native_mfem"
        New-Item -ItemType Directory -Force -Path $cppDeviceNativeRunDir | Out-Null
        Push-Location $cppDeviceNativeRunDir
        try {
            Invoke-Checked { & $deviceExe -n 16 -o 1 -b native_mfem } "C++ native device summary app"
        } finally {
            Pop-Location
        }

        $cppDeviceNativeCsv = Join-Path $cppDeviceNativeRunDir "metrics.csv"
        $cppDeviceNativeRows = Import-Csv $cppDeviceNativeCsv
        if ($cppDeviceNativeRows[0].status -ne "native_cpp_device_table") {
            throw "C++ native device summary status mismatch: got $($cppDeviceNativeRows[0].status)"
        }
        Assert-CsvFieldsClose $cppDeviceNativeRows[0] $deviceExpectedSummary 1.0e-12 `
            "C++ native device summary"
        $cppDeviceNativeSummary = [ordered]@{
            app = "dd_device"
            backend = "native_mfem"
            metrics_csv = $cppDeviceNativeCsv
            iv_forward_current_1v = [double]$cppDeviceNativeRows[0].iv_forward_current_1v
            cv_zero_bias_cqs = [double]$cppDeviceNativeRows[0].cv_zero_bias_cqs
            transient_terminal_current = [double]$cppDeviceNativeRows[0].transient_terminal_current
            checked_summary_fields = $deviceExpectedSummary.Count
            status = "passed"
        }

        $cppDeviceNativeTablesRunDir = Join-Path $artifactDir "cpp_device_native_tables"
        New-Item -ItemType Directory -Force -Path $cppDeviceNativeTablesRunDir | Out-Null
        Push-Location $cppDeviceNativeTablesRunDir
        try {
            Invoke-Checked { & $deviceExe -n 16 -o 1 -b native_table } "C++ native PN device full-table app"
        } finally {
            Pop-Location
        }

        $cppDeviceNativeTablesCsv = Join-Path $cppDeviceNativeTablesRunDir "metrics.csv"
        $cppDeviceNativeTablesRows = Import-Csv $cppDeviceNativeTablesCsv
        if ($cppDeviceNativeTablesRows[0].status -ne "native_cpp_device_tables") {
            throw "C++ native PN device full-table status mismatch: got $($cppDeviceNativeTablesRows[0].status)"
        }
        Assert-CsvFieldsClose $cppDeviceNativeTablesRows[0] $deviceExpectedSummary 1.0e-12 `
            "C++ native PN device full-table summary"

        $deviceCompareReport = Get-Content (Join-Path $artifactDir "legacy_device_compare.json") -Raw -Encoding UTF8 | ConvertFrom-Json
        $legacyDeviceRoot = [string]$deviceCompareReport.legacy_source_root
        $ivReport = Compare-NumericCsvTable `
            -ExpectedPath (Join-Path $legacyDeviceRoot "iv_curve.csv") `
            -ObservedPath (Join-Path $cppDeviceNativeTablesRunDir "iv_curve.csv") `
            -Label "C++ native PN IV table" `
            -Tolerance 1.0e-12
        $cvReport = Compare-NumericCsvTable `
            -ExpectedPath (Join-Path $legacyDeviceRoot "cv_curve.csv") `
            -ObservedPath (Join-Path $cppDeviceNativeTablesRunDir "cv_curve.csv") `
            -Label "C++ native PN CV table" `
            -Tolerance 1.0e-12
        $transientReport = Compare-NumericCsvTable `
            -ExpectedPath (Join-Path $legacyDeviceRoot "transient_current.csv") `
            -ObservedPath (Join-Path $cppDeviceNativeTablesRunDir "transient_current.csv") `
            -Label "C++ native PN transient table" `
            -Tolerance 1.0e-12
        $cppDeviceNativeTablesSummary = [ordered]@{
            app = "dd_device"
            backend = "native_table"
            metrics_csv = $cppDeviceNativeTablesCsv
            iv_csv = (Join-Path $cppDeviceNativeTablesRunDir "iv_curve.csv")
            cv_csv = (Join-Path $cppDeviceNativeTablesRunDir "cv_curve.csv")
            transient_csv = (Join-Path $cppDeviceNativeTablesRunDir "transient_current.csv")
            iv_rows = [int]$ivReport.rows
            cv_rows = [int]$cvReport.rows
            transient_rows = [int]$transientReport.rows
            max_abs_diff_vs_legacy = [Math]::Max($ivReport.max_abs_diff, [Math]::Max($cvReport.max_abs_diff, $transientReport.max_abs_diff))
            max_rel_diff_vs_legacy = [Math]::Max($ivReport.max_rel_diff, [Math]::Max($cvReport.max_rel_diff, $transientReport.max_rel_diff))
            abs_tolerance = 1.0e-12
            status = "passed"
        }

        $cppDeviceOperatorRunDir = Join-Path $artifactDir "cpp_device_native_operator_snapshot"
        New-Item -ItemType Directory -Force -Path $cppDeviceOperatorRunDir | Out-Null
        Push-Location $cppDeviceOperatorRunDir
        try {
            Invoke-Checked { & $deviceExe -n 160 -o 2 -b native_operator_snapshot } "C++ native PN operator snapshot app"
        } finally {
            Pop-Location
        }

        $cppDeviceOperatorCsv = Join-Path $cppDeviceOperatorRunDir "metrics.csv"
        $cppDeviceOperatorRows = Import-Csv $cppDeviceOperatorCsv
        if ($cppDeviceOperatorRows[0].status -ne "native_cpp_pn_operator_snapshot") {
            throw "C++ native PN operator snapshot status mismatch: got $($cppDeviceOperatorRows[0].status)"
        }
        $cppPNOperatorExpected = [ordered]@{
            total_dofs = [double]$pnOperatorSnapshot.space.total_dofs
            diffusion_fro_norm = [double]$pnOperatorSnapshot.operators.diffusion.fro_norm
            poisson_rhs_fro_norm = [double]$pnOperatorSnapshot.operators.poisson_rhs.fro_norm
            implicit_fro_norm = [double]$pnOperatorSnapshot.operators.implicit.fro_norm
            n0_norm2 = [double]$pnOperatorSnapshot.vectors.n0.norm2
            phi0_norm2 = [double]$pnOperatorSnapshot.vectors.phi0.norm2
            E0_norm2 = [double]$pnOperatorSnapshot.vectors.E0.norm2
            rhs0_norm2 = [double]$pnOperatorSnapshot.vectors.rhs0.norm2
            n1_norm2 = [double]$pnOperatorSnapshot.vectors.n1.norm2
            rhs1_norm2 = [double]$pnOperatorSnapshot.vectors.rhs1.norm2
            n2_norm2 = [double]$pnOperatorSnapshot.vectors.n2.norm2
            rhs2_norm2 = [double]$pnOperatorSnapshot.vectors.rhs2.norm2
            n3_norm2 = [double]$pnOperatorSnapshot.vectors.n3.norm2
            rhs3_norm2 = [double]$pnOperatorSnapshot.vectors.rhs3.norm2
            n_step_norm2 = [double]$pnOperatorSnapshot.vectors.n_step.norm2
            phi_step_norm2 = [double]$pnOperatorSnapshot.vectors.phi_step.norm2
            E_step_norm2 = [double]$pnOperatorSnapshot.vectors.E_step.norm2
            initial_right_contact = [double]$pnOperatorSnapshot.current.initial_right_contact
            step_right_contact = [double]$pnOperatorSnapshot.current.step_right_contact
        }
        $cppPNOperatorMaxAbs = 0.0
        foreach ($field in $cppPNOperatorExpected.Keys) {
            $property = $cppDeviceOperatorRows[0].PSObject.Properties[$field]
            if ($null -eq $property) {
                throw "C++ native PN operator snapshot missing field $field"
            }
            $observedValue = [double]$property.Value
            $expectedValue = [double]$cppPNOperatorExpected[$field]
            $absDiff = [Math]::Abs($observedValue - $expectedValue)
            $cppPNOperatorMaxAbs = [Math]::Max($cppPNOperatorMaxAbs, $absDiff)
            if ($absDiff -gt 1.0e-5) {
                throw "C++ native PN operator snapshot $field mismatch: got $observedValue expected $expectedValue abs_diff=$absDiff tolerance=1e-5"
            }
        }
        $cppDeviceOperatorSummary = [ordered]@{
            app = "dd_device"
            backend = "native_operator_snapshot"
            metrics_csv = $cppDeviceOperatorCsv
            total_dofs = [int]$cppDeviceOperatorRows[0].total_dofs
            diffusion_fro_norm = [double]$cppDeviceOperatorRows[0].diffusion_fro_norm
            rhs0_norm2 = [double]$cppDeviceOperatorRows[0].rhs0_norm2
            n_step_norm2 = [double]$cppDeviceOperatorRows[0].n_step_norm2
            step_right_contact = [double]$cppDeviceOperatorRows[0].step_right_contact
            checked_summary_fields = $cppPNOperatorExpected.Count
            max_abs_diff_vs_matlab_operator_snapshot = $cppPNOperatorMaxAbs
            abs_tolerance = 1.0e-5
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
    matlab_pn_operator_snapshot = $matlabPNOperatorSummary
    matlab_native_mfem = $matlabNativeSummary
    cpp_legacy_baseline = $cppSummary
    cpp_native_mfem = $cppNativeSummary
    cpp_device_legacy_baseline = $cppDeviceSummary
    cpp_device_matlab_bridge = $cppDeviceBridgeSummary
    cpp_device_native_mfem = $cppDeviceNativeSummary
    cpp_device_native_tables = $cppDeviceNativeTablesSummary
    cpp_device_native_operator_snapshot = $cppDeviceOperatorSummary
    status = "passed"
}

$summaryPath = Join-Path $artifactDir "validation_summary.json"
$summary | ConvertTo-Json -Depth 6 | Set-Content -Encoding UTF8 $summaryPath
& (Join-Path $PSScriptRoot "write_alignment_audit.ps1") -ValidationSummary $summaryPath
Write-Host "Legacy validation passed: $summaryPath"
