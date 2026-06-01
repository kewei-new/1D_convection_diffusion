param(
    [string]$ValidationSummary = "",
    [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$repoName = Split-Path $repoRoot -Leaf
$mfemRoot = Join-Path $repoRoot "mfem_dd"
if ([string]::IsNullOrWhiteSpace($ValidationSummary)) {
    $ValidationSummary = Join-Path $mfemRoot "artifacts\legacy_validation\validation_summary.json"
}
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $mfemRoot "artifacts\alignment_audit"
}
if (-not (Test-Path $ValidationSummary)) {
    throw "Validation summary not found: $ValidationSummary"
}

function Read-JsonIfPresent([string]$PathValue) {
    if ([string]::IsNullOrWhiteSpace($PathValue) -or -not (Test-Path $PathValue)) {
        return $null
    }
    return Get-Content $PathValue -Raw | ConvertFrom-Json
}

function Test-Passed($Node) {
    return ($null -ne $Node -and $Node.status -eq "passed")
}

function New-Check([string]$CaseName, [string]$Scope, [string]$Implementation,
                   [string]$Evidence, [bool]$Passed, [string]$Notes) {
    return [ordered]@{
        case = $CaseName
        scope = $Scope
        implementation = $Implementation
        evidence = $Evidence
        result = $(if ($Passed) { "passed" } else { "missing_or_failed" })
        notes = $Notes
    }
}

function New-Gap([string]$CaseName, [string]$Scope, [string]$Reason) {
    return [ordered]@{
        case = $CaseName
        scope = $Scope
        status = "open"
        reason = $Reason
    }
}

$summary = Get-Content $ValidationSummary -Raw | ConvertFrom-Json
$checks = @()
$gaps = @()

if ($repoName -like "1D_convection_diffusion") {
    $nativeMms = Read-JsonIfPresent $summary.matlab_native_mfem.compare_json
    $deviceCompare = Read-JsonIfPresent $summary.matlab_device_baseline.compare_json

    $checks += New-Check "dd1d_smooth_mms" "MATLAB legacy runtime full table" `
        "legacy_runtime" $summary.matlab_legacy_runtime.compare_json `
        (Test-Passed $summary.matlab_legacy_runtime) `
        "Old MATLAB solver is rerun through the unified API and checked against the stored table."
    $checks += New-Check "dd1d_smooth_mms" "MATLAB native full table" `
        "native_matlab_mfem" $summary.matlab_native_mfem.compare_json `
        ((Test-Passed $summary.matlab_native_mfem) -and $nativeMms.matches_full_table) `
        "Self-contained MATLAB-native IPDG/LDG/IMEX table matches the legacy table within tolerances."
    $checks += New-Check "dd1d_smooth_mms" "C++ baseline rows" `
        "legacy_baseline_adapter" "" `
        ((Test-Passed $summary.cpp_legacy_baseline) -and $summary.cpp_legacy_baseline.row_count -eq 5) `
        "C++ adapter emits all embedded legacy MMS rows."
    $checks += New-Check "dd1d_smooth_mms" "C++ native rows" `
        "native_cpp_mfem" "" `
        ((Test-Passed $summary.cpp_native_mfem) -and $summary.cpp_native_mfem.row_count -eq 5) `
        "C++ modal DG/IPDG/LDG/IMEX backend is compared against the MATLAB-native table."
    $checks += New-Check "dd_pn_device" "MATLAB physical outputs" `
        "native_matlab_mfem" $summary.matlab_device_baseline.compare_json `
        ((Test-Passed $summary.matlab_device_baseline) -and $deviceCompare.matches_baseline -and $deviceCompare.matches_full_outputs) `
        "The PN device solve is recomputed from the native MATLAB source mirror and compared to the legacy IV/CV/transient CSV outputs."
    $checks += New-Check "dd_pn_device" "C++ physical outputs" `
        "matlab_native_bridge" $summary.cpp_device_matlab_bridge.metrics_csv `
        (Test-Passed $summary.cpp_device_matlab_bridge) `
        "C++ device entry point invokes the MATLAB-native PN solve and emits selected matched metrics."

    $gaps += New-Gap "dd_pn_device" "C++ native physical solve" `
        "The C++ device path has a MATLAB-native bridge, but no native PN device solver is implemented."
    $gaps += New-Gap "legacy cleanup" "old directory deletion" `
        "No legacy directory is eligible for deletion until each target has legacy, MATLAB-native, and C++ evidence."
} elseif ($repoName -like "2D_convection_diffusion") {
    $nativeMms = Read-JsonIfPresent $summary.matlab_native_mfem.compare_json
    $deviceCompare = Read-JsonIfPresent $summary.matlab_device_baseline.compare_json

    $checks += New-Check "dd2d_smooth_mms" "MATLAB legacy runtime full table" `
        "legacy_runtime" $summary.matlab_legacy_runtime.compare_json `
        (Test-Passed $summary.matlab_legacy_runtime) `
        "Old MATLAB solver is rerun through the unified API and checked against the stored MMS table."
    $checks += New-Check "dd2d_smooth_mms" "MATLAB native full table" `
        "native_matlab_mfem" $summary.matlab_native_mfem.compare_json `
        ((Test-Passed $summary.matlab_native_mfem) -and $nativeMms.matches_full_table) `
        "Self-contained MATLAB-native source mirror matches the legacy MMS table within tolerances."
    $checks += New-Check "dd2d_smooth_mms" "C++ baseline rows" `
        "legacy_baseline_adapter" "" `
        ((Test-Passed $summary.cpp_legacy_baseline) -and $summary.cpp_legacy_baseline.row_count -eq 4) `
        "C++ adapter emits all embedded legacy MMS rows."
    $checks += New-Check "dd2d_smooth_mms" "C++ MATLAB-native bridge" `
        "matlab_native_bridge" "" `
        (Test-Passed $summary.cpp_matlab_native_bridge) `
        "C++ entry point invokes the MATLAB-native solver and checks the first full-precision native row."
    $checks += New-Check "dd_pn_device" "MATLAB native physical outputs" `
        "native_matlab_mfem" $summary.matlab_device_baseline.compare_json `
        ((Test-Passed $summary.matlab_device_baseline) -and $deviceCompare.matches_baseline -and $deviceCompare.matches_full_outputs) `
        "The selected device artifact is recomputed from the native MATLAB mirror and full coefficient arrays match."
    $checks += New-Check "dd_pn_device" "C++ physical outputs" `
        "legacy_baseline_adapter" "" `
        (Test-Passed $summary.cpp_device_legacy_baseline) `
        "C++ device app emits selected 2D device baseline metrics only."

    $gaps += New-Gap "dd2d_smooth_mms" "C++ native MMS solve" `
        "The C++ MMS path has a baseline adapter and MATLAB-native bridge, but no native C++ IPDG/LDG/IMEX implementation."
    $gaps += New-Gap "dd_pn_device" "C++ native physical solve" `
        "The C++ device path still emits selected baseline metrics; no native 2D device solver is implemented."
    $gaps += New-Gap "legacy cleanup" "old directory deletion" `
        "No legacy directory is eligible for deletion until each target has legacy, MATLAB-native, and C++ evidence."
} else {
    throw "Unsupported repository for alignment audit: $repoName"
}

$failed = @($checks | Where-Object { $_.result -ne "passed" })
$audit = [ordered]@{
    repo = $repoName
    generated_at = (Get-Date).ToString("s")
    validation_summary = (Resolve-Path $ValidationSummary).Path
    checks = $checks
    failed_checks = $failed
    remaining_gaps = $gaps
    complete_for_legacy_deletion = ($failed.Count -eq 0 -and $gaps.Count -eq 0)
    status = $(if ($failed.Count -eq 0) { "audited_with_open_gaps" } else { "failed" })
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$outPath = Join-Path $OutputDir "alignment_audit.json"
$audit | ConvertTo-Json -Depth 8 | Set-Content -Encoding UTF8 $outPath
Write-Host "Alignment audit written: $outPath"
if ($failed.Count -gt 0) {
    throw "Alignment audit has $($failed.Count) failed check(s)."
}
