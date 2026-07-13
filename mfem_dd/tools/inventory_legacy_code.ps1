param(
    [string[]]$SourceRoots = @(),
    [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
if ($SourceRoots.Count -eq 0) {
    $SourceRoots = @(
        (Get-ChildItem -Path $repoRoot -Directory | Where-Object {
            $_.Name -ne ".git"
        } | ForEach-Object FullName)
    )
}
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $repoRoot "mfem_dd\artifacts\code_inventory"
}
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

$extensions = @(".m", ".cpp", ".c", ".h", ".hpp", ".py", ".ps1")
$excludedPathParts = @("\\.git\\", "\\build\\", "\\artifacts\\", "\\__pycache__\\")

function Get-RelativePath([string]$Path) {
    return ($Path.Substring($repoRoot.Length).TrimStart([char[]]@('\', '/')) -replace '\\', '/')
}

function Get-Classification([System.IO.FileInfo]$File) {
    $relative = Get-RelativePath $File.FullName
    $name = $File.Name.ToLowerInvariant()
    $sample = ""
    try {
        $sample = (Get-Content -LiteralPath $File.FullName -TotalCount 240 -ErrorAction Stop) -join "`n"
    } catch {
        $sample = ""
    }
    $haystack = ("$relative`n$sample").ToLowerInvariant()
    $tags = [System.Collections.Generic.List[string]]::new()

    if ($haystack -match "\b(ipdg|sipg|interior.penalty)\b") { $tags.Add("space:IPDG/SIPG") }
    if ($haystack -match "\bldg\b|local.discontinuous.galerkin") { $tags.Add("space:LDG") }
    if ($haystack -match "\bdg\b|discontinuous.galerkin") { $tags.Add("space:DG-general") }
    if ($haystack -match "imex") { $tags.Add("time:IMEX-RK") }
    if ($haystack -match "tvd[_ -]?rk|ssp[_ -]?rk") { $tags.Add("time:TVD-RK") }
    if ($haystack -match "\brk[234]\b|runge.?kutta") { $tags.Add("time:RK") }
    if ($haystack -match "drift.?diffusion|\bdd\b|bipolar|pn.?junction|carrier") { $tags.Add("equation:DD") }
    if ($haystack -match "elliptic|poisson") { $tags.Add("equation:elliptic") }
    if ($haystack -match "parabolic|heat|convection.?diffusion") { $tags.Add("equation:parabolic/transport") }
    if ($haystack -match "periodic") { $tags.Add("bc:periodic") }
    if ($haystack -match "dirichlet") { $tags.Add("bc:Dirichlet") }
    if ($haystack -match "neumann") { $tags.Add("bc:Neumann") }
    if ($haystack -match "robin") { $tags.Add("bc:Robin") }
    if ($haystack -match "limiter|tvb|minmod") { $tags.Add("feature:limiter") }
    if ($haystack -match "assemble|matric") { $tags.Add("role:assembly") }
    if ($haystack -match "solve|gmres|cg\b|lu\\|\\|backslash") { $tags.Add("role:linear-or-pde-solve") }
    if ($haystack -match "projection|l2_projection") { $tags.Add("role:projection") }
    if ($haystack -match "quadrature|gauss|integral|weight") { $tags.Add("infra:quadrature") }
    if ($haystack -match "mesh|grid|rectangle|triang") { $tags.Add("infra:mesh") }
    if ($haystack -match "plot|visual|save_solution|postprocess|evaluate_.*error") { $tags.Add("infra:postprocess") }
    if ($haystack -match "main\.m$|test|smoke|example|run_") { $tags.Add("role:driver-or-test") }
    if ($haystack -match "config|case|startup|runtime|interface") { $tags.Add("infra:framework") }

    $primary = "unclassified"
    foreach ($prefix in @("equation:DD", "equation:elliptic", "equation:parabolic/transport",
                           "space:IPDG/SIPG", "space:LDG", "time:IMEX-RK", "time:TVD-RK",
                           "feature:limiter", "role:assembly", "infra:mesh", "infra:quadrature",
                           "infra:postprocess", "infra:framework")) {
        if ($tags -contains $prefix) { $primary = $prefix; break }
    }

    $normalizedName = $name -replace "\s*\(\d+\)", "" -replace "copy[_ -]?of[_ -]?", ""
    return [pscustomobject]@{
        relative_path = $relative
        source_root = ($relative -split "/")[0]
        extension = $File.Extension.ToLowerInvariant()
        bytes = $File.Length
        file_name = $File.Name
        normalized_name = $normalizedName
        primary_class = $primary
        tags = ($tags | Select-Object -Unique) -join ";"
    }
}

$files = foreach ($root in $SourceRoots) {
    if (-not (Test-Path -LiteralPath $root)) { continue }
    Get-ChildItem -LiteralPath $root -Recurse -File | Where-Object {
        $path = $_.FullName
        $extensions -contains $_.Extension.ToLowerInvariant() -and
        -not ($excludedPathParts | Where-Object { $path -like "*$_*" })
    }
}

$inventory = foreach ($file in $files) {
    $item = Get-Classification $file
    $hash = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    [pscustomobject]@{
        relative_path = $item.relative_path
        source_root = $item.source_root
        extension = $item.extension
        bytes = $item.bytes
        sha256 = $hash
        file_name = $item.file_name
        normalized_name = $item.normalized_name
        primary_class = $item.primary_class
        tags = $item.tags
    }
}

$inventory | Sort-Object relative_path | Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $OutputDir "source_inventory.csv")

$duplicates = foreach ($group in ($inventory | Group-Object sha256 | Where-Object Count -gt 1)) {
    $canonical = ($group.Group | Sort-Object relative_path | Select-Object -First 1).relative_path
    foreach ($item in $group.Group | Sort-Object relative_path) {
        [pscustomobject]@{ sha256 = $group.Name; copies = $group.Count; canonical_path = $canonical; relative_path = $item.relative_path; primary_class = $item.primary_class }
    }
}
$duplicates | Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $OutputDir "exact_duplicate_groups.csv")

$overlapCandidates = foreach ($group in ($inventory | Group-Object normalized_name | Where-Object Count -gt 1)) {
    $distinctHashes = @($group.Group.sha256 | Select-Object -Unique)
    if ($distinctHashes.Count -gt 1) {
        foreach ($item in $group.Group | Sort-Object relative_path) {
            [pscustomobject]@{ normalized_name = $group.Name; candidates = $group.Count; distinct_hashes = $distinctHashes.Count; relative_path = $item.relative_path; sha256 = $item.sha256; primary_class = $item.primary_class; tags = $item.tags }
        }
    }
}
$overlapCandidates | Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $OutputDir "overlap_candidates.csv")

$manualReview = @()
$manualReview += $inventory | Where-Object { $_.primary_class -eq "unclassified" } |
    ForEach-Object {
        [pscustomobject]@{ review_kind = "missing-classification"; relative_path = $_.relative_path; normalized_name = $_.normalized_name; primary_class = $_.primary_class; tags = $_.tags; decision = "pending"; question = "Assign numerical role and MFEM destination." }
    }
$manualReview += $overlapCandidates | ForEach-Object {
    [pscustomobject]@{ review_kind = "nonidentical-same-name"; relative_path = $_.relative_path; normalized_name = $_.normalized_name; primary_class = $_.primary_class; tags = $_.tags; decision = "pending"; question = "Compare against the selected canonical variant; record boundary, equation, limiter, or time-step difference." }
}
$manualReview | Sort-Object review_kind, normalized_name, relative_path |
    Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $OutputDir "manual_review_queue.csv")

$summary = [ordered]@{
    generated_at = (Get-Date).ToString("o")
    repository = Split-Path $repoRoot -Leaf
    source_roots = @($SourceRoots | ForEach-Object { Get-RelativePath $_ })
    source_files = @($inventory).Count
    exact_duplicate_files = @($duplicates).Count
    exact_duplicate_groups = @($duplicates | Group-Object sha256).Count
    overlap_candidate_files = @($overlapCandidates).Count
    manual_review_items = @($manualReview).Count
    by_primary_class = @($inventory | Group-Object primary_class | Sort-Object Name | ForEach-Object { [ordered]@{ class = $_.Name; files = $_.Count } })
    by_source_root = @($inventory | Group-Object source_root | Sort-Object Name | ForEach-Object { [ordered]@{ root = $_.Name; files = $_.Count } })
}
$summary | ConvertTo-Json -Depth 5 | Set-Content -Encoding utf8 (Join-Path $OutputDir "inventory_summary.json")
Write-Output ("Inventory written to {0}: {1} files, {2} exact duplicate groups, {3} overlap candidate files." -f $OutputDir, $summary.source_files, $summary.exact_duplicate_groups, $summary.overlap_candidate_files)
