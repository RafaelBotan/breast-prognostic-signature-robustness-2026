param(
    [Parameter(Mandatory = $true)]
    [string]$CohortDir,
    [string]$OutputDir = (Join-Path $PSScriptRoot 'outputs')
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $CohortDir -PathType Container)) {
    throw "Processed public-cohort directory not found: $CohortDir"
}

$env:COHORT_DIR = [IO.Path]::GetFullPath($CohortDir)
$env:STRESS_PROJECT_DIR = $PSScriptRoot
$env:SUBMITTED_AGGREGATE_DIR = Join-Path $PSScriptRoot 'data\reference_aggregate'
$env:REVISION_OUTPUT_DIR = [IO.Path]::GetFullPath($OutputDir)
$env:REVISION_CACHE_DIR = Join-Path $PSScriptRoot 'cache'
$env:FIGURE_OUTPUT_DIR = Join-Path $env:REVISION_OUTPUT_DIR 'figures'
$env:B_MATCH = '1000'

New-Item -ItemType Directory -Force -Path $env:REVISION_OUTPUT_DIR | Out-Null
New-Item -ItemType Directory -Force -Path $env:REVISION_CACHE_DIR | Out-Null
New-Item -ItemType Directory -Force -Path $env:FIGURE_OUTPUT_DIR | Out-Null

foreach ($script in @(
    'code\revision\reviewer1_revision_analysis.R',
    'code\revision\postprocess_matched_outputs.R',
    'code\revision\stress_e2_platform_revision.R',
    'code\revision\make_figure1_revision.R',
    'code\revision\make_figure2_revision.R'
)) {
    & Rscript (Join-Path $PSScriptRoot $script)
    if ($LASTEXITCODE -ne 0) { throw "Rscript failed: $script" }
}

Write-Output "Revision outputs written to $env:REVISION_OUTPUT_DIR"

