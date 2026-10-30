$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$taskPython = if (Test-Path -LiteralPath 'C:\Python314\python.exe') { 'C:\Python314\python.exe' } else { (Get-Command python -ErrorAction Stop).Source }
$logRoot = Join-Path $projectRoot '.code_publish'
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null
Push-Location -LiteralPath $projectRoot
try {
    & $taskPython -u tools/sync_annotation_github.py *> (Join-Path $logRoot 'last-sync.log')
    if ($LASTEXITCODE -ne 0) { throw 'GitHub code sync failed; see .code_publish/last-sync.log' }
} finally { Pop-Location }
