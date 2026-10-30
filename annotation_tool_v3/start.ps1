param([int]$Port = 8767, [switch]$Foreground, [switch]$NoOpen)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$serviceUrl = 'http://127.0.0.1:' + $Port

# Reuse the current server; launching this twice must not replace an active job.
try {
    $existing = Invoke-RestMethod -Uri ($serviceUrl + '/api/health') -TimeoutSec 2
} catch { $existing = $null }
if ($existing -and $existing.service -eq 'atomnmr' -and $existing.status -eq 'ok') {
    if (-not $NoOpen) { Start-Process $serviceUrl -WindowStyle Hidden }
    Write-Output ('AtomNMR already running: ' + $serviceUrl)
    return
}
# Older versions have no health endpoint. Recognize them without stopping them.
try {
    $legacy = Invoke-RestMethod -Uri ($serviceUrl + '/api/benchmark/runs') -TimeoutSec 2
} catch { $legacy = $null }
if ($legacy -and $legacy.PSObject.Properties.Name -contains 'runs') {
    if (-not $NoOpen) { Start-Process $serviceUrl -WindowStyle Hidden }
    Write-Output ('AtomNMR already running: ' + $serviceUrl)
    return
}

$pythonPath = if (Test-Path -LiteralPath 'C:\Python314\python.exe') { 'C:\Python314\python.exe' } else { (Get-Command python -ErrorAction Stop).Source }
if ($Foreground) {
    Push-Location -LiteralPath $projectRoot
    try {
        if ($NoOpen) { & $pythonPath -u -m annotation_tool_v3.server --port $Port }
        else { & $pythonPath -u -m annotation_tool_v3.server --port $Port --open }
    } finally { Pop-Location }
    return
}

$logDir = Join-Path $projectRoot 'data\assignment_database_v1\service_logs'
New-Item -ItemType Directory -Path $logDir -Force | Out-Null
$stamp = (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [guid]::NewGuid().ToString('N').Substring(0,8)
$outLog = Join-Path $logDir ($stamp + '.out.log')
$errLog = Join-Path $logDir ($stamp + '.err.log')
$process = Start-Process -FilePath $pythonPath -ArgumentList @('-u','-m','annotation_tool_v3.server','--port',"$Port") -WorkingDirectory $projectRoot -WindowStyle Hidden -RedirectStandardOutput $outLog -RedirectStandardError $errLog -PassThru
for ($attempt = 0; $attempt -lt 30; $attempt++) {
    try { $health = Invoke-RestMethod -Uri ($serviceUrl + '/api/health') -TimeoutSec 1 } catch { $health = $null }
    if ($health -and $health.service -eq 'atomnmr' -and $health.status -eq 'ok') {
        if (-not $NoOpen) { Start-Process $serviceUrl -WindowStyle Hidden }
        Write-Output ('AtomNMR started: ' + $serviceUrl + ' (PID ' + $health.pid + ')')
        return
    }
    $process.Refresh()
    if ($process.HasExited) { throw ('AtomNMR failed to start. Log: ' + $errLog) }
    Start-Sleep -Milliseconds 300
}
throw ('AtomNMR startup not confirmed. Check log before starting again: ' + $errLog)
