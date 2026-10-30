$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$taskPython = if (Test-Path -LiteralPath 'C:\Python314\python.exe') { 'C:\Python314\python.exe' } else { (Get-Command python -ErrorAction Stop).Source }
$logRoot = Join-Path $projectRoot '.code_publish'
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null
Push-Location -LiteralPath $projectRoot
try {
    # Task Scheduler does not inherit the Codex/terminal proxy environment.
    # Read the current user's system proxy each time, so FlClash port changes
    # do not require editing or storing any credential in this script.
    $proxySettings = Get-ItemProperty -LiteralPath 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Internet Settings' -ErrorAction SilentlyContinue
    if ($proxySettings.ProxyEnable -eq 1 -and $proxySettings.ProxyServer) {
        $proxyMap = @{}
        foreach ($entry in ([string]$proxySettings.ProxyServer).Split(';')) {
            if ($entry -match '^([^=]+)=(.+)$') { $proxyMap[$Matches[1]] = $Matches[2] }
            elseif ($entry.Trim()) { $proxyMap['default'] = $entry.Trim() }
        }
        foreach ($protocol in @('http','https')) {
            $selectedProxy = if ($proxyMap.ContainsKey($protocol)) { $proxyMap[$protocol] } else { $proxyMap['default'] }
            if ($selectedProxy) {
                if ($selectedProxy -notmatch '^[a-z]+://') { $selectedProxy = 'http://' + $selectedProxy }
                [Environment]::SetEnvironmentVariable(($protocol.ToUpperInvariant() + '_PROXY'), $selectedProxy, 'Process')
            }
        }
    }
    $outLog = Join-Path $logRoot 'last-sync.out.log'
    $errLog = Join-Path $logRoot 'last-sync.err.log'
    $combinedLog = Join-Path $logRoot 'last-sync.log'
    ('Started: ' + (Get-Date).ToString('o')) | Set-Content -LiteralPath $combinedLog -Encoding UTF8
    $taskFile = Join-Path $projectRoot 'tools\sync_annotation_github.py'
    $process = Start-Process -FilePath $taskPython -ArgumentList @('-u', ('"' + $taskFile + '"')) -WorkingDirectory $projectRoot -WindowStyle Hidden -Wait -PassThru -RedirectStandardOutput $outLog -RedirectStandardError $errLog
    foreach ($logFile in @($outLog,$errLog)) {
        if ((Get-Item -LiteralPath $logFile).Length -gt 0) { Get-Content -LiteralPath $logFile | Add-Content -LiteralPath $combinedLog -Encoding UTF8 }
    }
    ('Exit code: ' + $process.ExitCode) | Add-Content -LiteralPath $combinedLog -Encoding UTF8
    if ($process.ExitCode -ne 0) { throw 'GitHub code sync failed; see .code_publish/last-sync.log' }
} finally { Pop-Location }
