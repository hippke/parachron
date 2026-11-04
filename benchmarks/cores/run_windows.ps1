# Parachron core-scaling benchmark for Windows testers: installs Parachron into $HOME\parachron-test (a separate Python
# environment), downloads bench_cores.py, and runs the quick and then the full benchmark. See WINDOWS.md.
#
# Usage (PowerShell):
#   powershell -ExecutionPolicy Bypass -File .\run_windows.ps1            # quick + full test (about 10 minutes)
#   powershell -ExecutionPolicy Bypass -File .\run_windows.ps1 -Quick     # quick test only (about 2 minutes)
# Results: $HOME\parachron-test\parachron_bench_*.txt (and setup_log.txt if the installation fails).
param(
    [switch]$Quick,
    [string]$Package = "parachron",     # what pip installs (a file path for tests of unreleased versions)
    [string]$Script = "https://raw.githubusercontent.com/hippke/parachron/main/benchmarks/cores/bench_cores.py"
)

function Fail($msg) {
    Write-Host ""
    Write-Host "ERROR: $msg" -ForegroundColor Red
    Write-Host "Please send us the file $dir\setup_log.txt" -ForegroundColor Red
    Stop-Transcript | Out-Null
    exit 1
}

$dir = Join-Path $HOME "parachron-test"
New-Item -ItemType Directory -Force $dir | Out-Null
Set-Location $dir
Start-Transcript -Path (Join-Path $dir "setup_log.txt") | Out-Null
Write-Host "Parachron benchmark: working in $dir"

# --- Python: the py launcher of python.org, else python on the PATH (not the Microsoft Store placeholder)
$py = $null
if (Get-Command py -ErrorAction SilentlyContinue) { $py = @("py", "-3") }
elseif ((Get-Command python -ErrorAction SilentlyContinue) -and ((& python --version 2>&1) -match "Python 3")) { $py = @("python") }
if (-not $py) { Fail "Python 3 was not found. Install it from https://www.python.org/downloads/ (tick 'Add python.exe to PATH'), then open a new PowerShell window and run this script again." }
$pyexe = $py[0]; $pyargs = @($py | Select-Object -Skip 1)
Write-Host ("Python: " + (& $pyexe @pyargs --version 2>&1))

# --- separate environment and installation (compiles Parachron with clang-cl)
if (-not (Test-Path "venv\Scripts\python.exe")) {
    & $pyexe @pyargs -m venv venv
    if ($LASTEXITCODE -ne 0) { Fail "Could not create the Python environment." }
}
$vpy = Join-Path $dir "venv\Scripts\python.exe"
& $vpy -m pip install --upgrade pip
Write-Host "Installing Parachron (compiles for 2-5 minutes; warnings are normal) ..."
& $vpy -m pip install --upgrade --no-cache-dir $Package
if ($LASTEXITCODE -ne 0) { Fail "The installation of Parachron failed. Is 'C++ Clang tools for Windows' installed (see WINDOWS.md, step 3)?" }

# --- benchmark
& curl.exe -fsSL -o bench_cores.py $Script
if ($LASTEXITCODE -ne 0 -or -not (Test-Path bench_cores.py)) { Fail "Could not download the benchmark script." }
Stop-Transcript | Out-Null

Write-Host ""
Write-Host "Quick test (about 2 minutes) ..." -ForegroundColor Green
& $vpy bench_cores.py --quick
if ($LASTEXITCODE -ne 0) { Write-Host "The quick test failed; please send us the parachron_bench_*.txt file." -ForegroundColor Red; exit 1 }
if (-not $Quick) {
    Write-Host ""
    Write-Host "Full test (5-10 minutes). Please do not use the PC until it has finished." -ForegroundColor Green
    & $vpy bench_cores.py
}
Write-Host ""
Write-Host "Finished. Please send us all files parachron_bench_*.txt from $dir" -ForegroundColor Green
if (-not $env:CI) { explorer.exe $dir }
