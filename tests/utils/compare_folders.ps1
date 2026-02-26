# PowerShell script to compare the contents of moleditpy/src/moleditpy and moleditpy-linux/src/moleditpy_linux

$root1 = Join-Path (Get-Item .).FullName "moleditpy/src/moleditpy"
$root2 = Join-Path (Get-Item .).FullName "moleditpy-linux/src/moleditpy_linux"

if (!(Test-Path $root1)) {
    Write-Host "Error: $root1 not found." -ForegroundColor Red
    exit
}

if (!(Test-Path $root2)) {
    Write-Host "Error: $root2 not found." -ForegroundColor Red
    exit
}

Write-Host "Comparing SRC contents:" -ForegroundColor Cyan
Write-Host "1: $root1"
Write-Host "2: $root2"
Write-Host ""

# Get all relative paths in both folders
function Get-RelativePaths($root) {
    Get-ChildItem $root -Recurse | ForEach-Object {
        $_.FullName.Substring($root.Length + 1)
    }
}

$paths1 = Get-RelativePaths $root1
$paths2 = Get-RelativePaths $root2

# Compare existence
$diff = Compare-Object $paths1 $paths2

if ($null -eq $diff) {
    Write-Host "File and folder structure is identical." -ForegroundColor Green
} else {
    Write-Host "Structure differences found:" -ForegroundColor Yellow
    $diff | ForEach-Object {
        $side = if ($_.SideIndicator -eq "<=") { "Only in moleditpy" } else { "Only in moleditpy-linux" }
        Write-Host "  [$side] $($_.InputObject)"
    }
}

Write-Host "`nChecking for content differences in common files..." -ForegroundColor Cyan
$commonFiles = $paths1 | Where-Object { $paths2 -contains $_ }

$diffCount = 0
foreach ($relPath in $commonFiles) {
    $file1 = Join-Path $root1 $relPath
    $file2 = Join-Path $root2 $relPath
    
    # Skip directories
    if ((Get-Item $file1).PSIsContainer) { continue }

    $hash1 = Get-FileHash $file1
    $hash2 = Get-FileHash $file2
    
    if ($hash1.Hash -ne $hash2.Hash) {
        Write-Host " [DIFF] $relPath" -ForegroundColor Yellow
        $diffCount++
    }
}

if ($diffCount -eq 0) {
    Write-Host "All common files have identical content." -ForegroundColor Green
} else {
    Write-Host "`nTotal files with differences: $diffCount" -ForegroundColor Yellow
}

Write-Host "`nComparison complete." -ForegroundColor Cyan
