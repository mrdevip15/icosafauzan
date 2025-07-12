@echo off
echo ========================================
echo Monte Carlo Simulation Compiler Helper
echo ========================================
echo.
echo Since you don't have a local C compiler, here are your options:
echo.
echo Option 1: Use OnlineGDB (Recommended)
echo 1. Go to: https://www.onlinegdb.com/
echo 2. Copy the contents of icosa8.c
echo 3. Paste into the editor
echo 4. Click "Run" button
echo 5. When prompted, enter: 1000 5000 12345
echo.
echo Option 2: Install MinGW-w64
echo 1. Download from: https://www.mingw-w64.org/
echo 2. Add to PATH
echo 3. Run: gcc -o icosa8 icosa8.c -lm
echo.
echo Option 3: Use WSL (Windows Subsystem for Linux)
echo 1. Install Ubuntu from Microsoft Store
echo 2. Run: sudo apt install gcc
echo 3. Run: gcc -o icosa8 icosa8.c -lm
echo.
pause 