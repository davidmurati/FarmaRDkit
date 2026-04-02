@echo off
title RDKit Pharma Lab
color 0B
echo.
echo  =============================================
echo   RDKit Pharma Lab - Iniciando servidores...
echo  =============================================
echo.

REM Verificar Python
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python no encontrado. Instala Python 3.x primero.
    pause & exit /b 1
)

REM Verificar Node.js
node --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Node.js no encontrado. Instala Node.js primero.
    pause & exit /b 1
)

echo [1/2] Iniciando servidor Python (RDKit API en puerto 5000)...
start "Python RDKit API" cmd /k "python python\server.py"

echo [2/2] Esperando que Python arranque...
timeout /t 3 /nobreak >nul

echo [3/3] Iniciando servidor Node.js (Dashboard en puerto 3000)...
start "Node.js Dashboard" cmd /k "node app.js"

timeout /t 2 /nobreak >nul
echo.
echo  =============================================
echo   Abriendo navegador...
echo  =============================================
start http://localhost:3000

echo.
echo  Servidores en ejecucion:
echo    - Python API:   http://localhost:5000
echo    - Dashboard:    http://localhost:3000
echo.
echo  Cierra las ventanas de consola para detener los servidores.
