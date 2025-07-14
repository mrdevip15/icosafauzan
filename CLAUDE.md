# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Monte Carlo simulation project for studying 2D Ising models with 8-state Potts spins (icosahedral spins). The project implements statistical physics simulations to analyze critical phenomena and phase transitions in 2D spin systems.

## Core Components

### C Code (Simulation Engine)
- `icosa8.c` - Main Monte Carlo simulation using Wolff cluster algorithm
- `rn32.c` - Random number generator implementation  
- `cubic3lay.c` - Additional lattice model (untracked file)

### Python Analysis Pipeline
- `run_simulation.py` - Single simulation runner and data capture
- `run_multiple_simulations_parallel.py` - Parallel execution across multiple CPU cores
- `analyze_thermodynamics.py` - Comprehensive thermodynamic analysis and plotting
- `create_custom_range_plots.py` - Custom temperature range analysis

## Essential Commands

### Compilation
```bash
gcc -o icosa8 icosa8.c rn32.c -lm
```

### Python Environment Setup
```bash
pip install -r requirements.txt
```

### Run Complete Analysis Pipeline
```bash
python run_simulation.py
```

### Run Parallel Simulations (Production)
```bash
python run_multiple_simulations_parallel.py
```

### Manual Simulation
```bash
echo "8 8 1000 10000 12345" | ./icosa8.exe > output.txt
```

## Input Format
The C simulation expects 5 integer parameters via stdin:
- nx: lattice width
- ny: lattice height  
- nmcs1: equilibration steps
- nmcs2: measurement steps
- seed: random number generator seed

## Key Architecture Patterns

### Data Flow
1. C simulation generates raw thermodynamic data
2. Python scripts capture and parse simulation output
3. ThermodynamicAnalyzer class processes data into publication-ready plots
4. Multiple lattice sizes are analyzed for finite-size scaling

### Output Structure
- Raw simulation data in `simulation_runs/` directory structure
- Organized by lattice size (nx8_ny8, nx16_ny16, etc.)
- Multiple runs per configuration for statistical analysis
- Publication-quality PNG plots and detailed text reports

### Parallel Processing
The parallel simulation runner uses all available CPU cores to execute multiple Monte Carlo runs simultaneously, significantly reducing computation time for production analysis.

### Scientific Focus
This codebase is designed for academic/research use in statistical physics, specifically for studying:
- Critical temperature determination using Binder cumulant method
- Critical exponent extraction (β, γ)
- Finite-size scaling analysis
- Universality class identification for 2D Ising-like systems

## Dependencies
- C compiler with math library (-lm)
- Python 3 with numpy, matplotlib, pandas, scipy, seaborn
- Executable must be named `icosa8.exe` for Python scripts to work