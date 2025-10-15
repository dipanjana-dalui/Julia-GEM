# Julia-GEM 
Julia - [**G**]illespie [**E**]co-evolutionary [**M**]odel

GEMs represent a significant advancement in eco-evolutionary modeling by integrating the effects of trait variation and fitness into ecological interactions of multiple species without relying on fitness gradient assumptions. With GEMs one can incorporate evolutionary feedback directly into ecological dynamics.

## Coded by
Julia Version
- Dipanjana Dalui (ddalui2@unl.edu)
  
MATLAB versions (archived elsewhere)
- John DeLong (jpdelong@unl.edu)



## Installation
For now, Julia-GEM is not a package, and you cannot install it on Julia using the package manager.
However, you may install it as a local package on your system.




VSCode (optional, recommended)

## Quick Start Guide
#### Prerequisites
You will need Julia v1.10 or above. We recommend using
an IDE like VSCode (https://code.visualstudio.com/download)

1. All scripts are saved under src directory.

2. All major function scripts are inside directory src/function.

3. After setting your environment, run install_pkgs.jl the first time to install all needed packages. Subsequently, load all packages  simultaneously for use in the GEM-run file with the include() command.

4. All defintion and config files are prefixed with the model name
    bdLM: birth-death Logistoc model
    2spp: 2 species prey predator model
Make sure to load the correct defintion file with the correct config file.

5. New function names should be added to the GEM_function.jl scipt for repeated use or 
included  manually in the GEM_run file.

6. If loading as a local package, you only need the scripts under src/function. 
However, you will still need to use the same defition and config formats. 

7. Multithreading: set required number of threads BEFORE launching Julia. 
Easiest way to set threads to n (under Files for Windows, under Code for MACOS)
                    Settings > 
                    search "julia threads" > 
                    Julia:Num Threads > 
                    Edit in settings.json > 
                    "julia.NumThreads": n
