Before starting GEMs, please read this document.
You will need Julia v1.10 or above. We recommend using
an IDE like VSCode (https://code.visualstudio.com/download)

1. All scripts are saved under src directory.
2. All major function scripts are inside directory src/function.
3. After setting your environment, run install_pkgs.jl the first 
time to install all needed packages. Subsequently, load all packages 
simultaneously for use in the GEM-run file with the include() command.
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