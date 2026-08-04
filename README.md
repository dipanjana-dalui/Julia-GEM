# GEMs

Before starting with GEMs, please read this document in full.

## Requirements

- **Julia v1.10 or above**
- We recommend using an IDE such as [VS Code](https://code.visualstudio.com/download)

## Getting Started

1. All scripts are saved under the `src` directory.
2. All major function scripts are inside `src/function`.
3. After setting up your environment, run `install_pkgs.jl` **once** the first time to install all required packages. Subsequently, load all packages simultaneously for use in the GEM-run file with the `include()` command.

## Model Definition and Config Files

All definition and config files are prefixed with the model name:

| Prefix | Model |
|--------|-------|
| `bdLM` | Birth-death logistic model |
| `2spp` | Two-species predator-prey model |

> **Note:** Make sure to load the correct definition file together with the matching config file.

## Adding New Functions

New function names should be added to `GEM_function.jl` for repeated use, or included manually in the GEM-run file.

## Using as a Local Package

If loading this as a local package, you only need the scripts under `src/function`. However, you will still need to use the same definition and config formats.

## Multithreading

Set the required number of threads **before** launching Julia.

The easiest way to set the thread count to `n` in VS Code:

1. Open Settings (`Files` menu on Windows, `Code` menu on macOS)
2. Search for `julia threads`
3. Under **Julia: Num Threads**, click **Edit in settings.json**
4. Add or update:

```json
"julia.NumThreads": n
```
