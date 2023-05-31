# SailorMoon.jl
Translunar trajectory design in Julia

### Dependencies

- Public-packages: `DifferentialEquations.jl`, `Plots.jl`
- In-house packages: [`AstrodynamicsBase.jl`](https://github.com/Yuricst/AstrodynamicsBase.jl), [`joptimise.jl`](https://github.com/Yuricst/joptimise)

### Connection to `joptimise.jl` (SNOPT) (Last Edit: 05/31/2023 Yuji)

The functionality of SNOPT is confirmed by using a slightly modified version of `joptimise.jl`, which is currently at `for_SailorMoon` branch in [this link](https://github.com/UzTak/joptimise/tree/for_SailorMoon). 
The environment/setup used here is as follows:

- OS: WSL (Ubuntu), windows 11 
- Use the **Fortran/C/C++** libraries but not **Fortran/C** libraries; the [download page](https://ccom.ucsd.edu/~optimizers/downloads/software/academic/?id=8c697396914c) of SNOPT mentions that Linux needs C++ files, suffixed with "XXX_cpp".

Few tips:

- Make sure to put the SNOPT files in joptimise/src/snopt.
- Small values of `lencw` (default `joptimise.jl` set it to 500) can cause a segmentation fault / memory error. It looks like somewhere like 5000 works well. 
- Change l.17 of `joptimise/src/snopt.jl` so that `_cpp` files are loaded. 


<p align="center">
    <img src="./etc/sailormoon.gif" width="550" title="smg">
</p>


