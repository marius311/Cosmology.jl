# Cosmology.jl

*Note: this is a alpha software.*

Cosmology.jl is a code written in [Julia](http://julialang.org/) to compute various cosmological quantities like angular diameter distances, (WIP:) matter power spectra, or (TODO:) the CMB anisotropy spectra.

Requirements:
* Julia 1.5 or higher

Installation (from the package prompt):

```julia
pkg> add https://github.com/marius311/Cosmology.jl
```

Example:

```julia
julia> using Cosmology

julia> p = new_params(H0=67, ωb=0.0225, ωc=0.12, mν=0.06eV, Nν_massive=1, Nν_massless=2.046, reltol=1e-4);

julia> rs(p,zdrag(p))/Mpc # sound horizon at baryon drag redshift
146.94167579179634

# etc..
```

Things currently handled correctly:
* Background 
    * Massive neutrinos
    * BBN (via same interpolation tables as CAMB)
    * Reionization (comes bundled with Recfast)


## Motivation

The field of Cosmology already has some pretty mature and widely used codes which perform these calculations, [CAMB](camb.info) and [CLASS](class-code.net). So why Cosmology.jl?

* Cosmology.jl can be **faster**. Julia is [JIT](https://en.wikipedia.org/wiki/Just-in-time_compilation) compiled to native code, so it can be just as fast as C or Fortran, or faster. As a simple example where the codes perform essentially identical mathematical computations, the current calculation of the angular diameter distance is **7x faster than CAMB**. This is encouraging and I see no reason more complex calculations can't have similar performance. 

* The code is far more **readable**. These calculations shouldn't be a black box only a few experts understand. Julia is built to write clear scientific code (for example, the aforemention angular diameter distance calculation is roughly **10x fewer lines of code** than either CAMB or CLASS's, and still faster). A few specific things that really help:

    * **Unicode variable names**. This may seem silly but it makes a big differnce in reading code when you're familiar with the equations but not with the code-base. And Julia has many more Unicode characters available than Python, so this works really well. E.g., a code snippet:
    
        ```julia
        function add_derived!()
            ργ₀ = (π²/15)*Tγ₀^4
            h² = (H0/100)^2
            ωk = Ωk*h²
            if mν == 0
                Nν_massless += Nν_massive
                Nν_massive = 0
            end
            ων = ρν(0)/ρx_over_ωx
            ωΛ = h² - ωk - ωb - ωc - ων - ργ₀/ρx_over_ωx
        end
        ```
    
    * **DAE solver**. Lots of these calculations involve solving ordinary differential equations (ODEs). Traditional ODE solvers require you to manipulate your differential equations into the form `y′ = F(y, x)`, which often puts them in an unfamiliar form far removed from the physically motivated way in which they are usually written down. We use instead a "differential algebraic solver" (DAE) which lets you write any general system of differential equations `F(y′,y,x)=0`. This (coupled with Unicode variable names) allows for writing very clear code, for example from the matter power spectrum calculation:
    
        ```julia
        @eqs begin
            #Photons 
            δγ′ == -4/3*θγ + 4ϕ′
            θγ′ == k²/4*δγ + k²*ψ + a*nₑ*σT*(θb-θγ)
            
            #CDM
            δc′ == -θc + 3ϕ′
            θc′ == -(a′/a)*θc + k²*ψ
            
            #Baryons
            δb′ == -θb + 3ϕ′
            θb′ == -(a′/a)*θb + cₛ²*k²*δb + a*nₑ*σT/R*(θγ-θb) + k²*ψ
            
            #Einstein equations
            k²*ϕ + 3*(a′/a)*(ϕ′+(a′/a)*ψ) == -4π*a^2*(ρc(z)*δc + ρb(z)*δb + (ργ(z)+ρν(z))*δγ)
        end
        ```
    
    * **Planck units** Absolutely every quantity in Cosmology.jl is in Planck units, so there's never any confusion, and equations are as simple as possible.  (Other codes don't do this in an attempt to keep native 32bit support, since 32bit floats don't have a big enough exponent for expressing cosmological quantities in Planck units; we drop native 32bit support, i.e. you can run on 32bit, but its slower)
    
* Opens up the possibility **GPU acceleration** and **automatic differentiation**. Julia has great support for automatic differentiation, which works through almost any arbitrary Julia code. Something like [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) should work almost out-of-the-box and provide exact gradients in O(1) times the cost of running the calculation itself. Not only can this be used to compute Fisher matrices or for analysis techniques that require likelihood gradients (HMC, VI, etc...), but exact gradients can also be used to speed up the internal ODE solvers too. Additionally, Julia has great GPU integration, and using something like [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) or [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) one could offload much of the computation to GPU almost transparently, and without needing to write anything beside Julia code. 
