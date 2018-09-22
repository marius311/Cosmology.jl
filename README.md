# Cosmology.jl

*Note: this is a alpha software.*

Cosmology.jl is a code written in [Julia](http://julialang.org/) to compute various cosmological quantities like angular diameter distances, (WIP:) matter power spectra, or (TODO:) the CMB anisotropy spectra.

Requirements:
* Julia 0.7 or higher

Installation,

```julia
pkg> add https://github.com/marius311/SelfFunctions.jl
pkg> add https://github.com/marius311/TypeDefaults.jl
pkg> add https://github.com/marius311/Cosmology.jl
```

Example, 

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

* Cosmology.jl can be **faster**. Don't be fooled by the fact that Julia reads like Python, Julia is [JIT](https://en.wikipedia.org/wiki/Just-in-time_compilation) compiled to native code, so it can be just as fast as C or Fortran, or faster. As a simple example where the codes perform essentially identical mathematical computations, the current calculation of the angular diameter distance is **5x faster** than CAMB. This is encouraging and I see no reason more complex calculations can't have similar performance. 

* The code is far more **readable**. These calculations shouldn't be a black box only a few experts understand. Julia is *built* to write clear scientific code (for example, the aforemention angular diameter distance calculation is roughly **10x fewer lines of code** than either CAMB or CLASS's, and still faster) Here's specifically some of what we make use of:

    * **Unicode variable names**. Julia lets you use special characters in variable names. Cosmology.jl uses this a lot so that variable names immediately look like the equations we are all familiar with. Here's a snippet some simple code as an example,
    
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
    
    * **@self macro** These calculations involve lots of parameters, usually stored in some kind of parameter type, call it `params`. We use a macro `@self` that lets us omit having to write `params.X` every time we need to access parameter `X`, or having to pass `params` around to sub-functions. Its not a huge deal, but does make code a lot less cluttered. For example, instead of this,
    
        ```julia
        function H0(params,z)
            return 8π/3*sqrt(params.ργ₀*(1+z)^4 + params.ρm₀*(1+z)^4 + ρν(params,z) + params.ρΛ)
        end
        ``` 
        
        we have this,
        
        ```julia
        @self Params function H0(z)
            return 8π/3*sqrt(ργ₀*(1+z)^4 + ρm₀*(1+z)^4 + ρν(z) + ρΛ)
        end
        ```
        
        Macros are just one of the many language features of Julia that lets us write great scientific code.
    
* Opens up the possibility of **CPU vectorization**, **GPU acceleration** and **automatic differentiation**. In Julia you can easily use the vectorization features of your processor by prepending code with the `@simd` macro. You can also use any of the GPU programming packages like [ArrayFire](https://github.com/JuliaComputing/ArrayFire.jl) to run computations on your GPU. Either of these things would take lots of complex C / Fortran code, but can work almost transparently with Julia. 

    Additionally, Julia has several [automatic differentiation](http://www.juliadiff.org/) (AD) packages. AD means that you write your function, say `f(x)`, and Julia automatically writes the corresponding *exact* derivative of this function `df(x)/dx`, for you. This can be done for many complex functions `f` even if they contain if statements, for loops, etc... This is different than computing finite differences; these are *exact* derivatives, and they take only O(1) times longer than the original function to compute. These can then be used in conjunction with more sophisticated ODE solvers or quadrature function to make the code even faster. 

    These three things aren't currently in use in Cosmology.jl, but potentially could be resulting in even more drastic speed improvements in the future. 
