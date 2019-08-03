@kwdef struct Params{T<:Real,Fxe,Fw,Fρextra,Fρνmassive}
    # primary parameters
    ωb::T = 0.0225
    ωc::T = 0.12
    H0::T = 67.0
    Nν_massive::T = 1.0
    Nν_massless::T = 2.046
    mν::T = 0.06
    Ωk::T = 0.0
    Tcmb::T = 2.7255
    Yp::T = NaN
    xe::Fxe = nothing
    w::Fw = (z->-1)
    ρextra::Fρextra = (z->0)
    
    # precomputable functions
    ρνmassive::Fρνmassive = ρνmassive
    
    # derived
    ρb₀::T=NaN; ρc₀::T=NaN; ρν₀::T=NaN; ργ₀::T=NaN; ρk₀::T=NaN; ρΛ₀::T=NaN
    ων::T=NaN; ωγ::T=NaN; ωk::T=NaN; ωΛ::T=NaN
    Ωb::T=NaN; Ωc::T=NaN; Ων::T=NaN; Ωγ::T=NaN; ΩΛ::T=NaN
    Tγ₀::T=NaN
    h²::T=NaN
    
    # accuracy
    reltol::T = 1e-4

end

# p broadcasts as a scalar  
Broadcast.broadcastable(p::Params) = Ref(p)

@self Params function init_background!()
    h² = (H0/100)^2
    #photons
    Tγ₀ = (Tcmb*Kelvin)
    ργ₀ = (π²/15)*Tγ₀^4
    ωγ = ργ₀/ρx_over_ωx
    Ωγ = ωγ/h²
    #baryons
    ρb₀ = ωb*ρx_over_ωx
    Ωb = ωb/h²
    #CDM
    ρc₀ = ωc*ρx_over_ωx
    Ωc = ωc/h²
    #Neutrinos
    if mν == 0
        Nν_massless += Nν_massive
        Nν_massive = 0
    end
    ρν₀ = ρν(0)
    ων = ρν₀/ρx_over_ωx
    Ων = ων/h²
    #Curvature
    ωk = Ωk*h²
    ρk₀ = ωk*ρx_over_ωx
    #Dark energy
    ΩΛ = 1 - Ωk - Ωb - Ωc - Ων - Ωγ
    ωΛ = ΩΛ*h²
    ρΛ₀ = ωΛ*ρx_over_ωx
end



function new_params(T=Float64;kwargs...)
    p = Params(;kwargs...)
    p = init_background!(p)
    p = init_bbn!(p)
    p = init_reio!(p)
    p
end
    

@self Params{T} function integrate(f, xmin, xmax) where {T}
    quadgk(f, convert(T,xmin), convert(T,xmax); rtol=reltol)[1]::T
end

@self Params{T} function find_zero(f, xmin, xmax) where {T}
    Roots.find_zero(f, (convert(T,xmin), convert(T,xmax)), Roots.A42(); xrtol=reltol)::T
end
