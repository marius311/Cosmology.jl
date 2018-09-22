@defaults mutable struct Params{T<:Real}
    #primary parameters
    ωb::T = 0.0225
    ωc::T = 0.12
    H0::T = 67
    Nν_massive::T = 1
    Nν_massless::T = 2.046
    mν::T = 0.06
    Ωk::T = 0
    Tcmb::T = 2.7255
    Yp::T
    xe::Function
    
    #derived
    ρb₀::T; ρc₀::T; ρν₀::T; ργ₀::T; ρk₀::T; ρΛ₀::T
    ων::T; ωγ::T; ωk::T; ωΛ::T
    Ωb::T; Ωc::T; Ων::T; Ωγ::T; ΩΛ::T
    Tγ₀::T
    h²::T
    
    #accuracy
    reltol::T = 1e-4

end

# p broadcasts as a scalar  
Broadcast.broadcastable(p::Params) = Ref(p)

@self Params function init_background!()
    h² = (H0/100)^2
    #photons
    Tγ₀ = (Tcmb*Kelvin)
    ργ₀ = (π²/15)*Tγ₀^4
    ωγ₀ = ργ₀/ρx_over_ωx
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
    p = Params{T}(;kwargs...)
    init_background!(p)
    init_bbn!(p)
    init_reio!(p)
    p
end
    

@self Params{T} function integrate(f, xmin, xmax) where {T}
    quadgk(f, convert(T,xmin), convert(T,xmax); rtol=reltol)[1]::T
end

@self Params{T} function find_zero(f, xmin, xmax) where {T}
    Roots.find_zero(f, (convert(T,xmin), convert(T,xmax)), Roots.A42(); xrtol=reltol)::T
end
