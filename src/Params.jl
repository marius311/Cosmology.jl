
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

# fall back which promotes everything to a common numeric type if we tried to
# construct with mixed types
function Params(args...)
    T = promote_type(map(typeof, filter(arg -> arg isa Number, args))...)
    Params(map(arg -> arg isa Number ? T(arg) : arg, args)...)
end


# p broadcasts as a scalar  
Broadcast.broadcastable(p::Params) = Ref(p)

function init_background!(𝕡)
    @set! 𝕡.h² = h² = (𝕡.H0/100)^2
    # photons
    @set! 𝕡.Tγ₀ = 𝕡.Tcmb * Kelvin
    @set! 𝕡.ργ₀ = (π²/15) * 𝕡.Tγ₀^4
    @set! 𝕡.ωγ = 𝕡.ργ₀ / ρx_over_ωx
    @set! 𝕡.Ωγ = 𝕡.ωγ / h²
    # baryons
    @set! 𝕡.ρb₀ = 𝕡.ωb * ρx_over_ωx
    @set! 𝕡.Ωb = 𝕡.ωb / h²
    # CDM
    @set! 𝕡.ρc₀ = 𝕡.ωc * ρx_over_ωx
    @set! 𝕡.Ωc = 𝕡.ωc / h²
    # Neutrinos
    if 𝕡.mν == 0
        @set! 𝕡.Nν_massless += 𝕡.Nν_massive
        @set! 𝕡.Nν_massive = 0
    end
    @set! 𝕡.ρν₀ = ρν(𝕡,0)
    @set! 𝕡.ων = 𝕡.ρν₀ / ρx_over_ωx
    @set! 𝕡.Ων = 𝕡.ων / h²
    # Curvature
    @set! 𝕡.ωk = 𝕡.Ωk * h²
    @set! 𝕡.ρk₀ = 𝕡.ωk * ρx_over_ωx
    # Dark energy
    @set! 𝕡.ΩΛ = 1 - 𝕡.Ωk - 𝕡.Ωb - 𝕡.Ωc - 𝕡.Ων - 𝕡.Ωγ
    @set! 𝕡.ωΛ = 𝕡.ΩΛ * h²
    @set! 𝕡.ρΛ₀ = 𝕡.ωΛ * ρx_over_ωx
end



function new_params(T=Float64; kwargs...)
    𝕡 = Params(;kwargs...)
    𝕡 = init_background!(𝕡)
    𝕡 = init_bbn!(𝕡)
    𝕡 = init_reio!(𝕡)
    𝕡
end
    
