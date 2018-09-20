module Cosmology

using DelimitedFiles
using Dierckx
using Interpolations
using Libdl
using MacroTools: splitdef, combinedef, postwalk, isexpr, @capture, isdef, splitarg
using Parameters
using PyCall
using QuadGK


export new_params, Params, add_derived!,
       ργ, ρν, ρc, ρ_species, ρx_over_ωx,
       Hubble, Θmc, Θs, D_prop, DA, rs, theta2hubble!, zstar_HS, 
       τ, τd, zdrag, rdrag


include("Units.jl")
include("SelfFunctions.jl")
include("TypeDefaults.jl")


const ρx_over_ωx = 3(100km/second/Mpc)^2/(8π)
const nfac = 7/8*(4/11)^(4/3)
const π² = π^2
const H0units = km/second/Mpc
const Hfac = sqrt(8π/3)


function __init__()
    global brentq = pyimport("scipy.optimize")[:brentq]
end

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


include("BBN.jl")
include("Recfast.jl")


@self Params function quad(f, xmin, xmax)
    quadgk(f, convert(Float64,xmin), convert(Float64,xmax); rtol=reltol)[1]::Float64
end

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
    # init_bbn!(p)
    init_reio!(p)
    p
end
    


# ----------------------
# Input energy densities
# ----------------------
# other quantities are computed based on these


"""Energy density in photons at redshift z"""
@self Params ργ(z) = ργ₀*(1+z)^4

"""Energy density in neutrinos at redshift z"""
@self Params function ρν(z) 
    ρν = Nν_massless*nfac*ργ₀*(1+z)^4
    if mν != 0
        ρν += Nν_massive*ρ_species(z,mν/Nν_massive)
    end
    ρν
end

"""Energy density in cold dark matter at redshift """
@self Params ρc(z) = ρc₀*(1+z)^3

"""Energy density in cold dark matter at redshift z"""
@self Params ρb(z) = ρb₀*(1+z)^3
    
"""Energy density at scale factor a in a thermal species with mass m and 2 d.o.f."""
@self Params function ρ_species(z, m)
    a = 1/(1+z)
    a′ = 1e-7
    Tν = Tγ₀*(4/11)^(1/3)
    mT² = (m/(Tν/a′))^2
    integrand(p) = p^2*sqrt((p*(a′/a))^2+mT²)/(exp(sqrt(p^2+mT²))+1)
    (1/π² * (Tν/a′)^4 * (a′/a)^3 * quad(integrand,0,Inf))::Float64 #todo: why is this type unstable? 
end

# ----


# ------------------------------------
# Compute various angles and distances
# ------------------------------------

"""Hubble constant at redshift z"""
@self Params Hubble(z) = Hfac*sqrt(ρx_over_ωx*((ωc+ωb)*(1+z)^3 + ωk*(1+z)^2 + ωΛ) + ργ(z) + ρν(z))

"""
Θs at the decoupling redshift calculated from zstar_HS. 

This is like CosmoMC's "theta_mc", except CosmoMC's also uses some additional
approximations which we don't use here. 
"""
@self Params Θmc() = Θs(zstar_HS())

"""Angular size of the sound horizon [rad] at redshift z"""
@self Params Θs(z) = rs(z) / DA(z)

"""Conformal time between two redshifts (positive if z2>z1)."""
@self Params η(z1,z2) = quad(z->1/Hubble(z), z1, z2)

"""Conformal time to redshift z (between -∞ and 0)"""
@self Params η(z) = η(z,0)

"""Proper comoving distance to redshift z."""
@self Params dist(z) = -η(z)

"""Comoving angular-diameter distance to redshift z."""
@self Params function DA(z)
    d = dist(z)
    K = -ωk*(100km/second)^2
    if K==0
        d
    elseif K<0 
        1/sqrt(-K)*sin(d*sqrt(-K))
    elseif K>0
        1/sqrt(K)*sinh(d*sqrt(K))
    end
end

"""Comoving sound horizon at redshift z"""
@self Params function rs(z)
    Rovera = 3*ωb*ρx_over_ωx/(4*ργ₀)
    quad(z′ -> 1/Hubble(z′)/sqrt(3*(1+Rovera/(1+z′))), z, Inf)
end

"""Optical depth between two redshifts"""
@self Params function τ(z1, z2)
    σT*(ωb*ρx_over_ωx)/mH*(1-Yp) * quad(z->xe(z)/Hubble(z)*(1+z)^2, z1, z2)
end

"""Optical depth to redshift z"""
@self Params τ(z) = τ(xe, 0, z)

"""Baryon-drag optical depth between two redshifts"""
@self Params function τd(z1, z2)
    Rovera = 3*ωb*ρx_over_ωx/(4*ργ₀)
    σT*(ωb*ρx_over_ωx)/mH*(1-Yp)/Rovera * quad(z->xe(z)/Hubble(z)*(1+z)^3, z1, z2)
end

"""Baryon-drag optical depth to redshift z"""
@self Params τd(z) = τd(0, z)

"""Baryon-drag reshift (i.e. z such that τd(z)==1)"""
@self Params zdrag() = brentq(z->τd(z)-1, 800, 1400)

"""Comoving sound horizon at baryon-drag redshift"""
@self Params rdrag() = rs(zdrag())


# -------------
# Miscellaneous
# -------------

"""Set Θmc (by adjusting H0 accordingly). Returns H0."""
@self Params theta2hubble!(Θ) = brentq((H0′->(H0=H0′; add_derived!(); Θmc()-Θ)), 20, 200, rtol=reltol)

"""Redshift at decoupling using fitting formula from Hu & Sugiyama """
@self Params zstar_HS() = 1048*(1+0.00124*ωb^(-0.738))*(1+(0.0783*ωb^(-0.238)/(1+39.5*ωb^0.763))*(ωb+ωc+ων)^(0.560/(1+21.1*ωb^1.81)))


end
