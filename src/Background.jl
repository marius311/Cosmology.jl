# ----------------------
# Input energy densities
# ----------------------
# other quantities are computed based on these


"""Energy density in photons at redshift z"""
@self Params ργ(z) = ργ₀*(1+z)^4

"""Energy density in photons at redshift z"""
@self Params ρΛ(z) = ρΛ₀*(1+z)^(3*(1+w(z)))

"""Energy density in neutrinos at redshift z"""
@self Params ρν(z) = begin
    ρν = Nν_massless*(7/8*(4/11)^(4/3))*ργ₀*(1+z)^4
    if mν != 0
        ρν += Nν_massive*ρ_species(z,mν)
    end
    ρν
end

"""Energy density in cold dark matter at redshift """
@self Params ρc(z) = ρc₀*(1+z)^3

"""Energy density in cold dark matter at redshift z"""
@self Params ρb(z) = ρb₀*(1+z)^3
    
"""Energy density at scale factor a in a thermal species with mass m and 2 d.o.f."""
@self Params ρ_species(z, m) = begin
    a = 1/(1+z)
    a′ = 1e-7
    Tν = Tγ₀*(4/11)^(1/3)
    mT² = (m/(Tν/a′))^2
    integrand(p) = p^2*√((p*(a′/a))^2+mT²)/(exp(√(p^2+mT²))+1)
    (1/π² * (Tν/a′)^4 * (a′/a)^3 * integrate(integrand,0,Inf))
end

# ----


# ------------------------------------
# Compute various angles and distances
# ------------------------------------

"""Hubble constant at redshift z"""
@self Params Hubble(z) = √(8π/3)*√(ρk₀*(1+z)^2 + ρb(z) + ρc(z) + ρΛ(z) + ργ(z) + ρν(z) + ρextra(z))

"""
Θs at the decoupling redshift calculated from zstar_HS. 

This is like CosmoMC's "theta_mc", except CosmoMC's also uses some additional
approximations which we don't use here. 
"""
@self Params Θmc() = Θs(zstar_HS())

"""Angular size of the sound horizon [rad] at redshift z"""
@self Params Θs(z) = rs(z) / DA(z)

"""Conformal time between two redshifts (positive if z2>z1)."""
@self Params η(z1,z2) = integrate(z->1/Hubble(z), z1, z2)

"""Conformal time to redshift z (between -∞ and 0)"""
@self Params η(z) = η(z,0)

"""Comoving distance to redshift z."""
@self Params DC(z) = -η(z)

"""Comoving angular-diameter distance to redshift z."""
@self Params DA(z) = begin
    d = DC(z)
    K = -ωk*(100km/second)^2
    if K==0
        d
    elseif K<0 
        1/√(-K)*sin(d*√(-K))
    elseif K>0
        1/√(K)*sinh(d*√(K))
    end
end

"""Comoving luminosity distance to redshift z"""
@self Params DL(z) = (1+z)^2 * DA(z)

"""Derivative of comoving sound horizon at redshift z"""
@self Params ⅆrs_ⅆz(z) = begin
    R = 3ρb₀/(4ργ₀*(1+z))
    1/Hubble(z)/√(3*(1+R))
end

"""Comoving sound horizon at redshift z"""
@self Params rs(z) = integrate(z->ⅆrs_ⅆz(z), z, Inf)

"""Comoving sound horizon integrated over the visibility function"""
@self Params rs_vis() = integrate(z -> ⅆτ_ⅆz(z) * exp(-τ(z)) * rs(z), 0, Inf)

"""Derivative of square diffusion damping scale at redshift z"""
@self Params ⅆrd²_ⅆz(z) = begin
    R = 3ρb₀/(4ργ₀*(1+z))
    π^2/(σT*ρb₀/mH*(1-Yp)*xe(z)*Hubble(z)*(1+z)^2) * (R^2 + 16/15*(1+R))/(6*(1+R)^2)
end

"""Diffusion damping scale at redshift z"""
@self Params rd(z) = √(integrate(z->ⅆrd²_ⅆz(z), z, Inf))

"""Comoving damping scale integrated over the visibility function."""
@self Params rd_vis(r) = r*√(-log(integrate(z -> ⅆτ_ⅆz(z) * exp(-τ(z)) * exp(-(rd(z)/r)^2), 0, 1e4)))

"""Optical depth between two redshifts"""
@self Params ⅆτ_ⅆz(z) =  σT*ρb₀/mH*(1-Yp) * xe(z)/Hubble(z)*(1+z)^2

"""Optical depth between two redshifts"""
@self Params τ(z1, z2) = integrate(z->ⅆτ_ⅆz(z), z1, z2)

"""Optical depth to redshift z"""
@self Params τ(z) = τ(0, z)

"""Recombination redshift (i.e. z such that τ(z)==1)"""
@self Params zstar() = find_zero(z->τ(z)-1, 800, 1400)

"""Baryon-drag optical depth between two redshifts"""
@self Params τdrag(z1, z2) = begin
    R_div_a = 3ρb₀/4ργ₀
    σT*ρb₀/mH*(1-Yp)/R_div_a * integrate(z->xe(z)/Hubble(z)*(1+z)^3, z1, z2)
end

"""Baryon-drag optical depth to redshift z"""
@self Params τdrag(z) = τdrag(0, z)

"""Baryon-drag reshift (i.e. z such that τdrag(z)==1)"""
@self Params zdrag() = find_zero(z->τdrag(z)-1, 800, 1400)

"""Comoving sound horizon at baryon-drag redshift"""
@self Params rdrag() = rs(zdrag())



# -------------
# Miscellaneous
# -------------

"""Set Θmc (by adjusting H0 accordingly). Returns H0."""
@self Params theta2hubble!(Θ) = find_zero((H0′->(H0=H0′; add_derived!(); Θmc()-Θ)), 20, 200)

"""Redshift at decoupling using fitting formula from Hu & Sugiyama """
@self Params zstar_HS() = 1048*(1+0.00124*ωb^(-0.738))*(1+(0.0783*ωb^(-0.238)/(1+39.5*ωb^0.763))*(ωb+ωc+ων)^(0.560/(1+21.1*ωb^1.81)))
