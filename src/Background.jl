# ----------------------
# Input energy densities
# ----------------------
# other quantities are computed based on these


"""Energy density in photons at redshift z"""
ργ(𝕡, z) = 𝕡.ργ₀ * (1+z)^4

"""Equation of state parameter"""
w(𝕡, z) = 𝕡.w(z)

"""Energy density in dark energy at redshift z"""
ρΛ(𝕡, z) = 𝕡.ρΛ₀ * (1+z)^(3*(1+w(𝕡,z)))

"""Energy density in cold dark matter at redshift z"""
ρc(𝕡, z) = 𝕡.ρc₀ * (1+z)^3

"""Energy density in baryons at redshift z"""
ρb(𝕡, z) = 𝕡.ρb₀ * (1+z)^3

"""Extra customizable energy"""
ρextra(𝕡, z) = 𝕡.ρextra(z)

"""Energy density in massless+massive neutrinos at redshift z"""
ρν(𝕡, z) = ρνmassless(𝕡,z) + ρνmassive(𝕡,z)

"""Energy density in massless neutrinos at redshift z"""
ρνmassless(𝕡, z) = 𝕡.Nν_massless * (7/8*(4/11)^(4/3)) * 𝕡.ργ₀ * (1+z)^4

"""Energy density in massive neutrinos at redshift z"""
function ρνmassive(𝕡, z)
    if 𝕡.mν == 0
        0
    else
        a = 1/(1+z)
        a′ = 1e-7
        Tν = 𝕡.Tγ₀*(4/11)^(1/3)
        mT² = (𝕡.mν/(Tν/a′))^2
        integrand(p) = p^2*√((p*(a′/a))^2+mT²)/(exp(√(p^2+mT²))+1)
        (𝕡.Nν_massive * 1/π² * (Tν/a′)^4 * (a′/a)^3 * integrate(𝕡,integrand,0,Inf))
    end
end

function precompute_ρνmassive(𝕡, p, z = 10 .^ range(-3,10,length=128))
    spl = Spline1D(log.(z), log.(ρνmassive.(𝕡,p,z)))
    @set p.ρνmassive = z -> exp(spl(log(z)))
end

# ----



# ------------------------------------
# Compute various angles and distances
# ------------------------------------

"""Hubble constant at redshift z"""
Hubble(𝕡, z) = √(8π/3)*√(𝕡.ρk₀*(1+z)^2 + ρb(𝕡,z) + ρc(𝕡,z) + ρΛ(𝕡,z) + ργ(𝕡,z) + ρν(𝕡,z) + ρextra(𝕡,z))

"""
θs at the decoupling redshift calculated from zstar_HS. 

This is like CosmoMC's "theta_mc", except CosmoMC's also uses some additional
approximations which we don't use here. 
"""
θmc(𝕡) = θs(𝕡,zstar_HS(𝕡))

"""Angular size of the sound horizon [rad] at redshift z"""
θs(𝕡, z) = rs(𝕡,z) / DA(𝕡,z)

"""Conformal time between two redshifts (positive if z2>z1)."""
η(𝕡, z1, z2) = integrate(𝕡, z->1/Hubble(𝕡,z), z1, z2)

"""Conformal time to redshift z (between -∞ and 0)"""
η(𝕡, z) = η(𝕡, z, 0)

"""Comoving distance to redshift z."""
DC(𝕡, z) = -η(𝕡, z)

"""Comoving angular-diameter distance to redshift z."""
function DA(𝕡, z)
    d = DC(𝕡, z)
    K = -𝕡.ωk*(100km/second)^2
    if K==0
        d
    elseif K<0 
        1/√(-K)*sin(d*√(-K))
    elseif K>0
        1/√(K)*sinh(d*√(K))
    end
end

"""Comoving luminosity distance to redshift z"""
DL(𝕡, z) = (1+z)^2 * DA(𝕡,z)

"""Derivative of comoving sound horizon at redshift z"""
function ⅆrs_ⅆz(𝕡, z)
    R = 3𝕡.ρb₀/(4𝕡.ργ₀*(1+z))
    1/Hubble(𝕡,z)/√(3*(1+R))
end

"""Comoving sound horizon at redshift z"""
rs(𝕡, z) = integrate(𝕡, z->ⅆrs_ⅆz(𝕡,z), z, Inf)

"""Comoving sound horizon integrated over the visibility function"""
rs_vis(𝕡) = integrate(𝕡, z->ⅆτ_ⅆz(𝕡,z)*exp(-τ(𝕡,z))*rs(𝕡,z), 0, Inf)

"""Derivative of square diffusion damping scale at redshift z"""
function ⅆrd²_ⅆz(𝕡, z)
    R = 3𝕡.ρb₀ / (4𝕡.ργ₀*(1+z))
    π^2/(σT*𝕡.ρb₀/mH*(1-𝕡.Yp)*xe(𝕡,z)*Hubble(𝕡,z)*(1+z)^2) * (R^2 + 16/15*(1+R))/(6*(1+R)^2)
end

"""Diffusion damping scale at redshift z"""
rd(𝕡, z) = √(integrate(𝕡, z->ⅆrd²_ⅆz(𝕡,z), z, Inf))

"""Comoving damping scale integrated over the visibility function."""
rd_vis(𝕡, r) = r*√(-log(integrate(𝕡, z->ⅆτ_ⅆz(𝕡,z)*exp(-τ(𝕡,z))*exp(-(rd(𝕡,z)/r)^2), 0, 1e4)))

"""Optical depth between two redshifts"""
ⅆτ_ⅆz(𝕡, z) =  σT*𝕡.ρb₀/mH*(1-𝕡.Yp) * xe(𝕡,z)/Hubble(𝕡,z)*(1+z)^2

"""Optical depth between two redshifts"""
τ(𝕡, z1, z2) = integrate(𝕡, z->ⅆτ_ⅆz(𝕡,z), z1, z2)

"""Optical depth to redshift z"""
τ(𝕡, z) = τ(𝕡, 0, z)

"""Recombination redshift (i.e. z such that τ(z)==1)"""
zstar(𝕡) = find_zero(𝕡, z->τ(𝕡,z)-1, 800, 1400)

"""Ionization fraction at redshift z"""
xe(𝕡,z) = 𝕡.xe(z)

"""Baryon-drag optical depth between two redshifts"""
function τdrag(𝕡, z1, z2)
    R_div_a = 3𝕡.ρb₀/4𝕡.ργ₀
    σT*𝕡.ρb₀/mH*(1-𝕡.Yp)/R_div_a * integrate(𝕡, z->xe(𝕡,z)/Hubble(𝕡,z)*(1+z)^3, z1, z2)
end

"""Baryon-drag optical depth to redshift z"""
τdrag(𝕡, z) = τdrag(𝕡, 0, z)

"""Baryon-drag reshift (i.e. z such that τdrag(z)==1)"""
zdrag(𝕡) = find_zero(𝕡, z->τdrag(𝕡,z)-1, 800, 1400)

"""Comoving sound horizon at baryon-drag redshift"""
rdrag(𝕡) = rs(𝕡, zdrag(𝕡))


# -------------
# Miscellaneous
# -------------

"""Set θmc (by adjusting H0 accordingly). Returns H0."""
theta2hubble!(𝕡,θ) = find_zero(𝕡, (H0′->(@set!(𝕡.H0=H0′); init_background!(𝕡); θmc(𝕡)-θ)), 20, 200)

"""Redshift at decoupling using fitting formula from Hu & Sugiyama """
zstar_HS(𝕡) = 1048*(1+0.00124*𝕡.ωb^(-0.738))*(1+(0.0783*𝕡.ωb^(-0.238)/(1+39.5*𝕡.ωb^0.763))*(𝕡.ωb+𝕡.ωc+𝕡.ων)^(0.560/(1+21.1*𝕡.ωb^1.81)))
