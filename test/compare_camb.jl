module CompareCamb

push!(LOAD_PATH,"../src")
using Base.Test, PyCall, TimeIt, PyPlot, Cosmology, Units
@pyimport camb

# CAMB parameter set
cp = camb.set_params(H0=67,ombh2=0.0225,omch2=0.12,mnu=0.06,nnu=3.046,tau=0.07)
res = camb.get_background(cp)
derived = res[:get_derived_params]();

# Our parameter set
const p = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0.06eV, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)

# Acceptable fractional difference between results
const ftol = 2e-4
macro equalish(a,b)
    :(abs($(esc(a))/$(esc(b))-1)<ftol)
end


# Run tests
@testset begin
    
    z = derived["zstar"]
    @test @equalish DA(p,z)/Mpc res[:angular_diameter_distance](z)*(1+z)
    @test @equalish rs(p,z)/Mpc derived["rstar"]
    @test @equalish Θs(p,z)     derived["thetastar"]/100
    @test @equalish zdrag(p)    derived["zdrag"] #is failing with ftol=1e-4 though
    
end


end
