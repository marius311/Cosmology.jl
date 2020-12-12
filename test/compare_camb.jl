
using Test, PyCall, BenchmarkTools, PyPlot, Cosmology, Libdl, Printf
@pyimport camb

# CAMB parameter set
camb_res     = camb.get_background(camb.set_params(H0=67, ombh2=0.0225, omch2=0.12, mnu=0.06, nnu=3.046, tau=0.07, YHe=0.25))
camb_derived = camb_res[:get_derived_params]();

# Our parameter set
p = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0.06eV, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)


# Run correctness tests
# ---------------------

@testset begin
    
    z = camb_derived["zstar"]
    @test DA(p,z)/Mpc ≈ camb_res[:angular_diameter_distance](z)*(1+z)   rtol = 1e-4
    @test rs(p,z)/Mpc ≈ camb_derived["rstar"]                           rtol = 1e-4
    @test θs(p,z)     ≈ camb_derived["thetastar"]/100                   rtol = 1e-4
    @test zdrag(p)    ≈ camb_derived["zdrag"]                           rtol = 1e-4
    
end


# Run performance tests
# ---------------------

# set mν=0 for speed comparison (we can't directly compare when mν!=0 since CAMB
# uses some approximations that we don't)
cp = camb.get_background(camb.set_params(H0=67, ombh2=0.0225, omch2=0.12, mnu=0, nnu=3.046, tau=0.07, YHe=0.25))
p = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)

tcamb = @belapsed cp.angular_diameter_distance(1100.) # todo: ccall directly into the camb lib to remove any (small) Python overhead
tus = @belapsed DA($p,1100.)

println("Performace test (speedup compared to CAMB)")
@printf "DA - %.2fx\n" tcamb/tus
