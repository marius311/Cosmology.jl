
using Test, PyCall, BenchmarkTools, Cosmology, Libdl, Printf, Libdl
@pyimport camb
@assert VersionNumber(camb.__version__) >= v"1.1.3" "Need CAMB >= 1.1.3 for this test."
push!(Libdl.DL_LOAD_PATH,dirname(camb.camb[:camblib][:_name]))
py"""
from ctypes import byref, cast, c_void_p
"""

# CAMB parameter set
camb_res     = camb.get_background(camb.set_params(H0=67, ombh2=0.0225, omch2=0.12, mnu=0.06, nnu=3.046, tau=0.07, YHe=0.25))
camb_derived = camb_res[:get_derived_params]();

# Our parameter set
𝕡 = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0.06eV, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)


# Run correctness tests
# ---------------------

@testset begin
    
    z = camb_derived["zstar"]
    @test DA(𝕡,z)/Mpc ≈ camb_res[:angular_diameter_distance](z)*(1+z)   rtol = 1e-4
    @test rs(𝕡,z)/Mpc ≈ camb_derived["rstar"]                           rtol = 1e-4
    @test θs(𝕡,z)     ≈ camb_derived["thetastar"]/100                   rtol = 1e-4
    @test zdrag(𝕡)    ≈ camb_derived["zdrag"]                           rtol = 1e-4
    
end


# Run performance tests
# ---------------------

# set mν=0 for speed comparison (we can't directly compare when mν!=0 since CAMB
# uses some approximations that we don't)
cp = camb.get_background(camb.set_params(H0=67, ombh2=0.0225, omch2=0.12, mnu=0, nnu=3.046, tau=0.07, YHe=0.25))
𝕡 = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)

tus = @belapsed DA($𝕡, 1100.)
# ccall directly into the Fortran library to avoid any overhead
tcamb = @belapsed ccall(
    (:__results_MOD_cambdata_angulardiameterdistance, "camblib"), 
    Float64, 
    (Ptr{Cvoid}, Ref{Float64},),
    $(py"cast(byref($cp.fortran_self), c_void_p)"),
    1100.
)

println("Performace test (speedup compared to CAMB)")
@printf "DA - %.2fx\n" tcamb/tus
