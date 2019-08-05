module CompareCamb

push!(LOAD_PATH,"../src")
using Base.Test, PyCall, TimeIt, PyPlot, Cosmology, Units
@pyimport camb
push!(Base.Libdl.DL_LOAD_PATH,dirname(camb.camb[:camblib][:_name]))


# CAMB parameter set
camb_res     = camb.get_background(camb.set_params(H0=67,ombh2=0.0225,omch2=0.12,mnu=0.06,nnu=3.046,tau=0.07,YHe=0.25))
camb_derived = camb_res[:get_derived_params]();

# Our parameter set
const p = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0.06eV, Nν_massive=1, Nν_massless=2.046, reltol=1e-4)

# Acceptable fractional difference between us and CAMB
const ftol = 1e-4
macro equalish(a,b)
    :(abs($(esc(a))/$(esc(b))-1)<ftol)
end

# Run correctness tests
# ---------------------

@testset begin
    
    z = camb_derived["zstar"]
    @test @equalish DA(p,z)/Mpc camb_res[:angular_diameter_distance](z)*(1+z)
    @test @equalish rs(p,z)/Mpc camb_derived["rstar"]
    @test @equalish θs(p,z)     camb_derived["thetastar"]/100
    @test @equalish zdrag(p)    camb_derived["zdrag"] #is failing with ftol=1e-4 though
    
end


# Run performance tests
# ---------------------

# set mν=0 for speed comparison (we can't directly compare when mν!=0 since CAMB
# uses some approximations that we don't)
camb.get_background(camb.set_params(H0=67,ombh2=0.0225,omch2=0.12,mnu=0,nnu=3.046,tau=0.07,YHe=0.25))
const p2 = new_params(;H0=67, ωb=0.0225, ωc=0.12, mν=0., Nν_massive=1, Nν_massless=2.046, reltol=1e-4)

#todo: the function name assumes we compiled CAMB with gfortan:
tcamb = @timeit ccall((:__modelparams_MOD_angulardiameterdistance, "camblib"), Float64, (Ref{Float64},), 1100.)
tus = @timeit DA(p2,1100.)

println("Performace test (speedup compared to CAMB)")
@printf "DA - %.1fx\n" tcamb/tus

end
