include("../deps/deps.jl")


@self Params function init_reio!()
    z, xedat = get_xe(Ωb, Ωc, ΩΛ, H0, Tcmb, Yp)
    itp = Spline1D(reverse(z), reverse(xedat))
    xe = z->itp(Float64(z))
end


"""
    get_xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
           HOinp::Float64, Tnow::Float64, Yp::Float64, 
           Hswitch::Int64=1, Heswitch::Int64=6, 
           Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)
           
Wrapper of RECFAST Fortran code with parameters as defined in that code.         

Returns tuple of (z's, xe's)
"""
function get_xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
                HOinp::Float64, Tnow::Float64, Yp::Float64;
                Hswitch::Int64=1, Heswitch::Int64=6, 
                Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)
            
    xe = Array{Float64}(undef,Nz)
    ccall( 
        (:get_xe_, librecfast), Nothing, 
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
         Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        OmegaB, OmegaC, OmegaL, HOinp, Tnow, Yp, Hswitch, Heswitch, Nz, zstart, zend, xe
    )
    range(zstart,stop=zend,length=Nz+1)[2:end], xe
end
