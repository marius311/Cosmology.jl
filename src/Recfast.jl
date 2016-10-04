module Recfast

using Cosmology, SelfFunctions, Interpolations

export init_reio!

push!(Base.Libdl.DL_LOAD_PATH,dirname(Base.source_path()))


@self Params function init_reio!()
    zdat, xedat = get_xe(Ωb, Ωc, ΩΛ, H0, Tcmb, Yp)
    xe = (z::Float64) -> interpolate((reverse!(zdat),),reverse!(xedat),Gridded(Linear()))[z]
end


doc"""
    get_xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
           HOinp::Float64, Tnow::Float64, Yp::Float64, 
           Hswitch::Int64=0, Heswitch::Int64=0, 
           Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)
           
Wrapper of RECFAST Fortran code with parameters as defined in that code.         

Returns tuple of (z's, xe's)
"""
function get_xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
                HOinp::Float64, Tnow::Float64, Yp::Float64, 
                Hswitch::Int64=0, Heswitch::Int64=0, 
                Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)
            
    xe = Array{Float64}(Nz)
    ccall( 
        (:get_xe_, "recfast"), Void, 
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
         Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        OmegaB, OmegaC, OmegaL, HOinp, Tnow, Yp, Hswitch, Heswitch, Nz, zstart, zend, xe
    )
    collect(linspace(zstart,zend,Nz+1))[2:end], xe
end


end
