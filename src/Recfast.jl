module Recfast

function xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
            HOinp::Float64, Tnow::Float64, Yp::Float64, 
            Hswitch::Int64, Heswitch::Int64, 
            Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)
            
    xe = Array{Float64}(Nz)
    ccall( 
        (:get_xe_, "recfast"), Void, 
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
         Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        OmegaB, OmegaC, OmegaL, HOinp, Tnow, Yp, Hswitch, Heswitch, Nz, zstart, zend, xe
    )
    xe
end

push!(Base.Libdl.DL_LOAD_PATH,".")
for i=1:10
    xe(0.05,0.25,0.7,70.,2.72,0.25,0,1)
end

end
