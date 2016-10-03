module Recfast

function xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64, 
            HOinp::Float64, Tnow::Float64, Yp::Float64, 
            Hswitch::Int64, Heswitch::Int64)
            
    ccall( 
        (:__recfast_MOD_get, "recfast"), Void, 
        (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}),
        &OmegaB, &OmegaC, &OmegaL, &HOinp, &Tnow, &Yp, &Hswitch, &Heswitch
    )

end

push!(Base.Libdl.DL_LOAD_PATH,".")
xe(0.05,0.25,0.7,70.,2.72,0.25,0,0)

end
