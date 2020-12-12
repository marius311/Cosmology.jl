module Mpk

using Cosmology: Hubble, Params, ργ, ρν, η, ρc, ρb
using ODESolve: @eqs, odesolve
using PhysicalConstants: mp, σT
using NamedArrays

export solve_boltz

@self Params function solve_boltz(ks::Array, zs::Array)

    function F(k, lna, y, dy_dlna)

        a = exp(lna)
        z = 1/a - 1
        H = Hubble(z)
        k² = k^2
        
        δγ,  θγ,  δc,  θc,  δb,  θb,  ϕ  = y
        δγ′, θγ′, δc′, θc′, δb′, θb′, ϕ′ = dy_dlna*a*H
        a′ = a^2 * H
        
        ψ = ϕ
        
        R = 3ρb(z)/4ργ(z)
        cₛ² = 1/(3(1+R))
        nₑ = ρb(z)/mp
        
        @eqs begin
            #Photons 
            δγ′ == -4/3*θγ + 4ϕ′
            θγ′ == k²/4*δγ + k²*ψ + a*nₑ*σT*(θb-θγ)
            
            #CDM
            δc′ == -θc + 3ϕ′
            θc′ == -(a′/a)*θc + k²*ψ
            
            #Baryons
            δb′ == -θb + 3ϕ′
            θb′ == -(a′/a)*θb + cₛ²*k²*δb + a*nₑ*σT/R*(θγ-θb) + k²*ψ
            
            #Einstein equations
            k²*ϕ + 3*(a′/a)*(ϕ′+(a′/a)*ψ) == -4π*a^2*(ρc(z)*δc + ρb(z)*δb + (ργ(z)+ρν(z))*δγ)
        end
        
    end
    
    vars =        [ :δγ, :θγ, :δc, :θc, :δb, :θb, :ϕ ]
    y₀   = Float64[  4,   0,   3,   0,   3,   0,  -2 ]
    y′₀  = Float64[  0,   0,   0,   0,   0,   0,   0 ]
    
    a₀ = 1e-7
    as = 1./(1+zs)
    
    vars′ = [symbol(string(v)"′") for v in vars]
    soln = NamedArray(zeros(Float64,2*length(vars),length(ks),length(zs)), ([vars; vars′],ks,zs), ("var","k","z"))
    
    for (i,k) in enumerate(ks)
        y, y′ = odesolve((args...)->F(k,args...),y₀,y′₀,log([a₀,as...]))
        soln[1:length(vars),    i,:] = y[2:end,:]'
        soln[length(vars)+1:end,i,:] = y′[2:end,:]'
    end
    
    soln
    
end


end
