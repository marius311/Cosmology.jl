module Mpk

using With
using Cosmology: Hubble, Params, ργ, ρν, η, ρc, ρb
using ODESolve: @eqs, odesolve
using PhysicalConstants: mp, σT

export mpk

@with Params function mpk(k,zs::Array)

    function F(lna, y, dy_dlna)

        a = exp(lna)
        z = 1/a - 1
        H = Hubble(z)
        k² = k^2
        
        δγ,  θγ,  δc,  θc,  δb,  θb,  ϕ  = y
        δγ′, θγ′, δc′, θc′, δb′, θb′, ϕ′ = dy_dlna*a*H
        a′ = a^2 * H
        
        ψ = ϕ
        
        R = (3ρb(z))/(4ργ(z))
        cₛ² = 1/(3(1+R))
        nₑ = ρb(z)/mp
        
        @eqs begin
            #Photons 
            δγ′ == -4/3*θγ + 4ϕ′
            θγ′ == k²/4*δγ + k²*ψ
            
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
    
    #             δγ  θγ  δc  θc  δb  θb  ϕ
    y₀ =  Float64[4,  0,  3,  0,  3,  0, -2]
    y′₀ = Float64[0,  0,  0,  0,  0,  0,  0]
    
    a₀ = 1e-7
    as = 1./(1+zs)
    
    # println(F(a₀,y₀,y′₀))
    
    y, dy = odesolve(F,y₀,y′₀,log([a₀,as...]))
    
end


@with Params mpk(k,z::Float64) = mpk(k,[z])


end
