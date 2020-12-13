
function solve_boltz(𝕡::Params{T}, ks) where {T}

    function F(dy_dlna, y, k, lna)

        a = exp(lna)
        z = 1/a - 1
        H = Hubble(𝕡,z)
        k² = k^2
        
        δγ,  θγ,  δc,  θc,  δb,  θb,  ϕ  = y
        δγ′, θγ′, δc′, θc′, δb′, θb′, ϕ′ = dy_dlna*a*H
        a′ = a^2 * H
        
        ψ = ϕ
        
        R = 3ρb(𝕡,z)/4ργ(𝕡,z)
        cₛ² = 1/(3(1+R))
        nₑ = ρb(𝕡,z)/mp
        
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
            k²*ϕ + 3*(a′/a)*(ϕ′+(a′/a)*ψ) == -4π*a^2*(ρc(𝕡,z)*δc + ρb(𝕡,z)*δb + (ργ(𝕡,z)+ρν(𝕡,z))*δγ)
        
        end
        
    end
    
    vars =                             [ :δγ, :θγ, :δc, :θc, :δb,  :θb, :ϕ ]
    y₀   = ComponentArray((;(vars .=> T[   4,   0,   3,   0,   3,   0,  -2 ])...))
    y′₀  = ComponentArray((;(vars .=> T[   0,   0,   0,   0,   0,   0,   0 ])...))
    
    a₀ = 1e-7
    
    tmap(ks) do k
        solve(DAEProblem(F, y′₀, y₀, (log(a₀),0), k), IDA(), reltol=1e-4)
    end
    
end