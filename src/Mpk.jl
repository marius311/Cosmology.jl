module Mpk

using With
using Cosmology: Hubble, Params, ργ, ρν, η
using ODESolve: @eqs, odesolve

const i = im

@with Params function mpk(k,zs::Array)
    
    function F(η, y, dy)

        Θ₀, Θ₁, δ, v, Φ, a = y
        dΘ₀dη, dΘ₁dη, dδdη, dvdη, dΦdη, dadη = dy
        
        local z = real(1/a-1)
        
        #see e.g. Dodelson, Modern Cosmology, pg 186
        H = Hubble(z)
        @eqs begin
            dΘ₀dη + k*Θ₁ == -dΦdη
            dΘ₁dη - k*Θ₀/3 == -k*Φ/3
            dδdη + i*k*v == -3dΦdη
            dvdη + (dadη/a)*v == i*k*Φ
            (k^2)*Φ + 3*a*H*(dΦdη+a*H*Φ) == 4*π*a^2*(δ*ρc₀*a^-3 + (ργ(z)+ρν(z))*Θ₀)
            (dadη/a^2) == H
        end
        
    end
    
    a₀ = 1e-5
    z₀ = 1/a₀-1
    η₀ = η(z₀)
    H = Hubble(z₀)
    
    #                Θ₀, Θ₁,   δ,   v,     Φ, a
    y₀ =  Complex128[1,  0,    1/3, 0,     2, a₀]
    dy₀ = Complex128[0,  -k/3, 0,   2*i*k, -3/2*a₀*H, a₀^2*H]
    
    ηs = [η(z) for z in zs]
    y, dy = odesolve(F,y₀,dy₀,[η₀,ηs...],diffstates=fill(true,6))
    
end


end
