module Mpk

using With
using Cosmology: Hubble, Params, ργ, ρν, η
using ODESolve: @eqs, odesolve

export mpk

const i = im

@with Params function mpk(k,zs::Array)

    function F(lna, y, dy_dlna)

        a = exp(lna)
        z = 1/a - 1
        H = Hubble(z)
        
        Θ₀, Θ₁, δ, v, Φ = y
        dΘ₀dη, dΘ₁dη, dδdη, dvdη, dΦdη = dy_dlna*a*H
        dadη = a^2 * H
        
        #see e.g. Dodelson, Modern Cosmology, pg 186
        @eqs begin
            dΘ₀dη + k*Θ₁ == -dΦdη
            dΘ₁dη - k*Θ₀/3 == -k*Φ/3
            dδdη + i*k*v == -3dΦdη
            dvdη + (dadη/a)*v == i*k*Φ
            (k^2)*Φ + 3*(dadη/a)*(dΦdη+(dadη/a)*Φ) == 4*π*a^2*(δ*ρc₀*a^-3 + 4(ργ(z)+ρν(z))*Θ₀)
        end
        
    end
    
    a₀ = 1e-7
    z₀ = 1/a₀-1
    η₀ = η(z₀)
    H = Hubble(z₀)

    #                Θ₀, Θ₁, δ, v, Φ
    y₀ =  Complex128[1,  0,  3, 0, 2]
    dy₀ = Complex128[0,  0,  0, 0, 0]
    
    as = 1./(1+zs)
    y, dy = odesolve(F,y₀,dy₀,log([a₀,as...]))
    
end


end
