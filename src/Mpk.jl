module Mpk

using With
using Cosmology: Hubble, Params, ργ, ρν
using Sundials

const i = 0+im


function odesolve(F,y0,dy0,ts)

    n = length(y0)
    y0 = [real(y0); imag(y0)]
    dy0 = [real(dy0); imag(dy0)]
    
    function Fwrap(t,y,dy,eqs)
        y′  =  y[1:n] + i *y[n+1:end]
        dy′ = dy[1:n] + i*dy[n+1:end]
        
        eqs′ = F(t,y′,dy′)
    
        eqs[1:n]     = real(eqs′)
        eqs[n+1:end] = imag(eqs′)
    end
        
    Sundials.idasol(Fwrap,y0,dy0,ts)
end

macro eqs(ex)
    @assert isa(ex,Expr) && ex.head == :block
    
    rex = :([])
    
    for (i,arg) in enumerate(ex.args)
        if ~isa(arg,LineNumberNode)
            println(isa(arg,LineNumberNode),arg.head,"\n")
            @assert (isa(arg,Expr) && (arg.head == :comparison))
            push!(rex.args,:($(arg[1]) - $arg[3]))
        end
    end
    
    :($(esc(rex)))
end

@with Params function mpk(k,z)
    
    function F(η, y, dy)

        Θ₀, Θ₁, δ, v, Φ, a = y
        dΘ₀dη, dΘ₁dη, dδdη, dvdη, dΦdη, dadη = dy
        
        z = real(1/a-1)
        
        #see e.g. Dodelson, Modern Cosmology, pg 186
        @eqs begin
            dΘ₀dη + k*Θ₁ == -dΦdη
            dΘ₁dη - k*Θ₀/3 == -k*Φ/3
            dδdη + i*k*v == -3dΦdη
            dvdη + (dadη/a)*v == i*k*Φ
            (k^2)*Φ + 3(dadη/a)*(dΦdη+(dadη/a)*dΦdη) == 4*π*a^2*(δ*ρc₀*a^-3 + (ργ(z)+ρν(z))*Θ₀)
            (dadη/a^2) == Hubble(z)
        end
        
    end
    
    a₀ = 1e-6
    z₀ = 1/a₀-1
    η₀ = 
    
    y₀ =  [1, 0, 1,   0, 1, a₀]
    dy₀ = [0, 0, 0,   0, 0, a₀^2*Hubble(1/a₀-1)]
    
    odesolve(F,y₀,dy₀,[z₀,z])
    


end


end
