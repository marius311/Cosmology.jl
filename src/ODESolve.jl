"""

Provides a DAE (differential algebraic equation) solver, optionally with complex
variables, based on Sundials' IDASOL, along with a macro `@eqs` to make writing
the target function easier to read. 

For example, to solve the system of equations:

z'(t) + x(t) = 0
z(t) = x'(t)

we could write the target function as,

```julia
function F(t,y,dydt)
    x,z = y
    dxdt,dydt = dydt
    @eqs begin
        dzdt + x == 0
        z == dxdt
    end
end
```

If we impose x(0) = z'(0) = 0 and z(0) = x'(0) = 1, such that the solutions are
x(t) = sin(t) and z(t) = cos(t), we can then invoke the solver by:

```julia

y0 = [0,1]
dydt0 = [1,0]
ts = linspace(0,4π) #solve for two periods
odesolve(F,y0,dydt0,ts)
```

"""
module ODESolve

using Sundials

function odesolve(F,y0::Array{Complex128},dy0::Array{Complex128},ts; kwargs...)

    n = length(y0)
    y0 = [real(y0); imag(y0)]
    dy0 = [real(dy0); imag(dy0)]
    
    function Fwrap(t,y,dy,eqs)
        y  =  y[1:n] + im *y[n+1:end]
        dy = dy[1:n] + im*dy[n+1:end]
        
        _eqs = F(t,y,dy)
    
        eqs[1:n]     = real(_eqs)
        eqs[n+1:end] = imag(_eqs)
    end
        
    y,dy = Sundials.idasol(Fwrap,y0,dy0,ts; kwargs...)
    (y[:,1:n] + im*y[:,n+1:end]), (dy[:,1:n] + im*dy[:,n+1:end])

end


function odesolve(F,y0::Array{Float64},dy0::Array{Float64},ts; kwargs...)
    
    function Fwrap(t,y,dy,eqs)
        eqs[1:end] = F(t,y,dy)
    end
        
    Sundials.idasol(Fwrap,y0,dy0,ts; kwargs...)

end


macro eqs(ex)
    @assert isa(ex,Expr) && ex.head == :block
    
    rex = :([])
    
    for (i,arg) in enumerate(ex.args)
        @assert isa(arg,Expr)
        if (arg.head != :line)
            @assert (arg.head == :comparison) && (arg.args[2] == :(==))
            push!(rex.args,:($(arg.args[1]) - $(arg.args[3])))
        end
    end
    
    :($(esc(rex)))
end




end
