
macro eqs(ex)
    @assert isa(ex,Expr) && ex.head == :block
    tmps = []
    map!(ex.args, ex.args) do arg
        if isa(arg,Expr)
            if @capture(arg, lhs_ == rhs_)
                tmp = gensym()
                push!(tmps, tmp)
                :($tmp = $lhs - $rhs)
            else
                arg
            end
        else
            arg
        end
    end
    esc(quote
        $ex
        [$(tmps...)]
    end)
end

function integrate(𝕡::Params{T}, f, xmin, xmax) where {T}
    quadgk(f, convert(T,xmin), convert(T,xmax); rtol=𝕡.reltol)[1]::T
end

function find_zero(𝕡::Params{T}, f, xmin, xmax) where {T}
    Roots.find_zero(f, (convert(T,xmin), convert(T,xmax)), Roots.A42(); xrtol=𝕡.reltol)::T
end
