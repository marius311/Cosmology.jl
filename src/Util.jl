
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
