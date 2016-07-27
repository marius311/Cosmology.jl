module TypeDefaults

export @with_defaults

# Transforms :(a::b) -> :a, from 
decolon2(a::Expr) = (@assert a.head==:(::);  a.args[1])
decolon2(a::Symbol) = a

# Returns the name of the type as Symbol
function typename(typedef::Expr)
    if isa(typedef.args[2], Symbol)
        return typedef.args[2]
    elseif isa(typedef.args[2].args[1], Symbol)
        return typedef.args[2].args[1]
    elseif isa(typedef.args[2].args[1].args[1], Symbol)
        return typedef.args[2].args[1].args[1]
    else
        error("Could not parse type-head from: $typedef")
    end
end


macro with_defaults(typ)
    
    defaults = Dict()
    
    # go through checking for fields with a default value, 
    # strip the defaults from the definition and store it in `defaults`
    args = typ.args[3].args
    for i in 1:length(args)
        ex = args[i]
        if isa(ex,Expr) && (ex.head == :(=)) && !(isa(ex.args[1], Expr) && ex.args[1].head==:call)
            args[i] = ex.args[1]
            defaults[decolon2(ex.args[1])] = ex.args[2]
        end
    end
    
    # add our own constructor which accepts keywords and 
    # sets the defaults we collected earlier
    constructor = quote
        function $(typename(typ))(;kwargs...)
            self = new()
            kwargs = Dict(kwargs)
            for (k,v) in kwargs; self.(k) = v; end
            for (k,v) in $defaults; if !haskey(kwargs,k); (self.(k) = v); end; end
            self
        end
    end
    push!(typ.args[3].args, constructor)
        
    :($(esc(typ)))
end


# function show(io::IO, p::Params)
#     show(io,[k=>p.(k) for k=fieldnames(p)])
# end


end
