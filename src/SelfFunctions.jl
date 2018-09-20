
struct SelfFunction{F<:Function}
    name::Symbol
    typ::Type
    f::F
end
Base.show(io::IO, sf::SelfFunction) = print(io, "$(sf.name) (self function of type $(sf.typ))")
@inline (sf::SelfFunction)(args...; kwargs...) = sf.f(args...; kwargs...)
@inline selfcall(f::SelfFunction, t, args...; kwargs...) = f.f(t, args...; kwargs...)
@inline selfcall(f, t, args...; kwargs...) = f(args...; kwargs...)


macro self(typ, funcdef)
    fields = @eval fieldnames($typ)
    sfuncdef = splitdef(funcdef)
    
    insert!(sfuncdef[:args],1,:(self::$typ))
    
    function visit(ex; inside_func_args=false, localvars=[])
        rvisit(ex; kwargs...) = visit(ex; localvars=localvars, kwargs...)
        if ex isa Symbol
            # replace `x` with `self.x` where needed 
            ex in fields && !(ex in localvars) && !inside_func_args ? :(self.$ex) : ex
        elseif ex isa Expr
            if isexpr(ex,:kw)
                # in f(x=x) only the second `x` should (possibly) become self.x
                Expr(:kw, ex.args[1], rvisit(ex.args[2]))
            elseif @capture(ex, ((f_(args__; kwargs__)) | (f_(args__; kwargs__)::T_) | (f_(args__)) | (f_(args__)::T_)))
                # pass `self` argument implicitly when calling other "self" functions
                ex = :(selfcall($(rvisit(f)), self, $(rvisit.(args == nothing ? [] : args)...);  $(rvisit.(kwargs == nothing ? [] : kwargs)...)))
                T == nothing ? ex : :($ex::$T)
            elseif isdef(ex)
                # inner function definition, need to be careful about scope here
                sdef = splitdef(ex)
                map!(x->rvisit(x; inside_func_args=true), sdef[:args], sdef[:args])
                func_args = map(first,map(splitarg,sdef[:args]))
                sdef[:body] = visit(sdef[:body], localvars=[localvars; func_args])
                combinedef(sdef)
            else
                # recurse
                Expr(ex.head, map(rvisit, ex.args)...)
            end
        else
            ex
        end
    end

    sfuncdef[:body] = visit(sfuncdef[:body])
    fname = sfuncdef[:name]
    sfuncdef[:name] = Symbol("self_", sfuncdef[:name])

    quote
        $(esc(combinedef(sfuncdef)))
        Base.@__doc__ const $(esc(fname)) = SelfFunction($(QuoteNode(fname)), $(esc(typ)), $(esc(sfuncdef[:name])))
    end
    
end
