module With

export register_with, @with

with_functions = Dict()
with_types = Dict()


#maybe this will get into Julia at some point? 
type_name(dt::DataType) = dt.name.name



function register_with(typ::Type, funcs=nothing)
    typname = type_name(typ)
    with_types[typname] = fieldnames(typ)
    if funcs != nothing
        with_functions[typname] = union(get(with_functions,typname,[]),funcs)
    end
    return
end


macro with(typ, func)
    
    if ! ((func.head == :function) || (func.head == :(=)))
        throw("@with used incorrectly. it should be '@with Type function ... end' or '@with Type f() = ...'")
    end
        
    #check we've registered this type
    if isa(typ,Expr) && (typ.head == :(::))
        self, typ = typ.args
    else
        self = gensym()
    end
    
    if typ in keys(with_types)  
        fields = with_types[typ]
    else
        throw(ArgumentError("Need to call register_with($typ) before using @with $typ."))
    end
    
    # add "self" as a first positional argument
    args = func.args[1].args
    i = ((length(args) >= 2) && isa(args[2],Expr) && (args[2].head == :parameters)) ? 3 : 2
    insert!(args,i,:($self::$typ))
        
    # recurse through the expression tree and rename X to self.X if X is a
    # fieldname of mytype (unless X refers to a functinon keyword arg, in which
    # case leave it)   also change f(...) to f(self,....) if f is a registered
    # "with" function
    function visit(ex)
        if typeof(ex) == Expr
            if (ex.head == :kw)
                ex.args[2:end] = map(visit,ex.args[2:end])
            else
                ex.args = map(visit,ex.args)
            end
            if (ex.head == :call) && (ex.args[1] in get(with_functions,typ,[]))
                insert!(ex.args,2,:($self))
            end            
        elseif (typeof(ex) == Symbol) & (ex in fields)
            return :($self.$ex)
        end
        ex
    end
    func.args[2] = visit(func.args[2])
        
    :($(esc(func)))
end



end
