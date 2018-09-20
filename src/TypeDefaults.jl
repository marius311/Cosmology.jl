
macro defaults(typ)
    
    @capture(typ, (mutable struct X_; body__; end) | (struct X_; body__; end)) || throw("Unrecognized type definition")
    
    # go through checking for fields with a default value, 
    # strip the defaults from the definition and store it in `defaults`
    defaults = Dict()
    map!(ex -> begin
        if @capture(ex, field_ = default_value_) && !isexpr(field,:call)
            @capture(field, (field_name_::T_) | (field_name_))
            defaults[field_name] = default_value
            field
        else
            ex
        end
    end, body, body)
    
    # add our own constructor which accepts keywords and 
    # sets the defaults we collected earlier
    @capture(X, ((T_{Vs__} <: S_) | T_{Vs__} | T_ ))
    type_args = map(ex->(@capture(ex, (V_ <: U_) | V_); V), Vs==nothing ? [] : Vs)
    constructor_name = Vs==nothing ? T : (:($T{$(type_args...)}))
    constructor = quote
        function $constructor_name(;kwargs...) where {$(type_args...)}
            self = new{$(type_args...)}()
            kwargs = Dict(kwargs)
            for (k,v) in kwargs
                setfield!(self,k,convert(typeof(getfield(self,k)),v))
            end
            for (k,v) in $defaults
                if !haskey(kwargs,k)
                    setfield!(self,k,convert(typeof(getfield(self,k)),v))
                end
            end
            self
        end
    end
    
    if @capture(typ, (struct _ __ end))
        esc(:(struct $X; $([body; constructor]...); end))
    else
        esc(:(mutable struct $X; $([body; constructor]...); end))
    end
end
