module BBN

using Cosmology, SelfFunctions, Dierckx

export init_bbn!

@self Params function init_bbn!()
    dat = readdlm(joinpath(dirname(@__FILE__),"..","dat","BBN_full_alterBBN_880.1.dat"))
    ωbs, ΔNeffs = unique(dat[:,1]), unique(dat[:,3])
    itp = Spline2D(ωbs, ΔNeffs, reshape(dat[:,5],(length(ωbs),length(ΔNeffs))))
    Yp = itp(ωb,Nν_massive+Nν_massless-3.046)
end

end
