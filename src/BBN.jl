
function init_bbn!(𝕡)
    dat = readdlm(joinpath(dirname(@__FILE__),"..","dat","BBN_full_alterBBN_880.1.dat"),comments=true)
    ωbs, ΔNeffs = unique(dat[:,1]), unique(dat[:,3])
    itp = Spline2D(ωbs, ΔNeffs, reshape(dat[:,5],(length(ωbs),length(ΔNeffs))))
    @set! 𝕡.Yp = itp(𝕡.ωb, 𝕡.Nν_massive+𝕡.Nν_massless-3.046)
end
