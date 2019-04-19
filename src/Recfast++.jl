@self Params{T} function init_reio!() where {T}

    Nz = 10000
    zstart = 1e4
    zend = 1e-3
    zarr, Xe_H, Xe_He, Xe, TM = (Array{Cdouble}(undef,Nz) for i=1:5);
    params = Cdouble[Nz, zstart, zend, Yp, Tcmb, Ωc+Ωb+Ων, Ωb, ΩΛ, Ωk, H0/100, Nν_massless+Nν_massive, 1.14, 0, 0]

    ccall(
        (:Xe_frac, "/home/marius/src/Recfast++/libRecfast++.so"), Cint,
        (Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint),
        params, zarr, Xe_H, Xe_He, Xe, TM, 1
    )
    
    itp = Spline1D(reverse(zarr), reverse(Xe))
    xe = z->itp(T(z))

end
