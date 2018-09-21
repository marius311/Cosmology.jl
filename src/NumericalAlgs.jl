
# adapted from https://github.com/scipy/scipy/blob/8dba340293fe20e62e173bdf2c10ae208286692f/scipy/optimize/Zeros/brentq.c
# Roots.jl really needs to have this... 

function brentq(f, xa::T, xb::T) where {T}
    xtol = 2e-12
    rtol = 4eps()
    maxiter = 100

    xpre, xcur = xa, xb
    xblk, fblk, spre, scur = zero(T), zero(T), zero(T), zero(T)

    fpre = f(xpre)
    fcur = f(xcur)
    if (fpre*fcur > 0)
        throw("not a bracketing interval")
    end
    if (fpre == 0)
        return xpre
    end
    if (fcur == 0)
        return xcur
    end

    for i=0:maxiter
        if (fpre*fcur < 0)
            xblk = xpre
            fblk = fpre
            spre = scur = xcur - xpre
        end
        if (abs(fblk) < abs(fcur))
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre
        end

        delta = (xtol + rtol*abs(xcur))/2
        sbis = (xblk - xcur)/2
        if (fcur == 0 || abs(sbis) < delta)
            return xcur
        end

        if (abs(spre) > delta && abs(fcur) < abs(fpre))
            if (xpre == xblk)
                # /* interpolate */
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                # /* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre))
            end
            if (2*abs(stry) < min(abs(spre), 3*abs(sbis) - delta))
                # /* good short step */
                spre = scur
                scur = stry
            else
                # /* bisect */
                spre = sbis
                scur = sbis
            end
        else
            # /* bisect */
            spre = sbis
            scur = sbis
        end

        xpre,fpre = xcur, fcur
        if (abs(scur) > delta)
            xcur += scur
        else
            xcur += (sbis > 0 ? delta : -delta)
        end

        fcur = f(xcur)
    end
    return xcur
end
