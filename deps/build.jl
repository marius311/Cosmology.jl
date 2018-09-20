using BinDeps

@BinDeps.setup

librecfast = library_dependency("librecfast")

provides(SimpleBuild,
    (@build_steps begin
        ChangeDirectory(joinpath(BinDeps.depsdir(librecfast),"src","recfast"))
        `gfortran -ffixed-line-length-none -fPIC -shared -g -O0 recfast.for -o ../../usr/lib/librecfast.so`
    end), librecfast, os = :Unix)
    
@BinDeps.install Dict([(:librecfast, :librecfast)])
