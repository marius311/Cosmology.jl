using BinDeps

@BinDeps.setup

librecfast = library_dependency("librecfast")

provides(SimpleBuild,
    (@build_steps begin
        ChangeDirectory(joinpath(BinDeps.depsdir(librecfast),"src","recfast"))
        `make`
    end), librecfast, os = :Unix)
    
@BinDeps.install Dict([(:librecfast, :librecfast)])
