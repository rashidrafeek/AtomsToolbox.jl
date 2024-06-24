using Documenter, AtomsToolbox, AtomsBase

makedocs(;
    sitename="AtomsToolbox.jl",
    # modules = [AtomsToolbox]
)

deploydocs(
    repo = "github.com/rashidrafeek/AtomsToolbox.jl.git",
)
