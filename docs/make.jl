using Documenter, Powerlog

push!(LOAD_PATH,"../src/")

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/mwsohn/Powerlog.jl.git",
    julia = "0.6",
    osname = "linux"
)
