using Documenter
using EDO

makedocs(
    modules = [EDO],
    sitename = "EDO.jl",
    authors = "Saloua Naama, Mohamed El Waghf et Rachid ELMontassir",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
            "Accueil" => "index.md",
            "Sujet TP ordre "=>"sujet_ordre.md",
            "Résultats TP Ordre "=>"resultats_sujet_ordre.md"
            ]
    )

deploydocs(repo = "github.com/mathn7/EDO.git")
