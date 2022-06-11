using PGFPlotsX

t = @pgf Table({x = "μ_0", y = "β_0", z = "N", "col sep" = "comma"}, "$(label).csv")
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "death rate \$\\mu\$",
        ylabel = "birth rate \$\\beta\$",
        title = "population at equilibrium \$ N \$",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
            "mesh/rows" = length(βs),
            "mesh/cols" = length(μs),
        },
        t
    )
)
pgfsave("fig/heatmap-bm-$(label).tex", axis)

t = @pgf Table({x = "μ_0", y = "β_0", z = "∂ᵤN", "col sep" = "comma"}, "$(label).csv")
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "death rate \$\\mu\$",
        ylabel = "birth rate \$\\beta\$",
        title = "\$ \\partial_{\\mu}N \$",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
            "mesh/rows" = length(βs),
            "mesh/cols" = length(μs),
        },
        t
    )
)
pgfsave("fig/diff-heatmap-bm-$(label).tex", axis)

t = @pgf Table({x = "μ_0", y = "β_0", z = "∂ᵤlnN", "col sep" = "comma"}, "$(label).csv")
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "death rate \$\\mu\$",
        ylabel = "birth rate \$\\beta\$",
        title = "\$ \\partial_{\\mu} \ln N \$",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
            "mesh/rows" = length(βs),
            "mesh/cols" = length(μs),
        },
        t
    )
)
pgfsave("fig/diff-normal-heatmap-bm-$(label).tex", axis)
