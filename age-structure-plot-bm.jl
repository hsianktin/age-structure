using PGFPlotsX
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{colormaps}")
t = @pgf Table({x = "μ_0", y = "β_0", z = "N", "col sep" = "comma"}, "$(label).csv")
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "death rate \$\\mu\$",
        ylabel = "birth rate \$\\beta\$",
        title = "population at equilibrium \$ N \$",
        legend_pos="south east",
        view = (0, 90),
        colorbar,
        "colormap/temp",
    },
    Plot3(
        {
            surf,
            shader = "flat",
            "mesh/rows" = length(βs),
            "mesh/cols" = length(μs),
        },
        t
    ),
    LegendEntry("$type")
)
pgfsave("fig/heatmap-bm-$(label).tex", axis)

# identify nulcline
N̄ = maximum(df.N)
β⁺s = Float64[]
μ⁺s = Float64[]
for β in βs
    temp_df = df[df.β_0 .== β,:]
    for i in 2:length(temp_df.μ_0)-1
        if temp_df.N[i+1] < temp_df.N[i] && (temp_df.N[i] > temp_df.N[i-1]) 
            push!(μ⁺s,temp_df.μ_0[i])
            push!(β⁺s,β)
        end
    end
end
nullcline = DataFrame(
    μ = μ⁺s,
    β = β⁺s,
    N = [N̄ for i in 1:length(μ⁺s)],
)
CSV.write("fig/nullcline-$(label).csv", nullcline)

t = @pgf Table({x = "μ", y = "β", z = "N", "col sep" = "comma"}, "nullcline-$(label).csv")
plot_line = @pgf Plot3(
    {
        black,
        thick,
        dashed,
    },
    t
)
legend_entry = LegendEntry("nullcline")
push!(axis,plot_line)
push!(axis,legend_entry)
pgfsave("fig/nullcline-$(label).tex", axis)

# t = @pgf Table({x = "μ_0", y = "β_0", z = "∂ᵤN", "col sep" = "comma"}, "$(label).csv")
# @pgf axis = Axis(
#     {
#         width = "3.4in",
#         height = "3.4in",
#         xlabel = "death rate \$\\mu\$",
#         ylabel = "birth rate \$\\beta\$",
#         title = "\$ \\partial_{\\mu}N \$",
#         view = (0, 90),
#         colorbar,
#         "colormap/jet",
#     },
#     Plot3(
#         {
#             surf,
#             shader = "flat",
#             "mesh/rows" = length(βs),
#             "mesh/cols" = length(μs),
#         },
#         t
#     )
# )
# pgfsave("fig/diff-heatmap-bm-$(label).tex", axis)

# t = @pgf Table({x = "μ_0", y = "β_0", z = "∂ᵤlnN", "col sep" = "comma"}, "$(label).csv")
# @pgf axis = Axis(
#     {
#         width = "3.4in",
#         height = "3.4in",
#         xlabel = "death rate \$\\mu\$",
#         ylabel = "birth rate \$\\beta\$",
#         title = "\$ \\partial_{\\mu} \\ln N \$",
#         view = (0, 90),
#         colorbar,
#         "colormap/jet",
#     },
#     Plot3(
#         {
#             surf,
#             shader = "flat",
#             "mesh/rows" = length(βs),
#             "mesh/cols" = length(μs),
#         },
#         t
#     )
# )
# pgfsave("fig/diff-normal-heatmap-bm-$(label).tex", axis)
