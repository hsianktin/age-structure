using PGFPlotsX
function type2latex(type)
    if type == "semi-constant"
        return raw"$k=\theta(x'-a)\theta(a-x)$"
    elseif type == "semi-x′"
        return raw"$k=x'\theta(x'-a)\theta(a-x)$"
    elseif type == "semi-(x′-x)"
        return raw"$k=(x'-x)\theta(x'-a)\theta(a-x)$"
    else
        error("unknown type")
    end
end
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
        "colormap/hot2",
        "colorbar style"=@pgf {width="0.2cm"}
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
    LegendEntry(type2latex(type))
)

# identify nullcline on the left
N̄ = maximum(df.N)
β⁺s = Float64[]
μ⁺s = Float64[]
for β in βs
    temp_df = df[df.β_0 .== β,:]
    for i in 1:length(temp_df.μ_0)-1
        # if the population is increasing from i to i+1 but decreasing from i-1 to i or i == 1
        if temp_df.N[i+1] > temp_df.N[i] && (i == 1 || temp_df.N[i] < temp_df.N[i-1] )
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
sort!(nullcline,[:β], rev=true)
# identify nulcline on the right
N̄ = maximum(df.N)
β⁺s = Float64[]
μ⁺s = Float64[]
for β in βs
    temp_df = df[df.β_0 .== β,:]
    for i in 2:length(temp_df.μ_0)-1
        # if the population is decreasing from i to i+1 but increasing from i-1 to i
        if temp_df.N[i+1] < temp_df.N[i] && (temp_df.N[i] > temp_df.N[i-1]) 
            push!(μ⁺s,temp_df.μ_0[i])
            push!(β⁺s,β)
        end
    end
end
nullcline_r = DataFrame(
    μ = μ⁺s,
    β = β⁺s,
    N = [N̄ for i in 1:length(μ⁺s)],
)
nullcline = vcat(nullcline,nullcline_r)
CSV.write("fig/nullcline-$(label).csv", nullcline)

t = @pgf Table({x = "μ", y = "β", z = "N", "col sep" = "comma"}, "nullcline-$(label).csv")
plot_line = @pgf Plot3(
    {
        black,
        thick,
        dashed,
    },
    t,
    "\\closedcycle"
)
legend_entry = LegendEntry("nullcline")
push!(axis,plot_line)
push!(axis,legend_entry)
pgfsave("fig/nullcline-$(label).tex", axis)

# t = @pgf Table({x = "μ_0", y = "N", "col sep" = "comma"}, "$(label).csv")
# @pgf axis = Axis(
#     {
#         width = "3.4in",
#         height = "3.4in",
#         xlabel = "death rate \$\\mu\$",
#         ylabel = "population at equilibrium \$ N \$",
#         legend_pos="south east",
#     },
#     Plot(
#         t
#     ),
#     LegendEntry("$type")
# )
# pgfsave("fig/plot-$(label).tex", axis)

t = @pgf Table({x = "μ_0", y = "β_0", z = "Λ", "col sep" = "comma"}, "$(label).csv")
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "death rate \$\\mu\$",
        ylabel = "birth rate \$\\beta\$",
        title = "principle eigenvalue \$\\lambda_0\$",
        legend_pos="south east",
        view = (0, 90),
        colorbar,
        "colormap/hot2",
        "colorbar style"=@pgf {width="0.2cm"}
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
    LegendEntry(type2latex(type))
)
pgfsave("fig/eigenvals-$(label).tex", axis)
