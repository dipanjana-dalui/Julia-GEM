"""

Here we provide a standalone validation for two major function,
    draw_new_trait and pick_indiv.
These are important functions as they handle the heritability and 
the coefficient of variation respectively.

Using a scratch example first: 
Instead of running the full Gillespie simulation loop, this builds a small
dummy population directly, calls draw_new_traits many times across a range
of parent trait values and heritability (h2) settings, logs every call to a
DataFrame, and plots parent trait vs. offspring trait.



What to expect if draw_new_traits is working correctly:
  - h2 = 0.0  -> offspring trait clusters around the population mean (red dotted line),
                 basically ignoring the parent's trait value
  - h2 = 1.0  -> offspring trait tracks the parent's trait almost exactly (black dashed y=x line)
  - h2 in between -> offspring trait is a blend of the two, with more scatter


What to expect if pick_indiv is working correctly:
  - as the CV increases, the distrubution of the parameters should be wider.

"""

using DataFrames, Statistics, Random, Distributions, Plots

include("src/functions/DrawNewTrait.jl")
  

# ============================================================== #
#  Build a small example population for any one state
# ============================================================== #
Random.seed!(42)

no_param   = 1
no_columns = 3
state      = 1  


pop_traits = 0.3 .+ 0.001 .* randn(100)          # 100 fake individuals
genotypes  = rand([0.0, 1.0], 100)
x_dist     = hcat(fill(Float64(state),100), pop_traits, genotypes)

pop_mean_true = mean(pop_traits)

println("dummy population trait mean: ", round(pop_mean_true, digits=2))

 
 
h2_values            = 0.0:0.25:1.0
parent_trait_values  = 0.2 .+ 0.15 .* rand(100)     # range of possible dummy parent phenotypes values
n_reps               = 30           # repeat draws at each point to see stochastic spread

results = DataFrame(h2=Float64[], parent_trait=Float64[],
                    pop_mean=Float64[], offspring_trait=Float64[])

for h2_val in h2_values
    h2 = fill(h2_val, 3, 1)   # dummy h2 matrix; only h2[state, 1] is read by DrawNewTraits
    for pt in parent_trait_values
        parent_traits = [pt, 1.0]   # [quant trait, genotype] -- genotype just passes through
        for _ in 1:n_reps
            out = draw_new_trait(x_dist, parent_traits, h2, no_param, no_columns, state, 1)
            push!(results, (h2_val, pt, pop_mean_true, out.offspring_traits[1]))
        end
    end
end

results

plt_h = Plots.plot(layout=(1, length(h2_values)), size=(1500, 500), legend=false)

for (i, h2_val) in enumerate(h2_values)
    sub = results[results.h2 .== h2_val, :]
    Plots.scatter!(plt_h[i], sub.parent_trait, sub.offspring_trait,
             xlabel = "parent trait",
             ylabel = i == 1 ? "offspring trait" : "",
             title = "h2 = $h2_val", markersize=2, alpha=0.35, color=:steelblue,
             label = i == 1 ? "offspring data" : "",
             legend = i == 1 ? :topleft : false,
             left_margin = i == 1 ? 8Plots.mm : 2Plots.mm,
             bottom_margin = 6Plots.mm)
    Plots.plot!(plt_h[i], sub.parent_trait, sub.parent_trait,
          color=:black, linestyle=:dash,
          label = i == 1 ? "y = x" : "")
    hline!(plt_h[i], [pop_mean_true], color=:red, linestyle=:dot,
           label = i == 1 ? "population mean" : "")
end

display(plt_h)


Plots.savefig(plt_h, "drawnewtrait_valid.svg")

#========================================================================#
# 			Pick individual validation
#========================================================================#
include("src/functions/PickIndiv.jl")  

Random.seed!(42)

x_mean = 0.30    # target mean
x_std  = 0.03    # target std deviation
N      = 100000      # number of samples


cv_values = [ 0.1, 0.25, 0.5, 0.75, 1.0]
colors = [:steelblue, :darkorange, :seagreen, :purple, :crimson]

plt2 = Plots.plot(size=(900, 600), xlabel="trait value", ylabel="count",
                  title="pick_indiv samples across CV values", legend=:topright)

for (i, cv) in enumerate(cv_values)
    new_x_std = x_std * cv
    x_out = pick_indiv(x_mean, new_x_std, N)

    Plots.histogram!(plt2, x_out, bins=50, alpha=0.2, color=colors[i],
               label="cv = $cv")
end

Plots.vline!(plt2, [x_mean], color=:black, linestyle=:dash, linewidth=2, label="target mean")

display(plt2)

Plots.savefig(plt2, "pickindiv_valid.svg")
