# ========================================== #
#		        PLOT FUNCtions               #
# ========================================== #

function Pop_Plot(pop_time::DataFrame, i::Int64)
    pop_time = pop_time[pop_time.state_ID .== i, :]
    rep_col_name = names(pop_time)[4:size(pop_time)[2]]
    pop_stack = stack(pop_time, 
                    rep_col_name)
    pop_plot = data(pop_stack) * mapping(
        :time,
        :value,
        #color = :replicate,
        group = :variable,
        row = :GEM_ver => nonnumeric 
        ) *
        visual(Lines, color = :grey)
    draw(pop_plot, axis=(xlabel="time", ylabel="population abundances"))    
end

# =====================================================

function Trait_Plot(mediandf::DataFrame, vardf::DataFrame, 
    spID::Int64, trait_to_plot::String )

    #pick the right GEM version number and state ID
    mediandf = mediandf[mediandf.state_ID .== spID, :]
    vardf = vardf[vardf.state_ID .== spID, :]

    #pick the trait wanted
    mediantemp = mediandf[:, [:time, :rep, :GEM_ver]]
    mediantemp2 = mediandf[:,trait_to_plot]
    median2plot = hcat(mediantemp, mediantemp2)

    #grab the right var
    vartemp2 = vardf[:,trait_to_plot]
    
    dftemp = hcat(median2plot, vartemp2, makeunique=true)
    rename!(dftemp, :x1 => :median, :x1_1 => :var)

    upper_bound = dftemp.median .+ dftemp.var
    lower_bound = dftemp.median .- dftemp.var
    df2plot = hcat(dftemp, upper_bound, lower_bound, makeunique = true)
    rename!(df2plot,:x1 => :upper_bound, :x1_1 => :lower_bound )

    median_plot = data(df2plot) * mapping(
        :time,
        :median,
        #color = :rep,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    ) *
        visual(Lines, color = :grey)

    var_plot = data(df2plot) * mapping(
        :time,
        :upper_bound,
        :lower_bound,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    )    * visual(Band, color = :lightgrey, alpha=0.5)

    tempplot = data(df2plot) * (median_plot + var_plot)
    finalplot = draw(tempplot, axis=(xlabel="time", ylabel="median trait value"))

    return finalplot
end

# =====================================================

function Geno_Freq_Plot(freqdf::DataFrame, 
                        spID::Int64, 
                        trait_to_plot::String )

    #pick the right GEM version number and state ID
    freqdf = freqdf[freqdf.state_ID .== spID, :]
   
    #pick the trait wanted
    freqtemp = freqdf[:, [:time, :rep, :GEM_ver]]
    freqtemp2 = freqdf[:,trait_to_plot]
    freq2plot = hcat(freqtemp, freqtemp2)

    
    rename!(freq2plot, :x1 => :freq)

    freq_plot = data(freq2plot) * mapping(
        :time,
        :freq,
        #color = :rep,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    ) * visual(Lines, color = :grey)
        

    finalplot = draw(freq_plot, axis=(xlabel="time", ylabel="genotype frequency"))

    return finalplot
end