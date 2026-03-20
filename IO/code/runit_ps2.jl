#Author: Jordi Torres Vallverdú
# Start Date: 12/03/2025

using JLD2
using Statistics
using Random

############ Main Function ############

mutable struct ModelParams
    max_firms::Int
    kmax::Int
    start_firms::Int

    entry_low::Float64
    entry_high::Float64
    scrap_val::Float64
    entry_at::Int
    beta::Float64
    delta::Float64

    inv_mult::Float64

    intercept::Float64
    fixed_cost::Float64
    gamma::Float64

    tol::Float64

    profit_done::Bool
    eql_done::Bool
    prefix::String

    ds_wstart::Vector{Float64}
    ds_nsimx::Int
end

#### Model Parameters ####

max_firms = 3
kmax = 19
start_firms = 1

entry_low = 0.15
entry_high = 0.25
scrap_val = 0.1
entry_at = 4
beta = 0.925
delta = 0.7

inv_mult = 3

intercept = 3
fixed_cost = 0.2
gamma = 1

tol = 0.1

profit_done = false
eql_done = false
prefix = "cc"

ds_wstart = vcat(entry_at + 2, zeros(max_firms - 1))
ds_nsimx = 10000

c = ModelParams(
    max_firms,
    kmax,
    start_firms,
    entry_low,
    entry_high,
    scrap_val,
    entry_at,
    beta,
    delta,
    inv_mult,
    intercept,
    fixed_cost,
    gamma,
    tol,
    profit_done,
    eql_done,
    prefix,
    ds_wstart,
    ds_nsimx
)

############ RUN FUNCTIONS IN ORDER ############
include("states_ps2.jl") #this just runs my version of decode encode
include("static_profit_ps2.jl") #static profits
include("eqlma_ps2.jl") # policy function and fixed point
include("ds_ma_ps2.jl") # simulation code





############ RUN SIMULATION ############


Random.seed!(123) #so that results don't change in export... but not needed really. Just annoying to rewrite the tex... 
eql_ma(c)
avg_firms_base, avg_inv_base = ds_ma(c)

#Low entry cost
c.entry_low = 0.01
c.entry_high = 0.11

Random.seed!(123)
eql_ma(c)
avg_firms_low, avg_inv_low = ds_ma(c)

# Export LaTeX table
open("results.tex", "w") do f
    write(f, "\\begin{table}[h]\n")
    write(f, "\\centering\n")
    write(f, "\\begin{tabular}{lcc}\n")
    write(f, "\\hline\n")
    write(f, " & High Entry Cost & Low Entry Cost \\\\\n")
    write(f, " & [0.15, 0.25] & [0.01, 0.11] \\\\\n")
    write(f, "\\hline\n")
    write(f, "Average Active Firms & $(round(avg_firms_base, digits=4)) & $(round(avg_firms_low, digits=4)) \\\\\n")
    write(f, "Average Total Investment & $(round(avg_inv_base, digits=4)) & $(round(avg_inv_low, digits=4)) \\\\\n")
    write(f, "\\hline\n")
    write(f, "\\end{tabular}\n")
    write(f, "\\caption{Simulation Results: Baseline vs Low Entry Cost}\n")
    write(f, "\\label{tab:results}\n")
    write(f, "\\end{table}\n")
end

println("Results exported to results.tex")
println("EOF_____________________________:)")

###################EOF!