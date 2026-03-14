# Author: Jordi Torres Vallverdú
# Start Date: 12/03/2025

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


######## Baseline ########

#static_profit(c)
#eql_ma(c)
#ds_ma(c, "baseline.tex")


######## Low Entry Cost ########

c.entry_low = 0.01
c.entry_high = 0.11

#static_profit(c)
#eql_ma(c)
#ds_ma(c, "low_entry_cost.tex")

####FUNCTION 1: STATIC PROFITS


#This generates the states and the index. Probably is not so efficient as the matlab computation, but it is much more simple to compute and to understand.
function states(kmax)
    states = Dict()
    id = 1
    for w1 in 0:kmax
        for w2 in 0:w1
            for w3 in 0:w2
                states[(w1, w2, w3)] = id
                id += 1
            end
        end
    end
    return states
end

data_states = states(kmax)
println(data_states)

#Now we define the connection between the function of the states and the numeric index. It should be easy to check as it should be indexed easily. 
println(length(data_states))

#I guess this is a version of the decode/encode function, if needed afterwards.
index_to_state = Dict(v => k for (k, v) in data_states)
state_to_index = Dict(k => v for (k, v) in data_states)

println(index_to_state)

# Now we define cost function
function theta(gamma, w)
    cost = @. gamma * exp(-(w - 4.0)) #this should be in vectorized form for slight efficiency in terms of coding 
    return cost
end


#use an example. 
D = intercept
f = fixed_cost
gamma = 1
# Cournot equilibrium solver
function find_cournot_eq(cost, n_firms, p, D, f)

    n = n_firms

    while !((p - cost[n] >= 0) || (n == 1))
        n -= 1
        p = (D + sum(cost[1:n])) / (n + 1)
    end

    q = zeros(length(cost))   # quantities for all firms

    if p - cost[n] > 0
        q[1:n] .= p .- cost[1:n]
    end

    pstar = D - sum(q)

    profstar = (pstar .> cost) .* (pstar .- cost) .* q .- f

    return pstar, profstar, q
end


function ccprofit_all(D, f, index_to_state, gamma)
    S = length(index_to_state)
    mat_results = zeros(S, 3) # to store the results we are interested for each state. I will collect number of firms, profits, for state.

    for s in 1:S

        w1, w2, w3 = index_to_state[s]

        w = vcat(w1, w2, w3)              #vectorized form
        n_firms = sum(w .> 0)             #number of active firms, this method should be wlog. 

        cost = theta(gamma, w)

        if n_firms > 0
            # initial price guess using active firms only
            p = (D + sum(cost[1:n_firms])) / (n_firms + 1)

            pstar, profstar, q = find_cournot_eq(cost, n_firms, p, D, f)

            mat_results[s, :] .= profstar
        end


    end
    return mat_results
end


mat_results = ccprofit_all(D, f, index_to_state, gamma)

#This works quite well. Results also make sense. Find a way to represent them reasonbly well. 

#######################################################################################
############################# PART 2. 

#This starts from the matlab code eql_ma.m

rlnfirm = max_firms
#kmax, already defined 
stfirm = start_firms
x_entryl = entry_low
x_entryh = entry_high
phi = scrap_val
entry_k = entry_at
#beta, already defined. Note: this I should improve in the next iteration of the code, to make it robust to different versions. 
#delta 
a = inv_mult
#tol 

#Let's start with the main function of cal val. This is all testing so far.
function calval(nfirms)
    z1 = zeros(nfirms, 1) #lower bound states
    z2 = kmax * ones(nfirms, 1) #upper bound states.

    if nfirms>1
        if place ==1
            locmask= zeros


end








