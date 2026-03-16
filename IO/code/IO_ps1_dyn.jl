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


function state_helper(d, justone, z1, z2)

    #crucial part of the code, we need to reorder efficiencies after change in investment. 
    temp = hcat(d, justone)
    perm = sortperm(temp[:, 1], rev=true)
    temp = temp[perm, :]
    d_sorted = temp[:, 1] #key outcome

    e = d_sorted .- 1 #aggregate shock

    #check that is bounded
    e = max.(e, z1)
    d_sorted = min.(d_sorted, z2)

    #new place of the firm after the change of states
    pl1 = argmax(temp[:, 2])

    return e, d_sorted, pl1
end


#Let's start with the main function of cal val. This is all testing so far.
function calval(place, w, x, k, oldvalue, etable, multfac, two_n, kmax, nfirms, mask, delta, a)
    z1 = zeros(nfirms) #lower bound states
    z2 = fill(kmax, nfirms) #upper bound states.


    if nfirms > 1 ##should I condition on this right now?
        locmask = zeros(Int, nfirms, two_n)
        locmask[1:place-1, :] = mask[1:place-1, :]
        locmask[place+1:nfirms, :] = mask[place:nfirms-1, :] #Just inserts a row in mask where own firm is. This is much simpler than Matlab.
    else
        locmask = zeros(Int, 1, 1)
    end

    #modify investment and state. 
    x = copy(x)
    w = copy(w)

    x[place] = 0
    w[place] = k
    justone = zeros(nfirms)
    justone[place] = 1

    #probability of going up
    p_up = @. (a * x) / (1 + a * x)

    #initialize empty containers 
    val_up = 0.0
    val_stay = 0.0


    #initialize the loop to compute future values
    for i in 1:two_n #loop over the rivals?

        probmask = prod((locmask[:, i] .* p_up) .+ ((1 .- locmask[:, i]) .* (1 .- p_up)))
        ##############1. if firm does not move up. 

        d = w + locmask[:, i] #private shock

        e, d_sorted, pl1 = state_helper(d, justone, z1, z2)

        val_stay += (
            (1 - delta) * oldvalue[qencode(d_sorted, etable, multfac), pl1] +
            delta * oldvalue[qencode(e, etable, multfac), pl1]
        ) * probmask

        ##############2. if firm moves up 
        d_up = w + locmask[:, i]
        d_up[place] += 1 #just as if firm is successful too!

        e_up, d_up_sorted, pl1_up = state_helper(d_up, justone, z1, z2)

        val_up += (
            (1 - delta) * oldvalue[qencode(d_up_sorted, etable, multfac), pl1_up] +
            delta * oldvalue[qencode(e_up, etable, multfac), pl1_up]
        ) * probmask
    end

    return val_up, val_stay
end








