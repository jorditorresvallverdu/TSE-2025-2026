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
oldvalue = mat_results

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
function calcval(place, w, x, k,
    oldvalue,
    two_n, kmax, nfirms, mask,
    delta, a,
    state_to_index)

    z1 = zeros(Int, nfirms)
    z2 = fill(kmax, nfirms)

    # Build locmask
    if nfirms > 1
        locmask = zeros(Int, nfirms, two_n)
        locmask[1:place-1, :] = mask[1:place-1, :]
        locmask[place+1:nfirms, :] = mask[place:nfirms-1, :]
    else
        locmask = zeros(Int, 1, 1)
    end

    x = copy(x)
    w = copy(w)

    x[place] = 0
    w[place] = k

    justone = zeros(Int, nfirms)
    justone[place] = 1

    p_up = @. (a * x) / (1 + a * x)

    val_up = 0.0
    val_stay = 0.0

    for i in 1:two_n

        probmask = prod((locmask[:, i] .* p_up) .+ ((1 .- locmask[:, i]) .* (1 .- p_up)))

        # -------- STAY --------
        d = w + locmask[:, i]

        e, d_sorted, pl1 = state_helper(d, justone, z1, z2)

        idx_d = state_to_index[Tuple(sort(d_sorted, rev=true))]
        idx_e = state_to_index[Tuple(sort(e, rev=true))]

        val_stay += (
            (1 - delta) * oldvalue[idx_d, pl1] +
            delta * oldvalue[idx_e, pl1]
        ) * probmask

        # -------- MOVE UP --------
        d_up = w + locmask[:, i]
        d_up[place] += 1

        e_up, d_up_sorted, pl1_up = state_helper(d_up, justone, z1, z2)

        idx_d_up = state_to_index[Tuple(sort(d_up_sorted, rev=true))]
        idx_e_up = state_to_index[Tuple(sort(e_up, rev=true))]

        val_up += (
            (1 - delta) * oldvalue[idx_d_up, pl1_up] +
            delta * oldvalue[idx_e_up, pl1_up]
        ) * probmask
    end

    return val_up, val_stay
end


#######################################################################################
############################# PART 3. #################################################

#The idea here is to replicate the optimize function 
function unif_distribution_cdf(value, xh, xl)
    P = (value - xl) / (xh - xl)
    return P
end


function foc_investment(vj_up, vj_stay, a, beta)
    delta_cont = vj_up - vj_stay

    if delta_cont <= 0
        return 0.0 ##to makre sure we dont take sqrt of negative number
    end

    return (sqrt(a * beta * delta_cont) - 1) / a
end

function optimize(w, oldvalue, oldx, isentry, profit,
    wmax, two_n, kmax, nfirms, mask,
    phi, entry_k, beta, delta, a,
    state_to_index, index_to_state)

    # Decode state
    locw = collect(index_to_state[w])
    locwx = copy(locw)

    oval = copy(oldvalue[w, :])
    ox = copy(oldx[w, :])

    nval = zeros(nfirms)
    nx = zeros(nfirms)

    #################### 1. Exit condition ####################

    profit_w = profit[w, :]
    CV = profit_w .+ beta .* oval

    for s in 1:nfirms
        if CV[s] >= phi
            locwx[s] = locw[s]
        else
            locwx[s] = 0
        end
    end

    locwe = copy(locwx)

    #################### 2. Entry condition ####################

    p_entry = 0.0
    n_active = sum(locwe .> 0)

    if n_active == nfirms
        p_entry = 0.0
    else
        locwe[n_active+1] = entry_k
        locwe = sort(locwe, rev=true)
        idx_entry = state_to_index[Tuple(locwe)]
        p_entry = isentry[idx_entry]
    end

    #################### 3. Investment decision ####################

    for j in 1:nfirms
        if locwx[j] == 0
            nval[j] = phi
            nx[j] = 0.0
            continue
        end

        k = locwx[j]

        val_up_nentry, val_stay_nentry = calcval(
            j, locwx, ox, k,
            oldvalue,
            two_n, kmax, nfirms, mask,
            delta, a,
            state_to_index
        )

        val_up_entry, val_stay_entry = calcval(
            j, locwe, ox, k,
            oldvalue,
            two_n, kmax, nfirms, mask,
            delta, a,
            state_to_index
        )

        E_val_up = (1.0 - p_entry) * val_up_nentry + p_entry * val_up_entry
        E_val_stay = (1.0 - p_entry) * val_stay_nentry + p_entry * val_stay_entry

        xopt = foc_investment(E_val_up, E_val_stay, a, beta)
        xopt = max(0.0, xopt)

        p_up_j = (a * xopt) / (1 + a * xopt)

        profit_j = profit_w[j]
        Vj = profit_j + beta * (p_up_j * E_val_up + (1.0 - p_up_j) * E_val_stay) - xopt

        if Vj < phi
            locwx[j] = 0
            nval[j] = phi
            nx[j] = 0.0
        else
            nval[j] = Vj
            nx[j] = xopt
        end

        oval[j] = nval[j]
        ox[j] = nx[j]
    end

    return nval, nx
end


#######################################################################################
############################# PART 4. #################################################
# function contract 
function contract(oldvalue, oldx, profit,
    wmax, two_n, kmax, nfirms, mask,
    x_entryl, x_entryh, phi, entry_k, beta, delta, a,
    state_to_index, index_to_state)

    isentry = zeros(wmax)

    for w in 1:wmax
        locw = collect(index_to_state[w]) #this probably I need to adapt again. But will use an bottom up approach. 
        if (locw(nfirms) == 0)

            _, val_stay = calcval(place, w, x, k, oldvalue, two_n, kmax, nfirms, mask, delta, a, state_to_index)
            val = beta * val_stay
            isentry[w] = unif_distribution_cdf(val, x_entryh, x_entryl)

        end

    end

    #to bound it below 1 and above 0. 
    isentry = max.(isentry, 0.0)
    isentry = min.(isentry, 1.0)

    #Investment and value function. 
    newx = zeros(wmax, nfirms)
    newvalue = zeros(wmax, nfirms)

    for w in 1:wmax
        newx[w, :], newvalue[w, :] = optimize(w, oldvalue, oldx, isentry, profit,
            wmax, two_n, kmax, nfirms, mask,
            phi, entry_k, beta, delta, a,
            state_to_index, index_to_state)
    end

    return newvalue, newx, isentry

end

#######################################################################################
############################# PART 4. #################################################
function make_mask(nfirms)
    n = nfirms - 1
    two_n = 2^n

    if n == 0
        return zeros(Int, 1, 1)
    end

    return hcat([reverse(digits(i, base=2, pad=n)) for i in 0:two_n-1]...)
end

function eql_ma(c)

    #From params 

    rlnfirms = c.max_firms
    kmax = c.kmax
    stfirm = c.start_firms

    x_entryl = c.entry_low
    x_entryh = c.entry_high

    phi = c.scrap_val
    entry_k = c.entry_at
    beta = c.beta
    delta = c.delta
    a = c.inv_mult

    tol = c.tol

    #empty containers. 
    newvalue = zeros(wmax, nfirms)
    newx = zeros(wmax, nfirms)
    oldvalue = zeros(wmax, nfirms)
    oldx = zeros(wmax, nfirms)
    isentry = zeros(wmax)

    for nfirms in stfirm:rlnfirms

        #Number of combinations
        wmax = length(state_to_index) #this differs a bit. 

        #Load profit, should I call the function?
        #TO do

        #Number of rival actions 
        two_n = 2.0^(nfirms - 1)

        #Matrix with all possible actions of the other firms. 
        mask = make_mask(nfirms) #described above, same as dec2bin. 

        #dtable is what I already have, no? The rest I don't need at all. 

        oldvalue, oldx = initialize(dtable, nfirms, wmax, binom, newvalue, newx) #need to adapt initialize.

        ix = 1
        norm = tol + 1
        avgnorm = norm

        #Main loop 
        while ((norm > tol) && (avgnorm > 0.0001 * tol))

            newvalue, newx, isentry = contract(oldvalue, oldx, profit,
                wmax, two_n, kmax, nfirms, mask,
                x_entryl, x_entryh, phi, entry_k, beta, delta, a,
                state_to_index, index_to_state)

            norm = max(max(abs(oldvalue - newvalue)))
            avgnorm = mean(mean(abs(oldvalue - newvalue)))

            ix += 1
            oldx = newx
            oldvalue = newvalue

        end

        # Warning check
        flag = any(
            index_to_state[s][1] == kmax && newx[s, 1] > 0
            for s in 1:wmax
        )

        if flag
            println("Warning: Positive investment recorded at highest efficiency level.")
            println("Please consider increasing the maximum efficiency level (kmax).")
        end

        # Prising
        prising = a .* newx ./ (1 .+ a .* newx)

        # Save
        using JLD2
        @save "a.$(c.prefix)_markov$(nfirms).jld2" newvalue newx prising isentry

        c.eql_done = true
    end

end
