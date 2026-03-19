# Author: Jordi Torres Vallverdú
# Start Date: 12/03/2025

using JLD2
using Statistics

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
println(length(data_states))

#I guess this is a version of the decode/encode function, if needed afterwards.
index_to_state = Dict(v => k for (k, v) in data_states)
state_to_index = Dict(k => v for (k, v) in data_states)

# Now we define cost function
function theta(gamma, w)
    cost = @. gamma * exp(-(w - 4.0)) #this should be in vectorized form for slight efficiency in terms of coding
    return cost
end

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
function calcval(place, w, x, k,
    oldvalue,
    two_n, kmax, nfirms, mask,
    delta, a,
    state_to_index)

    # w and x are always length 3 (from the 3-tuple dictionary)
    # z1, z2 bounds must also be length 3
    z1 = zeros(Int, 3)
    z2 = fill(kmax, 3)

    # Build locmask: must have 3 rows to match w which is always length 3.
    # The rivals are the other nfirms-1 active firms, but inactive slots stay 0.
    # mask has (nfirms-1) rows; we embed it into a (3, two_n) matrix.
    locmask = zeros(Int, 3, two_n)
    if nfirms > 1
        # fill rows for rivals only (exclude the focal firm's row)
        rival_rows = vcat(1:place-1, place+1:nfirms)  # rows of the rivals in the nfirms world
        for (mi, ri) in enumerate(rival_rows)
            locmask[ri, :] = mask[mi, :]
        end
    end

    x = copy(x)
    w = copy(w)

    x[place] = 0
    w[place] = k

    justone = zeros(Int, 3)
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
            (1 - delta) * oldvalue[idx_d, min(pl1, nfirms)] +
            delta * oldvalue[idx_e, min(pl1, nfirms)]
        ) * probmask

        # -------- MOVE UP --------
        d_up = w + locmask[:, i]
        d_up[place] += 1

        e_up, d_up_sorted, pl1_up = state_helper(d_up, justone, z1, z2)

        idx_d_up = state_to_index[Tuple(sort(d_up_sorted, rev=true))]
        idx_e_up = state_to_index[Tuple(sort(e_up, rev=true))]

        val_up += (
            (1 - delta) * oldvalue[idx_d_up, min(pl1_up, nfirms)] +
            delta * oldvalue[idx_e_up, min(pl1_up, nfirms)]
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

    # Decode state — always a 3-tuple, always length 3
    locw = collect(index_to_state[w])
    locwx = copy(locw)

    oval = copy(oldvalue[w, 1:nfirms])
    ox = vcat(oldx[w, 1:nfirms], zeros(3 - nfirms))

    nval = zeros(nfirms)
    nx = zeros(nfirms)

    #################### 1. Exit condition ####################
    profit_w = profit[w, :]
    CV = profit_w .+ beta .* oval

    for s in 1:nfirms
        if CV[s] >= phi
            locwx[s] = locw[s]
        else
            locwx[s:end] .= 0  # equilibrium selection: if s exits, all weaker firms exit too
            break
        end
    end

    locwe = copy(locwx)

    #################### 2. Entry condition ####################
    p_entry = 0.0
    n_active = sum(locwe .> 0)

    if n_active >= nfirms  # market full, no entry possible
        p_entry = 0.0
    elseif n_active + 1 > 3  # safety: can't exceed dictionary tuple size
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
            oval[j] = phi
            ox[j] = 0.0
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
        locw = collect(index_to_state[w])  # always length 3

        if locw[nfirms] == 0

            # figure out where the potential entrant would sit
            n_active = sum(locw .> 0)
            place = n_active + 1

            # build the state with the entrant added at entry_k, then sort
            locw_entry = copy(locw)
            locw_entry[place] = entry_k
            locw_entry_sorted = sort(locw_entry, rev=true)

            # find where the entrant ended up after sorting (for calcval)
            pl_entry = findfirst(x -> x == entry_k, locw_entry_sorted)
            if pl_entry === nothing
                pl_entry = place
            end

            # entrant has no prior investment — length 3 to match w
            ox_entry = zeros(3)

            _, val_stay = calcval(pl_entry, locw_entry_sorted, ox_entry, entry_k,
                oldvalue, two_n, kmax, nfirms, mask, delta, a, state_to_index)

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
        newvalue[w, :], newx[w, :] = optimize(w, oldvalue, oldx, isentry, profit,
            wmax, two_n, kmax, nfirms, mask,
            phi, entry_k, beta, delta, a,
            state_to_index, index_to_state)
    end

    return newvalue, newx, isentry
end

#######################################################################################
############################# PART 5. #################################################
function make_mask(nfirms)
    n = nfirms - 1
    two_n = 2^n

    if n == 0
        return zeros(Int, n, two_n), two_n  # (0 x 1) — no rivals, no rows needed
    end

    return hcat([reverse(digits(i, base=2, pad=n)) for i in 0:two_n-1]...), two_n
end

# Translated from Matlab initialize.m --> THIS WAS DONE directly with Claude, due to lack of time. 
function initialize(nfirms, wmax, index_to_state, state_to_index, newvalue, newx)
    oldx = zeros(wmax, nfirms)
    oldvalue = zeros(wmax, nfirms)

    if nfirms == 1
        oldvalue[:, 1] = 1 .+ 0.1 .* (1:wmax)  # newvalue and newx don't exist for nfirms-1=0
    else
        for w in 1:wmax
            tuple_w = collect(index_to_state[w])

            # Initialize by mapping the current w to the corresponding nfirms-1 equilibrium
            # (ignore the last firm) — pad to 3-tuple for dictionary lookup
            sub_tuple = Tuple(vcat(tuple_w[1:nfirms-1], zeros(Int, 3 - (nfirms - 1))))
            n_sub = state_to_index[sub_tuple]
            oldvalue[w, 1:nfirms-1] = newvalue[n_sub, 1:nfirms-1]  # FIX: slice rhs too
            oldx[w, 1:nfirms-1] = newx[n_sub, 1:nfirms-1]  # FIX: slice rhs too

            # Initialize the last firm by ignoring the second last firm,
            # i.e., swap the last two firms and set the last firm to 0 state
            tuple_mod = copy(tuple_w)
            tuple_mod[nfirms-1] = tuple_w[nfirms]
            tuple_mod[nfirms] = 0
            n_mod = state_to_index[Tuple(sort(tuple_mod, rev=true))]
            oldvalue[w, nfirms] = oldvalue[n_mod, nfirms-1]
            oldx[w, nfirms] = oldx[n_mod, nfirms-1]
        end
    end

    return oldvalue, oldx
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

    D = c.intercept
    f = c.fixed_cost
    γ = c.gamma

    # Build state dictionaries
    local_states = states(kmax)
    local_i2s = Dict(v => k for (k, v) in local_states)
    local_s2i = Dict(k => v for (k, v) in local_states)
    wmax = length(local_i2s)

    # Compute static profits
    profit_full = ccprofit_all(D, f, local_i2s, γ)

    # Full containers — always (wmax, rlnfirms) so initialize can read correct slices
    newvalue = zeros(wmax, rlnfirms)
    newx = zeros(wmax, rlnfirms)
    isentry = zeros(wmax)

    for nfirms in stfirm:rlnfirms

        profit = profit_full[:, 1:nfirms]

        #Number of rival actions + mask
        mask, two_n = make_mask(nfirms)

        oldvalue, oldx = initialize(nfirms, wmax, local_i2s, local_s2i, newvalue, newx)

        ix = 1
        norm = tol + 1
        avgnorm = norm

        #Main loop — use temporaries so newvalue/newx stay (wmax, rlnfirms) between nfirms iterations
        while ((norm > tol) && (avgnorm > 0.0001 * tol))

            newvalue_n, newx_n, isentry = contract(oldvalue, oldx, profit,
                wmax, two_n, kmax, nfirms, mask,
                x_entryl, x_entryh, phi, entry_k, beta, delta, a,
                local_s2i, local_i2s)

            norm = maximum(abs.(oldvalue .- newvalue_n))
            avgnorm = mean(abs.(oldvalue .- newvalue_n))
            println("ix=$ix norm=$norm avgnorm=$avgnorm sample_val=$(newvalue_n[100,1])")

            ix += 1
            oldx = newx_n
            oldvalue = newvalue_n

        end

        # After convergence, store nfirms slice back into full containers for next initialize
        newvalue[:, 1:nfirms] = oldvalue
        newx[:, 1:nfirms] = oldx

        println("nfirms=$nfirms converged after $ix iterations. norm=$norm")

        # Warning check
        flag = any(
            local_i2s[s][1] == kmax && newx[s, 1] > 0
            for s in 1:wmax
        )
        if flag
            println("Warning: Positive investment recorded at highest efficiency level.")
            println("Please consider increasing the maximum efficiency level (kmax).")
        end

        # Prising
        prising = a .* newx ./ (1 .+ a .* newx)

        # Save
        @save "a.$(c.prefix)_markov$(nfirms).jld2" newvalue newx prising isentry

        c.eql_done = true
    end
    return newvalue, newx

end

#######################################################################################
############################# PART 6: DS_MA ############################################

function ds_ma(c)

    # load
    @load "a.$(c.prefix)_markov3.jld2" newvalue newx prising isentry

    # initialize
    w = [c.ds_wstart[1], c.ds_wstart[2], c.ds_wstart[3]]  # (6,0,0)

    # containers
    active_firms = zeros(c.ds_nsimx)
    total_inv = zeros(c.ds_nsimx)
    total_realized_inv = zeros(c.ds_nsimx)
    # main loop
    for t in 1:c.ds_nsimx
        s = state_to_index[Tuple(sort(Int.(w), rev=true))]

        local new_entrant_pos = 0
        # record BEFORE updating
        active_firms[t] = sum(w .> 0)
        total_realized_inv[t] = sum(newx[s, :]) #to distinguish from made

        # 1. EXIT
        w_new = copy(Int.(w))
        for j in 1:3
            if newvalue[s, j] <= c.scrap_val
                w_new[j:end] .= 0
                break
            end
        end

        # 2. ENTRY
        n_active = sum(w_new .> 0)
        if n_active < 3
            s_post_exit = state_to_index[Tuple(sort(w_new, rev=true))]
            if rand() < isentry[s_post_exit]
                w_new[n_active+1] = c.entry_at
                sort!(w_new, rev=true)
                new_entrant_pos = findfirst(x -> x == c.entry_at, w_new)

            end
        end

        # 3. INVESTMENT + 4. AGGREGATE SHOCK
        s_post_entry = state_to_index[Tuple(sort(w_new, rev=true))]

        total_inv[t] = sum(newx[s_post_entry, :]) #to distinguish from realized

        tau = [j == new_entrant_pos ? 0 : Int(rand() < prising[s_post_entry, j]) for j in 1:3]
        nu = rand() < c.delta ? 1 : 0

        # 5. UPDATE
        w_new = w_new .+ tau .- nu
        w_new = max.(w_new, 0)
        w_new = min.(w_new, c.kmax)
        sort!(w_new, rev=true)

        w = w_new
    end


    println("Average active firms: ", mean(active_firms))
    println("Average total investment: ", mean(total_inv))
    println("Realized Investment (-shocks):", mean(total_realized_inv))
end

newvalue, newx = eql_ma(c)


# After eql_ma(c)
# Pick a few interesting states to compare

#Some testing 

# State (6,0,0) - the simulation starting point
s = state_to_index[(6, 0, 0)]
println("State (6,0,0):")
println("  newvalue: ", newvalue[s, :])
println("  newx: ", newx[s, :])

# State (10,5,0) - two active firms
s2 = state_to_index[(10, 5, 0)]
println("State (10,5,0):")
println("  newvalue: ", newvalue[s2, :])
println("  newx: ", newx[s2, :])

# State (10,5,3) - three active firms
s3 = state_to_index[(10, 5, 3)]
println("State (10,5,3):")
println("  newvalue: ", newvalue[s3, :])
println("  newx: ", newx[s3, :])

println("profit (10,5,0): ", mat_results[state_to_index[(10, 5, 0)], :])
println("profit (10,5,3): ", mat_results[state_to_index[(10, 5, 3)], :])


w = [10, 5, 3]
cost = theta(1, w)
println("costs: ", cost)
p = (3 + sum(cost)) / 4  # initial price guess
println("initial price: ", p)
println("p - cost[3]: ", p - cost[3])

#FINAL RESULTS 
ds_ma(c) #baseline 


#CF




###Change the support of c. 