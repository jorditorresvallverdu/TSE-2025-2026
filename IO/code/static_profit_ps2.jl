
####FUNCTION 1: STATIC PROFITS

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



########################################EOF