
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



############################EOF