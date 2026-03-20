
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


################################EOF 

