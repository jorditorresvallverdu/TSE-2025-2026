
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


##################################################EOF