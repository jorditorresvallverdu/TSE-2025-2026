
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

###############################EOF