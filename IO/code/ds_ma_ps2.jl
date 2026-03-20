
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

        # record BEFORE updating
        active_firms[t] = sum(w .> 0)
        total_realized_inv[t] = sum(newx[s, :])

        # 1. EXIT
        w_new = copy(Int.(w))
        for j in 1:3
            if newvalue[s, j] <= c.scrap_val
                w_new[j:end] .= 0
                break
            end
        end

        # 2. INVESTMENT — post-exit state, already weights over entry probability
        s_post_exit = state_to_index[Tuple(sort(w_new, rev=true))]
        total_inv[t] = sum(newx[s_post_exit, :])
        tau = [Int(rand() < prising[s_post_exit, j]) for j in 1:3]

        # 3. ENTRY
        n_active = sum(w_new .> 0)
        if n_active < 3
            if rand() < isentry[s_post_exit]
                w_new[n_active+1] = c.entry_at
                sort!(w_new, rev=true)
            end
        end

        # 4. AGGREGATE SHOCK
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
    #println("Realized Investment (-shocks): ", mean(total_realized_inv))
    return mean(active_firms), mean(total_inv)
end

################################# EOF