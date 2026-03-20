####################

#Run main pieces
include("calcval_ps2.jl")
include("optimize_ps2.jl")
include("contract_ps2.jl")
include("initialize_ps2.jl")



#Run function 
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

            if nfirms == 1 && ix == 1
                println("\nTask 2 verification — value function after first iteration N=1:")
                println("Value at w=10: ", newvalue_n[10, 1])
                println("Value at w=19: ", newvalue_n[19, 1])

                # export directly here
                open("task2_values.tex", "w") do f
                    write(f, "After the first contraction iteration with \$N=1\$ and \$\\bar{w}=19\$:\n")
                    write(f, "\\begin{itemize}\n")
                    write(f, "\\item Value function at \$w=10\$: $(round(newvalue_n[10,1], digits=4))\n")
                    write(f, "\\item Value function at \$w=19\$: $(round(newvalue_n[19,1], digits=4))\n")
                    write(f, "\\end{itemize}\n")
                end
                println("Task 2 exported to task2_values.tex")
            end

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


###############################EOF