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



#Export, exercise 1. # Task 1 verification — decode table for N=3, wbar=3
small_states = states(3)
small_i2s = Dict(v => k for (k, v) in small_states)
n_states = length(small_i2s)
table = hcat([collect(small_i2s[i]) for i in 1:n_states]...)
println("\nTask 1 — Decode table for N=3, w̄=3:")
println(table)

# Export Task 1 table to latex
open("task1_table.tex", "w") do f
    write(f, "\\begin{table}[h]\n")
    write(f, "\\centering\n")
    write(f, "\\resizebox{\\textwidth}{!}{%\n")  # resize since table is wide
    write(f, "\\begin{tabular}{l$(repeat("c", n_states))}\n")
    write(f, "\\hline\n")
    write(f, "State code & " * join(1:n_states, " & ") * " \\\\\n")
    write(f, "\\hline\n")
    row_labels = ["\$w_1\$", "\$w_2\$", "\$w_3\$"]
    for row in 1:3
        write(f, row_labels[row] * " & " * join(table[row, :], " & ") * " \\\\\n")
    end
    write(f, "\\hline\n")
    write(f, "\\end{tabular}}\n")
    write(f, "\\caption{Decode table for \$N=3\$, \$\\bar{w}=3\$. Column \$i\$ represents the weakly descending 3-tuple encoded as state code \$i\$.}\n")
    write(f, "\\label{tab:decode}\n")
    write(f, "\\end{table}\n")
end
println("Task 1 table exported to task1_table.tex")


#########EOF
