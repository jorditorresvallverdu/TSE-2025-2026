##Goal: This is to answer the third problem set of Macro. 

# ============================================================
# Parameters
# ============================================================

using Plots, Statistics
# Discount and separation
const β = 0.99
const δ = 0.04

# Matching function
const α = 0.5
const ω = 0.5
const s = 0.1

# Vacancy cost
const c0 = 0.01
const c1 = 0.11

# Production function coefficients
const p1 = 0.2
const p2 = 1.0
const p3 = -4.0
const p4 = 8.0
const p5 = -2.2
const p6 = 7.0

# Home production
const b0 = 0.7

# Grids
const nx = 20
const ny = 20


#########################

#Define functions 

function p(x, y)
    p = p1 + p2 * x + p3 * y + p4 * (x^2) + p5 * (y^2) + p6 * x * y
    return p
end


function c(v) #doubt is why do we need vacacny post in this model??
    num = v^(1 + c1)
    den = 1 + c1
    tot = c0 * (num / den)
    return tot
end


#distribution / grids 
y_grid = range(0, 1, length=ny) |> collect #collect is just to call it a vector.
fy = fill(1 / ny, ny)

x_grid = range(0, 1, length=nx) |> collect
fx = [(nx - i)^3 for i in 0:(nx-1)]
fx = fx ./ sum(fx)


#define U, this should be fixed.

b = [b0 * maximum(p.(x_grid[i], y_grid)) for i in 1:nx]

p_matrix = [p(x_grid[i], y_grid[j]) for i in 1:nx, j in 1:ny]


function fixed_point_s(p_matrix, b)
    S = p_matrix .- b    # reshape b for correct broadcasting

    tol = 1e-6
    crit = 1.0

    while crit > tol
        S_new = p_matrix .- b .+ (1 - δ) * β .* max.(S, 0)
        crit = maximum(abs.(S_new .- S))    # S not S_not
        S = S_new
    end

    return S

end

# FIX 1: call fixed_point_s to get S globally (was never called before)
S = fixed_point_s(p_matrix, b)

# Plot 1: matching set — area where surplus is positive
match_set = S .> 0    # nx × ny boolean matrix

plot1 = heatmap(y_grid, x_grid, match_set,
    xlabel="Firm productivity (y)",
    ylabel="Worker ability (x)",
    title="Matching Set: S(x,y) > 0",
    color=:RdYlGn)

# Plot 2a: surplus as contour
plot2 = contourf(y_grid, x_grid, S,
    xlabel="Firm productivity (y)",
    ylabel="Worker ability (x)",
    title="Surplus S(x,y) — Contour",
    color=:viridis)

# Plot 2b: surplus as 3D surface
plot3 = surface(y_grid, x_grid, S,
    xlabel="y",
    ylabel="x",
    zlabel="S(x,y)",
    title="Surplus S(x,y) — 3D Surface",
    color=:viridis)


###Exercise 2. 


#Coding prior functions

# FIX 2: added S as argument so Exercise 4 can pass its own S
function vacancy_dist(q, searcher_x, S, c0_val=c0)

    v = zeros(ny)

    for j in 1:ny
        benefit = 0.0
        for i in 1:nx
            benefit += (S[i, j] > 0) * (1 - β * (1 - δ)) * S[i, j] * searcher_x[i]
        end
        v[j] = benefit > 0 ? ((q / c0_val) * benefit)^(1 / c1) : 0.0
    end

    V = sum(v)

    return v, V
end

function distributions(v, V, λ, s, S, δ)
    g = zeros(nx, ny)
    u = zeros(nx)
    e = zeros(nx, ny)

    for i in 1:nx
        sorted_j = sortperm(S[i, :])

        for j in sorted_j
            if S[i, j] <= 0
                g[i, j] = 0.0
                continue
            end

            inflow = λ * (v[j] / V) * (S[i, j] > 0)
            for k in 1:ny
                if S[i, k] > 0 && S[i, k] < S[i, j]
                    inflow += s * λ * (v[j] / V) * g[i, k]
                end
            end

            outflow = δ
            for k in 1:ny
                outflow += s * λ * v[k] / V * (S[i, k] > S[i, j])
            end
            g[i, j] = inflow / outflow
        end

        u[i] = fx[i] / (1 + sum(g[i, :]))
        e[i, :] = u[i] .* g[i, :]
    end

    return u, e
end

function update_lambda(u, e, α, ω, V)

    L = sum(u) + s * sum(e)

    λ_new = α * L^ω * V^(1 - ω) / L

    return λ_new
end


# FIX 3: added S and c0_val as arguments with defaults so existing call still works
function fixed_point(q0, λ0, S=S, c0_val=c0)

    ϵ = 1e-6
    crit = 1.0
    q = q0
    λ = λ0

    v = zeros(ny)
    V = 0.0
    u = zeros(nx)
    e = zeros(nx, ny)

    searcher_x = copy(fx)

    damp = 0.3
    iter = 0

    while crit > ϵ
        iter += 1
        if iter > 2000
            println("Warning: did not converge after 2000 iterations, crit=$crit")
            break
        end

        v, V = vacancy_dist(q, searcher_x, S, c0_val)

        u, e = distributions(v, V, λ, s, S, δ)

        λ_new = update_lambda(u, e, α, ω, V)

        L = sum(u) + s * sum(e)
        q_new = α * (L / V)^ω

        for i in 1:nx
            searcher_x[i] = (u[i] + s * sum(e[i, :])) / L
        end

        crit = max(abs(q_new - q), abs(λ_new - λ))

        q = damp * q_new + (1 - damp) * q
        λ = damp * λ_new + (1 - damp) * λ

    end

    println("Finished in $iter iterations, crit=$crit")
    return q, λ, v, V, u, e

end


q0 = 0.001
λ0 = 0.001

q, λ, v, V, u, e = fixed_point(q0, λ0)

println("Unemployment rate = $(sum(u))")
println("Total employment = $(sum(e))")

# Exercise 2: employment distribution
plot4 = heatmap(y_grid, x_grid, e,
    xlabel="Firm productivity (y)",
    ylabel="Worker ability (x)",
    title="Employment distribution e(x,y)",
    color=:viridis)

# Exercise 3
plot5 = bar(y_grid, v, xlabel="y", ylabel="v(y)", title="Vacancy distribution")
plot6 = bar(x_grid, u, xlabel="x", ylabel="u(x)", title="Unemployment distribution")


# Exercise 4

function solve_model(p1_new, c0_new)

    # FIX 5: was p_new, should be p1_new
    p_matrix_new = [p1_new + x_grid[i] + y_grid[j] for i in 1:nx, j in 1:ny]

    b_new = [b0 * maximum(p_matrix_new[i, :]) for i in 1:nx]

    S_new = fixed_point_s(p_matrix_new, b_new)

    # FIX 6: pass S_new and c0_new into fixed_point (now accepted as arguments)
    q, λ, v, V, u, e = fixed_point(q0, λ0, S_new, c0_new)

    Y = sum(e .* p_matrix_new)
    U = sum(u)

    return Y, U

end

#define target of calibration
Y_target = sum(e .* p_matrix)
U_target = sum(u)

#define grid 
grid_p1new = range(-2, 2, length=100)
grid_c0new = range(0.001, 0.1, length=100)

function calibration(Y_target, U_target)

    best_loss = Inf
    best_p1 = 0.0
    best_c0 = 0.0

    for p1_new in grid_p1new
        for c0_new in grid_c0new
            Y, U = solve_model(p1_new, c0_new)
            loss = (Y - Y_target)^2 / Y_target^2 + (U - U_target)^2 / U_target^2
            if loss < best_loss
                best_loss = loss
                best_p1 = p1_new
                best_c0 = c0_new
            end
        end
    end

    return best_p1, best_c0
end

best_p1, best_c0 = calibration(Y_target, U_target)



####Print tables for 2 and 4

# ============================================================
# Generate LaTeX tables
# ============================================================

# Run the linear model with best calibration
p_matrix_lin = [best_p1 + x_grid[i] + y_grid[j] for i in 1:nx, j in 1:ny]
b_lin = [b0 * maximum(p_matrix_lin[i, :]) for i in 1:nx]
S_lin = fixed_point_s(p_matrix_lin, b_lin)
q_lin, λ_lin, v_lin, V_lin, u_lin, e_lin = fixed_point(q0, λ0, S_lin, best_c0)

Y_base = sum(e .* p_matrix)
U_base = sum(u)
Y_lin = sum(e_lin .* p_matrix_lin)
U_lin = sum(u_lin)

# Table 1: Exercise 2 — Baseline equilibrium
tab1 = """
\\begin{table}[htbp]
\\centering
\\caption{Baseline Equilibrium (Complementarities)}
\\begin{tabular}{lc}
\\hline\\hline
Variable & Value \\\\
\\hline
Unemployment rate & $(round(U_base, digits=6)) \\\\
Total employment & $(round(1 - U_base, digits=6)) \\\\
Total output & $(round(Y_base, digits=4)) \\\\
Contact rate (\\lambda) & $(round(λ, digits=4)) \\\\
Vacancy filling rate (q) & $(round(q, digits=6)) \\\\
Total vacancies (V) & $(round(V, digits=4)) \\\\
Fraction S>0 & $(round(sum(S .> 0)/length(S), digits=4)) \\\\
\\hline\\hline
\\end{tabular}
\\end{table}
"""

# Table 2: Exercise 4 — Comparison
tab2 = """
\\begin{table}[htbp]
\\centering
\\caption{Comparison: Complementarities vs.\\ Linear Production}
\\begin{tabular}{lcc}
\\hline\\hline
 & Complementarities & Linear \\\\
\\hline
\$p_1\$ & $(p1) & $(round(best_p1, digits=4)) \\\\
\$c_0\$ & $(c0) & $(round(best_c0, digits=6)) \\\\
Unemployment rate & $(round(U_base, digits=6)) & $(round(U_lin, digits=6)) \\\\
Total output & $(round(Y_base, digits=4)) & $(round(Y_lin, digits=4)) \\\\
Contact rate (\$\\lambda\$) & $(round(λ, digits=4)) & $(round(λ_lin, digits=4)) \\\\
Vacancy filling rate (\$q\$) & $(round(q, digits=6)) & $(round(q_lin, digits=6)) \\\\
Total vacancies (\$V\$) & $(round(V, digits=4)) & $(round(V_lin, digits=4)) \\\\
Fraction \$S>0\$ & $(round(sum(S .> 0)/length(S), digits=4)) & $(round(sum(S_lin .> 0)/length(S_lin), digits=4)) \\\\
\\hline\\hline
\\end{tabular}
\\end{table}
"""

# Write to files
open("table_baseline.tex", "w") do f
    write(f, tab1)
end
open("table_comparison.tex", "w") do f
    write(f, tab2)
end

println("Tables saved to table_baseline.tex and table_comparison.tex")

####eof