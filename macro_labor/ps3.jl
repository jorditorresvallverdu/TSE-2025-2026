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

S = p_matrix .- b    # reshape b for correct broadcasting

tol = 1e-6
crit = 1.0

while crit > tol
    S_new = p_matrix .- b .+ (1 - δ) * β .* max.(S, 0)
    crit = maximum(abs.(S_new .- S))    # S not S_not
    S = S_new
end


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

# PATCH: added searcher_x argument (replaces fx as weight)
function vacancy_dist(q, searcher_x)

    v = zeros(ny)

    for j in 1:ny
        benefit = 0.0
        for i in 1:nx
            # PATCH: fx[i] -> searcher_x[i]
            benefit += (S[i, j] > 0) * (1 - β * (1 - δ)) * S[i, j] * searcher_x[i]
        end
        v[j] = benefit > 0 ? ((q / c0) * benefit)^(1 / c1) : 0.0
    end

    V = sum(v)

    return v, V
end

function distributions(v, V, λ, s, S, δ)
    g = zeros(nx, ny)
    u = zeros(nx)
    e = zeros(nx, ny)

    for i in 1:nx
        # PATCH: sort by surplus so inflows from worse matches are ready
        sorted_j = sortperm(S[i, :])

        for j in sorted_j
            if S[i, j] <= 0
                g[i, j] = 0.0
                continue
            end

            inflow = λ * (v[j] / V) * (S[i, j] > 0)
            # PATCH: add inflow from worse matches (on-the-job search)
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

        # u and e computed AFTER g[i,:] is fully filled
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


function fixed_point(q0, λ0)

    ϵ = 1e-6
    crit = 1.0
    q = q0
    λ = λ0

    # declare outside loop so they are accessible at return
    v = zeros(ny)
    V = 0.0
    u = zeros(nx)
    e = zeros(nx, ny)

    # PATCH: initialize searcher weights
    searcher_x = copy(fx)

    damp = 0.3  # PATCH: damping for stability
    iter = 0

    while crit > ϵ
        iter += 1
        if iter > 2000
            println("Warning: did not converge after 2000 iterations, crit=$crit")
            break
        end

        # PATCH: pass searcher_x
        v, V = vacancy_dist(q, searcher_x)

        u, e = distributions(v, V, λ, s, S, δ)

        λ_new = update_lambda(u, e, α, ω, V)

        # PATCH: correct formula q = m(L,V)/V = α*(L/V)^ω
        L = sum(u) + s * sum(e)
        q_new = α * (L / V)^ω

        # PATCH: update searcher weights
        for i in 1:nx
            searcher_x[i] = (u[i] + s * sum(e[i, :])) / L
        end

        crit = max(abs(q_new - q), abs(λ_new - λ))

        # PATCH: damped update
        q = damp * q_new + (1 - damp) * q
        λ = damp * λ_new + (1 - damp) * λ

    end

    println("Finished in $iter iterations, crit=$crit")
    # Add this right after fixed_point converges, to diagnose:
    for j in 1:ny
        ben = 0.0
        for i in 1:nx
            ben += (S[i, j] > 0) * (1 - β * (1 - δ)) * S[i, j] * searcher_x[i]
        end
        println("y=$(round(y_grid[j],digits=2))  benefit=$ben  v=$(v[j])")
    end
    return q, λ, v, V, u, e


end


q0 = 0.001
λ0 = 0.001
damp = 0.1

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



[q, λ]

