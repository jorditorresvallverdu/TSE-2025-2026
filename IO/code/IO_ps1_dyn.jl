# Author: Jordi Torres Vallverdú
# Start Date: 12/03/2025

############ Main Function ############

mutable struct ModelParams
    max_firms::Int
    kmax::Int
    start_firms::Int

    entry_low::Float64
    entry_high::Float64
    scrap_val::Float64
    entry_at::Int
    beta::Float64
    delta::Float64

    inv_mult::Float64

    intercept::Float64
    fixed_cost::Float64
    gamma::Float64

    tol::Float64

    profit_done::Bool
    eql_done::Bool
    prefix::String

    ds_wstart::Vector{Float64}
    ds_nsimx::Int
end

    #### Model Parameters ####

    max_firms = 3
    kmax = 19
    start_firms = 1

    entry_low = 0.15
    entry_high = 0.25
    scrap_val = 0.1
    entry_at = 4
    beta = 0.925
    delta = 0.7

    inv_mult = 3

    intercept = 3
    fixed_cost = 0.2
    gamma = 1

    tol = 0.1

    profit_done = false
    eql_done = false
    prefix = "cc"

    ds_wstart = vcat(entry_at + 2, zeros(max_firms - 1))
    ds_nsimx = 10000


    c = ModelParams(
        max_firms,
        kmax,
        start_firms,
        entry_low,
        entry_high,
        scrap_val,
        entry_at,
        beta,
        delta,
        inv_mult,
        intercept,
        fixed_cost,
        gamma,
        tol,
        profit_done,
        eql_done,
        prefix,
        ds_wstart,
        ds_nsimx
    )


    ######## Baseline ########

    #static_profit(c)
    #eql_ma(c)
    #ds_ma(c, "baseline.tex")


    ######## Low Entry Cost ########

    c.entry_low = 0.01
    c.entry_high = 0.11

    #static_profit(c)
    #eql_ma(c)
    #ds_ma(c, "low_entry_cost.tex")

####FUNCTION 1: STATIC PROFITS


function static_profit(c)

    nfmax= c.max_firms
    kkmax= c.kmax

    #Set up binomial coefficients



end



#This generates the states and the index. Probably is not so efficient as the matlab computation, but it is much more simple to compute
function states(kmax)
    states= Dict()  
    id=1
    for w1 in 0:kmax
        for w2 in 0:w1
            for w3 in 0:w2
                states[(w1,w2,w3)]= id 
                id +=1    
            end 
        end 
    end 
    return states
end 
data_states= states(kmax)
println(data_states)
#Now we define the connection between the function of the states and the numeric index. It should be easy to check as it should be indexed easily. 
println(length(data_states))

#I guess this is a version of the decode function, if needed afterwards.
index_to_state = Dict(v => k for (k,v) in data_states)
state_to_index = Dict(k => v for (k,v) in data_states)

println(index_to_state)

function theta(gamma, w)
    cost= @. gamma * exp(-(w-4.0)) #this should be in vectorized form for slight efficiency in terms of coding 
    return cost
end 



function static_profit_computation(data_states, c)

#Declare parameter values. 

num_data_states= length(data_states)

    for states in 1:num_data_states

        #Computes marginal cost : c1, c2, c3 
        w1, w2, w3= index_to_state[states]
        c1= theta(c.gamma, w1)
        if w2>0
        c2= theta(c.gamma, w2)
        end 
        if w3>0
        c3= theta(c.gamma, w3)
        end 

        theta= vcat(c1, c2, c3)
        #Determine active firms
        n_firms= count(c1,c2,c3) #different from 0

        p= (D + () )


    end


    

end 

    D = intercept;  


        #Computes marginal cost : c1, c2, c3 
        w1, w2, w3 = index_to_state[10] #this is just an example

        w= vcat(w1, w2, w3) # just to provide it in vectorized form
        n_firms= count(w.>0) #count the elements of the vector that are larger than 0

        cost= theta(1, w)
        #Determine active firms

        #here we need to use the fact w1>w2>w3 -> theta1< theta2< theta3
        p= (D + sum(cost) )/(n_firms+1)

        #here we need to generate top profits. 
        function find_eq()
            while (p - cost[1] >= 0)                    

            #finish here the function...

