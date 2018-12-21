using Distributions
using Plots
using ProgressMeter
using DifferentialEquations
using Statistics
using DelimitedFiles

function ODEmodel(du,u,params,t)

    (xC,xD) = u
    (alphaC,alphaD,sigma,beta,kappa,n,dC,dD,eta,Nt) = params

    # Compute expected public good benefit
    y = xC/(xC+xD)

    pC = beta/n*(1.0 + (n-1)*y)
    pD = beta/n*(n-1)*y

    (rC,rD) = nonlinearPublicGoodForm(pC,pD,params)

    du[1] = rC*xC*(1.0 - (xC+xD)/eta) - dC*xC
    du[2] = rD*xD*(1.0 - (xC+xD)/eta) - dD*xD


end

function cumsum(x)
    y = Array{Float64}(undef,length(x))
    for n = 1 : length(x)
        y[n] = sum(x[1:n])
    end

    y = y./y[end]

    return y
end

function computePublicGoodBenefit(c,d,n,beta)

    if c-1 >= 0 && d >= 0 && n > 1
        ############ A cooperator sees n-1 other cells: ############
        dist = Hypergeometric(c-1,d,n-1)
        cdfdist = cdf.(Ref(dist),support(dist))

        # Find the number of cooperators in the group
        roll = rand()
        cInGroup = findall(x->(x>roll),cdfdist)[1]-1

        if d < n-1         # Less in population than group size
            cInGroup += n-d-1
        end
    
        # Public good benefit
        pC = beta/n*(1.0 + cInGroup)

        # println(cInGroup)

    else
        pC = beta/n
    end
    if c >=0 && d-1 >= 0 && n > 1
        ############ A defector sees n-1 other cells: ############
        # Probability of picking individuals without replacement
        dist = Hypergeometric(c,d-1,n-1)
        cdfdist = cdf.(Ref(dist),support(dist))

        # Find the number of cooperators in the group
        roll = rand()
        cInGroup = findall(x->(x>roll),cdfdist)[1]-1

        if d < n-1         # Less in population than group size
            cInGroup += n-d-1
        end

        # Public good benefit
        pD = beta/n*cInGroup
    else
        pD = 0.0
    end

    # println(cInGroup)

    return (pC,pD)

end

function nonlinearPublicGoodForm(pC,pD,params)
    
    (alphaC,alphaD,sigma,beta,kappa,n,dC,dD,eta) = params

    # Compute rC and rD
    rC = alphaC*(1+exp(sigma))/(1+exp(sigma-pC)) - kappa
    rD = alphaD*(1+exp(sigma))/(1+exp(sigma-pD))

    return (rC,rD)

end


function transitionRates(x,params)

    # Grab populations and parameters
    (c, d) = x
    (alphaC,alphaD,sigma,beta,kappa,n,dC,dD,eta,Nt) = params

    count = 0
    # pCvec = Vector{Float64}
    # pDvec = Vector{Float64}
    # t_10  = Vector{Float64}


    # Initialize the rates array
    rates = Array{Float64}(undef,4)

    cVec = Array{Float64}(undef,Nt); cVec[1] = c
    dVec = Array{Float64}(undef,Nt); dVec[1] = d
    pCvec = Array{Float64}(undef,Nt-1);
    pDvec = Array{Float64}(undef,Nt-1);
    t = Array{Float64}(undef,Nt); t[1] = 0.0
    
    for i = 2 : Nt

        # Compute public good benefit
        (pC,pD) = computePublicGoodBenefit(c,d,n,beta)
        # Compute rC and rD
        (rC,rD) = nonlinearPublicGoodForm(pC,pD,params)

        # print(rC," ",rD,"\n")

        # Rate function (includes all the possible transitions)
        rates[1] = max(0,rC*c*(1.0 - (c + d)/eta))      # Birth of a cooperator
        rates[2] = max(0,rD*d*(1.0 - (c + d)/eta))      # Birth of a defector
        rates[3] = dC*c                                 # Death of a cooperator
        rates[4] = dD*d                                 # Death of a defector

        # Calculate the ratetotal this will be used to compute time till next event
        # and which event occurs
        ratetotal = sum(rates)

        # println(rates)

        # Calculate time till next event
        roll = rand()
        t[i] = t[i-1] - log(roll)/ratetotal

        # Cumulative sum of the rates
        cumrates = cumsum(rates)

        
        # Now we roll to see which event will occur
        roll = rand()
        # print(roll,"\n")
        index = findall(x->(x>roll),cumrates)[1]

        decreaseC = [3;5;6]     # Events which lead to decrease in cooperators
        decreaseD = [4;7;8]     # Events which lead to decrease in defectors

        # Change the population size
        if any(decreaseC .== index)        # C -> C - 1
            c -= 1
            # println("c dies\n")
        elseif any(decreaseD .== index)    # D -> D - 1
            d -= 1
            # println("d dies\n")
        elseif index == 1                           # C -> C + 1
            c += 1
            # println("c born\n")
        else
            d += 1
            # println("d born\n")
        end

        cVec[i] = c
        dVec[i] = d

        pCvec[i-1] = pC
        pDvec[i-1] = pD

        if c*d == 0
            t    = t[1:i]
            cVec = cVec[1:i]
            dVec = dVec[1:i]
            pCvec = pCvec[1:i]
            pDvec = pDvec[1:i]
            break
        end

        if i % 10000 == 0
            # println(pC," ", pD,"\n")
        end

    end

return (t,cVec,dVec,pCvec,pDvec)

end

# (C0,D0)
C0 = 50; D0 = 50;
x = (C0,D0)

betaval = 5.0; n = 15; sigma = 2.0;
Nt = 10^5
plotsteps = Int(max(Nt/1000,100))

# # (alphaC,alphaD,sigma,beta,kappa,n,dC,dD,eta,Nt)
# params = (1.0,1.0,2.0,betaval,0.5,n,0.1,0.1,1000.0,Nt)

Nruns = 200

savematrix = Array{Any}(undef,Nruns)

params = (1.0,1.0,sigma,betaval,0.5,n,0.1,0.1,1000.0,Nt)

@showprogress 1 "Computing..." for i = 1 : Nruns
    
    (t,cooperators,defectors,pCvec,pDvec) = transitionRates(x,params)

    t = t[1:plotsteps:end];
    cooperators = cooperators[1:plotsteps:end]
    defectors = defectors[1:plotsteps:end]
    savematrix[i] = [t cooperators defectors]

end

outfile = string("/Users/gregorykimmel/Dropbox/02_nonlinearPGG/manuscript/",
"Figures/SIF1/code/manytrajectories_C0_",C0,"_D0_",D0,".txt")

writedlm(outfile,savematrix)

u0 = [float(C0),float(D0)]
tspan = (0.0,1000.0)

prob = ODEProblem(ODEmodel,u0,tspan,params)
sol = solve(prob,Tsit5())
plotting = false

Nt_ODE = length(sol.t)
xC      = Array{Float64}(undef,Nt_ODE)
xD      = Array{Float64}(undef,Nt_ODE)
y       = Array{Float64}(undef,Nt_ODE)
pC_ODE  = Array{Float64}(undef,Nt_ODE)
pD_ODE  = Array{Float64}(undef,Nt_ODE)
for i = 1 : Nt_ODE
    xC[i] = sol.u[i][1]
    xD[i] = sol.u[i][2]
    y[i]  = xC[i]/(xC[i]+xD[i])
    pC_ODE[i] = betaval/n*(1.0+y[i]*(n-1))
    pD_ODE[i] = betaval/n*y[i]*(n-1)
end

outfile = string("/Users/gregorykimmel/Dropbox/02_nonlinearPGG/manuscript/",
"Figures/SIF1/code/ODEtrajectory_C0_",C0,"_D0_",D0,".txt")

writedlm(outfile,[sol.t xC])

# neighborhoodSizeMax = 30

# cooperatorDeath = zeros(neighborhoodSizeMax)
# defectorDeath = zeros(neighborhoodSizeMax)
# coexistence = zeros(neighborhoodSizeMax)
# meancooperators = zeros(neighborhoodSizeMax)
# meandefectors = zeros(neighborhoodSizeMax)

# @showprogress 1 "Computing..."  for i = 1 : Nruns
#     for j = 1 : neighborhoodSizeMax
#         params = (1.0,1.0,sigma,betaval,0.5,j,0.1,0.1,1000.0,Nt)
#         (t,cooperators,defectors,pCvec,pDvec) = transitionRates(x,params)

#         # If cooperators have died
#         if cooperators[end] == 0
#             cooperatorDeath[j] += 1
#         elseif defectors[end] == 0
#             defectorDeath[j] += 1
#         else
#             coexistence[j] += 1
#             startingIndex = Int(0.9*length(t))
#             meancooperators[j] = mean(cooperators[startingIndex:end])
#             meandefectors[j] = mean(defectors[startingIndex:end])
#         end
#     end
# end

# for i = 1 : Nruns
#     for j = 1 : neighborhoodSizeMax
#         params = (1.0,1.0,sigma,betaval,0.5,j,0.1,0.1,1000.0,Nt)
#         (t,cooperators,defectors,pCvec,pDvec) = transitionRates(x,params)

#         # If cooperators have died
#         if cooperators[end] == 0
#             cooperatorDeath[j] += 1
#         elseif defectors[end] == 0
#             defectorDeath[j] += 1
#         else
#             coexistence[j] += 1
#             startingIndex = Int(0.9*length(t))
#         end
#         meancooperators[j] = mean(cooperators[startingIndex:end])
#         meandefectors[j] = mean(defectors[startingIndex:end])
#     end
# end