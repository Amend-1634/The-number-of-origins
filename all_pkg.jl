using Pkg
Pkg.add("Distributions") #all else sampling function
Pkg.add("PoissonRandom") #pois_rand()
Pkg.add("LazyArrays") #ApplyArray(vcat,)
Pkg.add("Multibreak") #@multibreak
Pkg.add("Plots")
Pkg.add("Distributions") #MvNormal
Pkg.add("StatsBase") #isposdef
Pkg.add("LinearAlgebra") #eigen, Hermitian
Pkg.add("Random") #Random.seed!()
Pkg.add("Combinatorics")
# using Base.Threads #@threads

#a. 3_*jl series (1d local migration model) generate result(time to fixation
 # and the asymptotic number of origins) as half as the result of Joejyn's,
  # cause not found (2020/4/13)

 #speculation: you may try rand(Poisson(lambda)) rather than pois_rand(lambda)
  #in sampling for mutation or genetic drift in function

using PoissonRandom
# Random.seed!(123);rand(pois_rand(10))
Random.seed!(123);pois_rand(10)
 #pois_rand() in PoissonRandom.jl reported to be faster than Poisson()
  #in Distributions.jl; however, this two don't generate the same output
   #while setting the sampling pool(Random.seed!())the same.

using Distributions
Random.seed!(123);rand(Poisson(10))

#b.
 #2_2gaussian_multinomial.jl (positive definite matrix)
 #2_2gaussian_multinomial2.jl (positive semi-definite matrix)

#c.out of memory problem when deme number>500
 #sparse matrix only reduced the memory to 3/4 as observed--insufficient

#d. in data analysis, it could be better to involve
 #Approximate Bayesian computation and carry out systematic simulation
  #rather than only extract certain combination of parameters in interest

#e. rounding issue when scaling the individuals(n) into frequency(x)
 #x=n/sum(n)
 #check num2strexp.m(Matlab codes not translated)


















# indexing D cell array
?

median and mean
guassian and multinomial (set the seed)

for i in 1 #when 1 instead of form 1:2 @threads cannot be used
    println("ok")
  end

##Can D[i][:,t]=x_aft store back correctly?
D=Array{Array{Float64,2},1}(undef,3)#deme cell array - each with their own frequency array
D[1]=rand(2,3)
D[1][2,:]=zeros(1,3)#works
D[1][:,3]=zeros(1,2)#works
D[1]

##nope, plot() cannot creat a new pane for later!
# using Plots.PlotMeasures #for mm
# for i in 1:length(legend)
#     plot(lab="")
#       # if i==1
#           scatter!(tsamp, avn_mtx[:,i], markersize=4, label=legend[i],
#          markershape=:rect, legend=:outertopright,
#            legendfontsize=8, grid=false,fg_legend = :transparent,
#            margin=5mm, size=(800,800), yerror=se_mtx[:,i]) # color="$(colors[i])",
#       # else
#       #     scatter!(tsamp, avn_mtx[:,i], label=legend[i],
#       #   yerror=se_mtx[:,i]) #  color="$(colors[i])",
#       # end
#       plot!(tsamp, avn_mtx[:,i],lab="")#color="$(colors[i])"
#  end
# xlabel!("Generation(t)")
# ylabel!("Average number of origins 𝛍(t) ")
