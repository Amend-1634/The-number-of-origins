#why my fig2_nonspatial_origins plot have higher plateau than 2019paper?
using Distributions
using Base.Threads

print(rand(Multinomial(100,[0.1,0.3,0.6])))

function simu_p(n)
    x=rand(n)
    x=x/sum(x)
    return x
end

##dif k (the number of categories of probability)
function multi(k)
 a=Array{Int}(undef,k,0)
 p=simu_p(k)
 for i in 1:1000
  b=rand(Multinomial(100,p))
  a=hcat(a,b)
 end
 return a
end

using Plots
# plot(a)

a=multi(4)#change input and see
histogram(a'/100,xticks=0:0.05:0.8,bar_width=0.005)
 #visually normal distribution

##dif @thread
# using Base.Threads
function multi(k)#k as the number of categories
 a=Array{Int}(undef,k,0)
 p=simu_p(k)
 @threads for i in 1:1000
  b=rand(Multinomial(100,p))
  a=hcat(a,b)
 end
 return a
end

a=multi(4)
histogram(a'/100,xticks=0:0.05:0.8,bar_width=0.005)
 #samplign seed change every time
pwd()
savefig("..\\photo\\multinomial_distribution_test")

#conclusion:
 #seed (sampling pool) change with replication, k. @threads condition too
 #visually normal distribution
