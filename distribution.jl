using Distributions
 a=rand(TDist(30),100)
 histogram(a,xlims=(-10,10),label="30")

 b=rand(TDist(20),100)
 histogram!(b,alpha=0.5)


a=[]
 using Base.Threads
 @threads for i in 1:3
     a[i]=rand(1)
 end
a


a=zeros(3,4)
 @threads for i in 1:3, j in 1:4 #not supported
     a[i,j]=rand(1)
 end
 a

a=zeros(3,4)
  @threads for i in 1:3 #ok
      for j in 1:4
          a[i,j]=rand(1)[1]
      end
  end
  a


A = Array{Array{Float64,2},1}(undef,10)#nDemes=10

using Distributions
    function b()
        jump=rand(TDist(1),10000)#12>10
        # jump=trunc.(jump)
        for j in 1:length(jump)
         # while jump[j]<0 && abs(jump[j]) >= 3 || jump[j]>0 && jump[j] > 7 || abs(jump[j])<1 #jump[j]==0
         while abs(jump[j])<1 #jump[j]==0
             jump[j]=4*rand(TDist(1),1)[1]
             # jump=trunc.(jump)
         end
        end
            return jump
    end
    jump=b()
    jump2=rand(TDist(1),10000)#nature t
    # jump3=rand(TDist(30),10000)#approximate Gaussian distribution
    jump4=rand(Normal(0,1),10000)#approximate Gaussian distribution
    jump5=rand(Binomial(1,0.5),10000)#one number...cannot compare

jump

# using Gadfly #disaster!!!!!!keep away
# plot(jump,bincount=10000,Geom.histogram)
using Plots
gr()
plot(jump)

# jump2=rand(TDist(1),10000)#𝝂=1 => extreme--far different from Gaussian
#                             #very centralized=>check with binomial
# histogram(jump2,normalize=:probability)

using Plots
using Distributions

X = rand(Normal(), 100)
histogram(X, bins=5)

using StatsPlots
    histogram(jump,alpha=0.3,xlims=(-25,25),label="limited t",normalize=:probability,bins=10000)
    histogram!(jump2,alpha=0.8,label="nature t",xlims=(-25,25),normalize=:probability,bins=10000)
    histogram!(jump4,alpha=0.5,label="Gaussian",normalize=:probability,bins=10000)
  # histogram!(jump3,alpha=0.7,label="Approx Gaussian by t",normalize=:probability)
  # histogram!(jump5,alpha=0.3,label="Binomial",normalize=:probability)

using PyPlot
PyPlot.plt[:hist]([jump],bins=8)#???are you ok

plot(jump, x="Price", Geom.histogram)#dataset("ggplot2", "diamonds")

histogram(rand(TDist(4),10000),xlims=(-15,15),normalize=:probability,alpha=0.2,bar_position=:stack,label="4")
 histogram!(rand(TDist(10),10000),normalize=:probability,alpha=0.5,label="10")
 histogram!(rand(TDist(30),10000),normalize=:probability,alpha=0.8,label="30")

 histogram!(rand(TDist(1),10000),normalize=:probability,alpha=0.8,label="1")#extremly
 #tail much fatter=>wider spread

 #lower 𝝂, more centralized



#𝝂=1,2,30
using Distributions
 histogram(rand(TDist(1),10000),xlims=(-15,15),normalize=:probability,alpha=0.8,bins=100000,label="v=1",grid=false)
 histogram!(rand(TDist(2),10000),normalize=:probability,alpha=0.5,label="v=2")
 histogram!(rand(TDist(30),10000),normalize=:probability,alpha=0.4,label="v=30")


#𝝂=30, original and elongated and reduced
histogram(rand(TDist(1),10000),xlims=(-15,15),normalize=:probability,alpha=0.8,bins=100000,label="Original",grid=false)
 histogram!(5 .* rand(TDist(1),10000),xlims=(-15,15),normalize=:probability,alpha=0.8,bins=800000,label="Elongated",grid=false)









#Centralization:
#fat-tailed=> Gaussian =>
    #if the population range issymetric or narrower
        #redistributed (limited) t distribution

nDemes=[10 50 100]
as nu

using Distributions
cdf(Uniform(0, 1), 0.6) - cdf(Uniform(0, 1), 0.2)
# MvNormal(covariance)#multivariate normal distribution
# cd(MvNormal(..),0.95)

cdf(Normal(0,1),1.29) #from z to the probability area to its left
cdf(Normal(0,1),0.05)





##
using Distributions

using Distributions
    jump2=rand(TDist(1),1000)#12>10
    jump3=rand(TDist(2),1000)#12>10
    jump4=rand(TDist(30),1000)#approximate Gaussian distribution
    jump5=rand(Normal(0,1),1000)#approximate Gaussian distribution

histogram(jump5,xlims=(-15,15), alpha=0.6,label="Normal",normalized=:pdf,bins=50)
    # histogram!(jump3,alpha=0.7,label="v=2",normalized=:pdf,bins=200)
    histogram!(jump2,alpha=0.4,label="v=1",normalized=:pdf,bar_position=:overlay,bins=5000)
    # histogram!(jump4,alpha=0.4,label="v=30",normalized=:pdf,bins=50)
    histogram!(jump,xlims=(-15,15), alpha=0.6,label="trunc",normalized=:pdf)




function x()
    jump=rand(TDist(1),1000)#12>10
    jump=trunc.(jump)
    for j in 1:length(jump)
     # while jump[j]<0 && abs(jump[j]) >= 3 || jump[j]>0 && jump[j] > 7 || jump[j]==0
     while  jump[j]==0
         jump[j]=4*rand(TDist(1),1)[1]
         jump=trunc.(jump)
     end
    end
    return jump
    end
    jump=x()
    jump
    jump2=rand(TDist(1),1000)#12>10  ##nature t dist, v=1
    jump3=rand(Binomial(1),1000)#binomial, local

histogram!(jump2,alpha=0.3,label="v=1",normalized=:pdf,bar_position=:overlay,bins=10000)

histogram(jump3,alpha=0.5,label="Binomial",normalized=:pdf)#cannot be normalized


# import matplotlib.pyplot as plt


using Plots, StatsPlots
histogram(jump,alpha=0.5,xlims=(-15,15),label="limited",normalized=:pdf)

# plot(jump,seriestype=:bar,density=true)

 # bins=:scott
 #,normed=1
 # bins=:scott
 # , weights=repeat(1:5, outer=200)


histogram(jump2,label="unlim",alpha=0.5,xlims=(-15,15))

histogram(jump2,xlims=(-25,25),alpha=1,label="v=1",normalized=:pdf,bar_position=:dodge)
histogram!(jump3,alpha=0.5,label="Approx Gaussian",normalized=:pdf)
histogram!(jump4,alpha=0.2,label="Gaussian",normalized=:pdf)
histogram!(jump5,alpha=0.2,label="Gaussian",normalized=:pdf)



sum(jump)
sum(jump2)

mean(jump)
mean(jump2)


##power law Paulose, 2019 used
#total populatio N
#N7   #N9
tot_freq
ave_origins

a=collect(1:24)
c=float.(-1 .- a)
b=a .* (3 .^ c)
using Plots
plot(b)

# J_((r))=μ〖 r〗^(-(1+μ))

3^(-4)
