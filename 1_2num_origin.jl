#time series sampling

using PoissonRandom
using Distributions #Multinomial()
using Random
using Combinatorics
using Base.Threads #@threads
using Plots
pwd()
cd()
include("")


function num_origin(s0_seq,twoNmu_seq,N,replicates) #subplot generator
#s0_seq:vector of selection coefficient(s0)--enable comparison
#twoNmu_seq:vector of population-scaled mutation rate (twoNmu)--enable comparison
#N: total population size
#replicates

  ndx_mtx= Array{Float64}(undef,2,0) #collect the paired s0 and twoNmu for later legend

  simu_mtx= Array{Float64}(undef,2001,0) #matrix for simulated average number of origins
   #2001:T+1;
  simu_mtx_se= Array{Float64}(undef,2001,0)#matrix for the standard error of
   #simulated average number of origins

  semi_mtx= Array{Float64}(undef,2001,0) #matrix for predicted average number of origins
  semi_mtx_se= Array{Float64}(undef,2001,0) #matrix for the se of predicted number of origins

 @threads for s0 in s0_seq #@threads could make it faster by apply multi-threads
    #but be careful if using @threads: a. s0_seq has to have >1 elements
                                      #b. the for loop is parallel instead of successive

    @threads for twoNmu  in twoNmu_seq

      simu_seq_pile=Array{Float64}(undef,2001,0)
      semi_seq_pile=Array{Float64}(undef,2001,0)
      @threads for j in 1:replicates
         K=1
         T=2000

         x, X, X_deter, K=allele_frequency(K,N,twoNmu,s0,T)#call the function
         #simu_seq  =concatenate=>  simu_seq_pile  =summarize replicates=>  simu_seq_mean
          #and simu_seq_se  =concatenate=>  simu_mtx and simu_mtx_se
         simu_seq=Float64[]
         semi_seq=Float64[]

          for k in 1:size(x)[2]#size(x)[2]: generations
             ns=Int(1e3)#sample size
             #Simulations
             xk = rand(Multinomial(Int(ns),x[:,k])) #sampling noise-same distribution with genetic drift
             simu=count(a->0,xk[1:end-1]) #sampled allele--simulation
             #threshold: when mutant allele frequency >= 1/ns can be sampled

             push!(simu_seq,simu)

             semi=twoNmu*log(1+X[k]*ns/twoNmu)
             push!(semi_seq,semi)

          end #for k in generations
         simu_seq_pile=hcat(simu_seq_pile,simu_seq)#row as generation, col as replicate
         semi_seq_pile=hcat(semi_seq_pile,semi_seq)

     end#for j in 1:replicates
     simu_seq_mean=sum(simu_seq_pile,dims=2)/replicates #dims=2: by row;
     simu_seq_se=std(simu_seq_pile,dims=2)/sqrt(replicates) #standard error

     semi_seq_mean=sum(semi_seq_pile,dims=2)/replicates
     semi_seq_se=std(semi_seq_pile,dims=2)/sqrt(replicates)

     ndx_mtx =hcat(ndx_mtx,[s0,twoNmu]) #index matrix, include s0 and twoNmu

     simu_mtx=hcat(simu_mtx, simu_seq_mean) #matrix containing time-series average number of origins from simulations
     simu_mtx_se=hcat(simu_mtx_se, simu_seq_se) #matrix containing the se of time-series number of origins from simulations

     semi_mtx=hcat(semi_mtx,semi_seq_mean) #matrix containing time-series average number of origins from semideterministic prediction
     semi_mtx_se=hcat(semi_mtx_se,semi_seq_se) #matrix containing the se of time-series number of origins from semideterministic prediction

    end#for s0
 end#for twoNmu
 return ndx_mtx,simu_mtx,simu_mtx_se,semi_mtx,semi_mtx_se
end#function num_origin()
