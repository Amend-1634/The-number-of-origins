using PoissonRandom
using Distributions #Multinomial()
using Random
using Combinatorics
# using Base.Threads #@threads
include("2_3allele_frequency_gaussian.jl")

function num_origin_gaussian(s0_seq,twoNmu_seq,N,replicates,T) #subplot generator
 #Input
   #s0_seq:  sequence(I mean vector) of initial selection coefficient in comparison
   #twoNmu_seq: sequence(I mean vector) of population-scaled mutation rate (2Nμ) in comparison
   #N: total population size
   #T: generations

  #Setting undefined arrays with certani size--to collect data within loops
  ndx_mtx= Array{Float64}(undef,2,0)

  #simulated data
  simu_mtx= Array{Float64}(undef,T+1,0) #mean number of origins (η)
  simu_mtx_se= Array{Float64}(undef,T+1,0) #standard error of num of origins
  #structure: simu=>simu_seq=>simu_seq_pile=>simu_mtx(_se)

  #semi-deterministic prediction data
  semi_mtx= Array{Float64}(undef,T+1,0) #mean η
  semi_mtx_se= Array{Float64}(undef,T+1,0) #se


 for s0 in s0_seq
     for twoNmu  in twoNmu_seq

      #array created to contain all replicates
      simu_seq_pile=Array{Float64}(undef,T+1,0)
      semi_seq_pile=Array{Float64}(undef,T+1,0)

      for j in 1:replicates
         K=1 #K: pre-existing allele number(assumes 1 wild type and 0 mutant)

         x=allele_frequency_gaussian(K,N,twoNmu,s0,T)[1]
          #[1]: only require the first output element
          #get allele frequency matrix from included function

         X=sum(x[1:end-1,:],dims=1) #total mutant allele frequency

         simu_seq=Float64[] #vector for established allele in one replicate
                                                                #in simulation
         semi_seq=Float64[] #                                    in prediction

          for k in 1:T+1 #for each generation
             ns=Int(1e3)#sample size

             #Simulations
             xk = rand(Multinomial(Int(ns),x[:,k])) #sampling noise-same
              #distribution with genetic drift
             simu=count(a->a>0,xk[1:end-1]) #origin(allele) number--simulation
             #threshold: when mutant allele frequency >= 1/ns can be sampled

             push!(simu_seq,simu)#add origin number in this generation into
              #the time-series vector
             semi=twoNmu*log(1+X[k]*ns/twoNmu)
             push!(semi_seq,semi)

          end #for k in generations

         #pile up the replicates time-series vectors
         simu_seq_pile=hcat(simu_seq_pile,simu_seq)#row-generation, column-replicate
         semi_seq_pile=hcat(semi_seq_pile,semi_seq)

     end#for j in 1:replicates

     #summurize replicates (by row)
     simu_seq_mean=sum(simu_seq_pile,dims=2)/replicates #mean in simulation
     simu_seq_se=std(simu_seq_pile,dims=2)/sqrt(replicates) #get se

     semi_seq_mean=sum(semi_seq_pile,dims=2)/replicates #mean in semideterministic prediction
     semi_seq_se=std(semi_seq_pile,dims=2)/sqrt(replicates)# se

     ndx_mtx =hcat(ndx_mtx,[s0,twoNmu]) #index matrix

     #horizontally concatenate time series summary vector
      #for each combination of s0 and twoNmu
     simu_mtx=hcat(simu_mtx,simu_seq_mean)
     simu_mtx_se=hcat(simu_mtx_se,simu_seq_se)

     semi_mtx=hcat(semi_mtx,semi_seq_mean)
     semi_mtx_se=hcat(semi_mtx_se,semi_seq_se)

    end#for s0
 end#for twoNmu
 return ndx_mtx,simu_mtx,simu_mtx_se,semi_mtx,semi_mtx_se
 #ndx_mtx: index matrix collecting combinations of s0 and twoNmu
 #simu_mtx: simulated time-series mean of number of origins
 #simu_mtx_se: simulated time-series se of number of origins
 #semi_mtx: predicted time-series mean of number of origins in semideterministic algorithm
 #semi_mtx_se:  predicted time-series se of number of origins in semideterministic algorithm

end#function num_origin_Gaussian()
