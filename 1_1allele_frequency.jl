using PoissonRandom
using Distributions
using Random
using LazyArrays

function allele_frequency(K,N,twoNmu,s0,T)
  # K - initial number of alleles (K-1 resistance alleles and 1 "wild type")
  # N - population size
  # twoNmu - population-scaled mutation rate
  # s0 - selection coefficients
  # T - number of generations

  n0=zeros(K-1,1)
  #if assuming no pre-existing allele at the beginning K=1
  #actually K is not nesccessary here unless we input the frequency of
   #pre-exsiting allele instead of using zero matrix

  mu=twoNmu/2/N #mutation rate per population per individual (2 here as haploid)
  #s0=twoNs/2/N #selection strength
  s_=ApplyArray(vcat,fill(s0,K-1),0) #vector of selection coefficients relative to the Kth allele
  x=zeros(K,T+1) #x: frequency (0 to 1)
  t=collect(1:T+1) #time vector
  n0=ApplyArray(vcat,n0,[N-sum(n0)])#N-sum(n0): wild type individuals
  x0=n0/N #scale individuals (n0) to frequency (x0)
  x[:,1]=x0 #define the frequency column of first generation
  xwt=x[end,1] #extract wild type frequency
  X_deter=Float64[]#For later collecting Deterministic growth of total frequency of all mutant alleles

  q=0 #pre-existing origins (mutants) number (to update m)
  for i in 2:T+1 #start from 2nd generation end in T+1th generation
    xwt=x[end,i-1]
    #Random.seed!(1234) #set seed to define the sampling pool if you want the sampling result the same
    m=pois_rand(N*mu*xwt) #mutant number arised in population this generation
    #N*mu*xwt  as parameter

      if q+m>K-1
         x_1=x[1:end-1,:] #bottleneck while profiling--the codes will be optimized in 2_*.jl or 3_*.jl...series
         x_2=zeros(q+m-K+1,T+1)
         x_3=collect(x[end,:]') #collect() transforms LinearAlgebra.Adjoint=>normal array
         x=ApplyArray(vcat,x_1,x_2,x_3)#optimizedd vcat(), vertically concatenate
         #better than x=[x1;x2;x3] #bottleneck

         s_=ApplyArray(vcat,s_[1:end-1],fill(s0,q+m-K+1),0)
          #amend selection coefficient vector correspondingly

         K=q+m+1 #q+m is the predicted total mutant in Poisson distribution
       end

     #vector to add mutant frequency (1/N) to frequency matrix
     mx=zeros(K,1)
      if m!=0
         mx[q+1:q+m]=fill(1/N,m)
         #each mutant must arise in population at frequency 1/N
         q=q+m
      end

     # ss=s_[1:end-1] .- s_'*x[:,i-1] #ss: adjustment parameter for probability
     #size K-1*1            1*K   K*1
     #size K-1*1                1(constant now)
     #then simplified to xwt*s0
     ss=xwt*s0

     X=x[1:K-1,i-1]#all mutants frequency in last generation
     xx=X+ss.*X#.* means corresponding multiplication
     #matrix multiplication: ss as vector(n,)*(n,1), same to 2d (1,n)*(n*1)

     #xx: probability of mutant allele in this Generation
     xx[xx.<0].=0
     xxK=1-sum(xx)

     #for wild type and correction (if not)
      if xxK<0
         xxK=0
         xx=xx/sum(xx)
      end

    #Genetic drift
    n = rand(Multinomial(Int(N),[xx;xxK])) #Multinomial random numbers
    #N as sampling times, [xx,xxK] probability of mutant and wild type
    #xx+xxK=1

    z=n/N+mx #z include all alleles' frequency
    #add mutant frequency (1/N) arised this generation to the frequeny matrix
     #after genetic drift

    if sum(z)>1
        z=z/sum(z) #scale
    end
    x[:,i]=z
    if isnan(sum(x[:,i]))
            print("Error")
            println(x[:,(i-10:i)])#10 generations backward to see what happened
            println(n)
            println([xx,xxK])
            break #end the for loop
    end

    #Deterministic growth--prediction of X
    X=[0]#no mutant allele assumed at the beginning
    new_Xd=((s0*X[1]+mu)*exp(1)^((s0+mu)*i)-mu*(1-X[1]))/
     ((s0*X[1]+mu)*exp(1)^((s0+mu)*i)+s0*(1-X[1])) #estimated total mutant frequency in Deterministic growth
     push!(X_deter,new_Xd) #add the prediction in each generation to the vector of Deterministic prediction vector

   end #end of for i in 2:T+1

   prepend!(X_deter,[0]) #prepend the Deterministic total mutant frequency in the first generation
   X=sum(x[1:end-1,:],dims=1)#column sum-total mutant frequency each generation
    #simulated total mutant frequency
   K=size(x)[1] #haplotype number
   return x, X, X_deter, K

   #x-allele frequency
   #X-simulated total mutant frequency
   #X_deter-predicted total mutant frequency according to deterministic growth
   #K-total mutant numbers (could be lost during drift)

 end #function allele_frequency()
