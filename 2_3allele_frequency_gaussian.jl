using Plots
using PoissonRandom
using Distributions
using Random
using LazyArrays #Apply(vcat,) #could perform with higher computational efficiency compared to vcat()
include("2_2gaussian_multinomial.jl")

function allele_frequency_gaussian(K,N,twoNmu,s0,T)
  #K: pre-existing allele number (including wild type allele and K-1 mutant allele)
  #N: total population size
  #twoNmu: population-scaled mutation rate (2Nμ)
  #s0: initial selection coefficient
  #T: generations (time)

  n0=zeros(K-1,1) #assuming no pre-existing allele at the beginning
  mu=twoNmu/2/N #mutation rate per population per individual
  #s0=twoNs/2/N #selection #commented out as I change input from twoNs to s0
  s_=ApplyArray(vcat,fill(s0,K-1),0) #vector of selection coefficients relative to each allele
  x=zeros(K,T+1)#K: initial number of each allele; x: frequency of each allele (0 to 1)
  t=collect(1:T+1)
  n0=ApplyArray(vcat,n0,[N-sum(n0)]) #vertically concatenate n0(as mutant
   #allele individuals) and wild type individuals [N-sum(n0)]
  x0=n0/N #individuals scaled into frequency
  x[:,1]=x0 #assign first generation allele frequency (mutant and wild type)
  xwt=x[end,1] #extract wild type frequency
  X_deter=Float64[] #vector prepared for Deterministic growth of total
   #frequency of all mutant alleles

  q=0 #pre-existing origins (mutants) number (for iterative m)
   #with K-1=0 mutant allele assumed pre-existing in default we delete
    #K in input in later 3*.jl and 4*.jl series to reduce the redundancy
  for i in 2:T+1

    xwt=x[end,i-1] #wild type frequency

    #using PoissonRandom
    #using Random
    #Random.seed!(1234) #set seed =>same sampling result
    m=pois_rand(N*mu*xwt) #mutant arised in population this generation
    #N*mu*xwt  as parameter

      if q+m>K-1
        #insert the mutant allele row between the previous mutant allele and wild type allele
         x_1=x[1:end-1,:] #bottleneck
         rows=q+m-K+1
         x_2=zeros(rows,T+1)
         x_3=collect(x[end,:]') #collect() transforms LinearAlgebra.Adjoint=>normal array
         x=ApplyArray(vcat,x_1,x_2,x_3)
         # x=[x1;x2;x3] #profiling bottleneck

         s_=ApplyArray(vcat,s_[1:end-1],fill(s0,q+m-K+1),0)#insert selection
          #coefficient for new mutant alleles
         K=q+m+1# q+m is the  total number of mutant allele
       end

     mx=zeros(K,1)#once mutant has arises it cannot arise by mutation again, and
      # so each allele has a special and independent row
      if m!=0
         mx[q+1:q+m]=fill(1/N,m) #mutant arised this generation
         #each mutant must arise in population at frequency 1/N
         q=q+m #iterate q
      end

    #adjust the selection coefficient according to allele frequency in last generation
     ss_c=s_'*x[:,i-1] #(1,K)*(K,1) => (1,1) constant
     #allele frequency (x) high=>selection advantage small=>no more advantage
     ss=s_[1:end-1] .- ss_c  #ss: adjustment parameter for probability
     #K-1*1                      1(constant now)

     X=x[1:K-1,i-1] #all mutants frequency in last generation
     xx=X+ss.*X #.* means elementwise multiplication
     #xx: probability of allele in this Generation

     xx[xx.<0].=0 #in case ss as negative made xx negative
     xxK=1-sum(xx)#correct wild type frequency (xxK here) correspondingly
      if xxK<0
         xxK=0
         xx=xx/sum(xx)
      end

    #Genetic drift
    #using Distributions
    #z as all mutant allele frequency
    z = gaussian_multinomial(N, [xx; xxK]);#input all all allele(mutants, wide type)

    z = z+mx  #mx as the matrix only contian new mutant allele frequency
    #after drift add the mutant freuqency arised in this generation to
     #the allele frequency matrix
    #as the new mutant will not join the genetic drift this generation (its
     #initial freuqncy has to be 1/N instead of being sampled as other number
      #in genetic drift)

    if sum(z)>1#correction by scaling
        z=z/sum(z)
    end
    x[:,i]=z
    if isnan(sum(x[:,i])) #report if there's NA in x
            print("Error")
            println(x[:,i-10:i])#10 generations backward to see what happened?
            println(n)
            println([xx,xxK])
            return "error"
    end

    #Deterministic growth--prediction of X
    X=[0]#no resistant allele assumed at the beginning(required in algorithm as below)
    new_Xd=((s0*X[1]+mu)*exp(1)^((s0+mu)*i)-mu*(1-X[1]))/
     ((s0*X[1]+mu)*exp(1)^((s0+mu)*i)+s0*(1-X[1]))
     push!(X_deter,new_Xd) #from 2nd generation to T+1 generation

   end #end of for i in 2:T+1
   prepend!(X_deter,[0])
    #prepend total frequency of mutant allele in 1st generation to vector
   X=sum(x[1:end-1,:],dims=1)#column sum-total mutant frequency each generation
   K=size(x)[1]#allele amount; K-1 mutants, 1 wild type
   return x,X,X_deter,K

   #K-final allele numbers
   #x-allele frequency
   #X-sum of mutant allele frequency
   #t-generation
 end #of function
