using PoissonRandom
using LazyArrays
using Multibreak
using Base.Threads
using Plots
include("2_2gaussian_multinomial.jl")

#scaling without tracking the migrant individuals
 # =>constant population size

#1d model: migration distribution as Student's t distribution
 # =>include local and non-local migration

function gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes, Gauss,nu)#N

    #Output
     #D: cell array containing frequency array(x) in each deme -
      #rows of x--independent origins/haplotypes,
       #with common ordering across all demes (same num of rows)
      #column of x--generations
     #X_D:  cell array containing total mutant frequency array(x)
      #in each deme, row and col of x setting same as above

    #Input:
    # N: population size of metapopulation
    # twoNmu: the population-scaled mutation rate
    # s0 = twoNs/2/NDeme/nDemes #selection coefficient, s_ to avoid name crash with s in Julia
    # mig = twoNmig/2/NDeme/nDemes # migration rate
    # T: number of generation
    # nDemes: the number of demes
    # Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)
    # nu: 𝝂 is the degree of freedom of the student t distribution deciding migration of individuals

    NDeme=round(Int,N/nDemes) #initial value as integer
    # NDeme: population size of each deme

    ## Population Parameters
    mu = twoNmu/2/NDeme/nDemes # mutation coefficient

    ## Preallocating dataframes
    migCell=Array{Array{Float64,1},1}(undef,nDemes) # stores how many migrants per deme
                               #one column is enough
    migLR=Array{Array{Float64,2},1}(undef,nDemes) #left right migration amount
    mutantCount = zeros(nDemes,T+1) # counts number of de novo mutants generated
     #per deme for each time
    migrantCount = zeros(nDemes,T+1) # counts number of migrants generated per deme
     #for each time
    ## Creating the initial frequency array, with each array representing one deme
    D=Array{Array{Float64,2},1}(undef,nDemes)#deme cell array - each with their own frequency array
    X_D=zeros(nDemes,T+1) #the matrix for the total frequency of all mutants
    X_D[:,1].=0 #no mutant in the first generation

    # X_D_deter= zeros(1,T+1) #total mutant frequency predicted by the semi-Deterministic algorithm
    # X_D_deter[:,1].=0 #no predicted mutant in the first generation

    x = zeros(1,T+1)
    x[end,1] = 1 #Initial frequency array of the wt (no mutants present yet)

    for i in 1:nDemes # populate deme cell array with initial wt-frequency of 1 @ T= 0
        D[i] = x
    end

 for t in 2:T+1 #Generations #T+1
    #println("generation######################",t)
     for i in 1:nDemes
        #println("deme################", i)

        ##Selection
        ## Filling up the frequency array
        x=D[i]
        ##mutation, selection, genetic drift
        xwt = x[end,t-1] #wild type frequency at t-1
        X = x[1:end-1,t-1] # mutant frequencies at t-1

        # ss=s_[1:end-1] .- s_'*x[:,t-1] #ss: adjustment parameter for
         #probability simplified to ss=xwt*s_
         #only rounding issue cause original or simplified appears to be
          #"not equal"
        #xx = X + xwt*s_*X # new mutant frequencies after selection coefficient
        #(K-1,1)  1*(1,K)*(K,1)
        #(K-1,1)    1(constant now)
        #(a,b) is describing the size of matrix (a row, b column)
        ss=xwt*s0 #selection modification

        xx=X+ss*X
          #vector + constant*vector

        #s_=[fill(s0,K-1),0] #0 for wild type #not needed after simplfication

        xx[xx.<0].=0 #in case ss<0 => xx<0
        xxK = 1-sum(xx) #Frequency of wild type (wt is the Kth allele)

        if xxK < 0 # in case wild type frequency is negative
            xxK = 0 # set wt frequency to 0
            xx = xx/sum(xx) #resets all mutant frequencies to sum to 1
        end
        p=[xx;xxK]

        if !isempty(p[p.<0])
            println("tsamp_p<0",p[p.<0])
        end

        if sum(p)!=1
            println("tsamp_sum(p)", sum(p))
        end

        ## Genetic drift （Gaussian approximation of multinomial sampling）
        #println("p for drift", [xx;xxK])
        if Gauss == 1 #true
            # println("[xx;xxK]",[xx;xxK])

            z=gaussian_multinomial(NDeme,[xx;xxK]) #z in frequency
            n_=round.(Int,z*NDeme)#n in individuals, in integer individuals

            # if sum(n_)!=NDeme
            #     #println("rounding in converting frequency to individuals", sum(n_))
            # end
            #within range (-2,2) so not taking any treatment (2020/02/27)

        elseif Gauss == 0 #false
            n_ = rand(Multinomial(Int(NDeme),[xx;xxK])) #Multinomial random numbers
        end

        x[:,t] = n_ #n outputed in individuals in the column of this generation
        n=x #represent matrix with only using columnn in individuals
        D[i]=n

        ## Mutation
        nwt = n[end,t] #wild type individuals in t generation
        m = pois_rand(mu*nwt) #number of de novo mutants generated
        mutantCount[i,t] = m #stores number of mutants generated into the array
         #i deme, t generation


        if m != 0 # with >= mutant alleles generated
            if  m > nwt
                m=nwt
            end

            xN = zeros(m,T+1) #rows of new mutant alleles
            for f in 1:nDemes #add xN to all demes
                x = D[f]
                x=ApplyArray(vcat,xN,x) #new mutant piled on top
                if f==i #frequency of the deme that mutant emerged as 1/NDeme
                    x[1:m,t] .=1 # 1 individual arising per mutant seeding into above array
                    x[end,t] = x[end,t] - m
                end
                D[f]=x
            end
        end

    end # for i in 1:nDemes


    ## Migration
     #println("Migration########generation $t")
     migIn=[] #so the arrray counting migration in will renew each generation
                #but persist in this generation to next nDemes

    nHaplotypes=size(D[1])[1] #row number as haplotype number
     migIn=zeros(nHaplotypes,nDemes) #to be filled by jumps in student t distributions

    ## Determining the number of migrants leaving the demes
    for i in 1:nDemes
        n_i = D[i] #extract from cell array, multi-step=>better computational efficiency
        n_B = n_i[:,t] #frequency column in t generation, i-th deme
        MIG = mig*n_B #expected number of migrants individuals #column vec
        migLV = pois_rand.(MIG) #actual number of migrants

        for g in 1:nHaplotypes # poisson may generate migration cases > total individuals
                                # while having a relatively large migration rate

            if migLV[g] > n_B[g] # looks for negative values and sets to zero
                migLV[g] = Int(n_B[g])
            end
        end

        ##determine direction and how many demes (distance) it jumped across
        idnot0=findall(!isequal(0),migLV)
        for id_h in idnot0 #id_h: index of haplotype with migration cases
                jump=rand(TDist(nu), migLV[id_h]) #\binu 𝝂 as the degree of freedom
                jump=trunc.(jump)#matrix of how many deme each individual of
                 #this haplotype has jumped through
                #- as toward left while + as right
                # println("jump ", jump)

                for j in 1:length(jump)
                    while jump[j]<0 && abs(jump[j]) >= i || jump[j]>0 && jump[j] > (nDemes-i) || jump[j]==0
                        jump_j=rand(TDist(nu), 1)[1]
                        jump[j]=trunc(jump_j)
                    end
                    #finish checking jump elements

                    col=Int(i+jump[j])
                    migIn[id_h,col] += 1 #add one individual to the
                     #corresponding haplotype(id_h) and nDemes(i) column

                end #for j in 1:length(jump)

            end #for id_h

        n_lv = n_B - migLV #after leaving; both column vector
        n_i[:,t] = n_lv # storing population after migrants leave in
        D[i] = n_i # stores "frequency" array back in to deme cell array

        migCell[i] = migLV # leaving array! not left-right array (migLR)
        migrantCount[i,t] = sum(migLV) #sum of migrants per deme
    end #for i in 1:nDemes

    # if sum(migrantCount[:,t])!=sum(migIn) #if migLV==migIn
        # println("migLV and In not equal########################")
        # println("migrantCount[:,t]", migrantCount[:,t])
        # println("sum(migIn)", sum(migIn))
    # end

    # #println("migIn_before entering  ",size(migIn))

    ## Accounting for migrants entering in the "frequency" array, D
    for i in 1:nDemes
        n_B=D[i][:,t]

        n_aft=n_B + migIn[:,i]

        x_aft = n_aft/sum(n_aft)
        D[i][:,t]=x_aft #store vector of individuals back to D cell array
    end

  end #for t = 2:T+1


  #Total frequency of all mutants
  for i in 1:nDemes
      X_D[i,:]=sum(D[i][1:end-1,:], dims=1) #store total frequency of mutants
  end
  X_D=mean(X_D,dims=1) #now the mean value of all the demes

  return D, X_D
  #D: cell array containing frequency array(x) in each deme -
   #rows of x--independent origins/haplotypes,
    #with common ordering across all demes (same num of rows)
   #column of x--generations
  #X_D:  cell array containing total mutant frequency array(x)
   #in each deme, row and col of x setting same as above
 end #function Gaussian_1d_localmig()
