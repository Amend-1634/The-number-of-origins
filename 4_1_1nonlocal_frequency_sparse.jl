using PoissonRandom
# using LazyArrays #ApplyArray(vcat,) sucks
using SparseArrays
using Multibreak
using Base.Threads
using Plots
include("2_2gaussian_multinomial.jl")

#constant population (and deme) size
#only difference froom "4_1_1nonlocal_frequency.jl" is
 #sparse array(only store non-zero element as entries) used to reduced the
  #memory allocation=>serve for large dataset

# function gaussian_1d_nonlocalmig(NDeme,twoNmu,twoNs,twoNm,T,nDemes, Gauss)
function gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes, Gauss,nu)#N

    #Output
     #D: cell array containing frequency array in each deme -
      #rows of this array correspond to independent origins/haplotypes,
       #with common ordering across all demes (same num of rows)
    #Input:
    # NDeme: population size of each deme
    # nDemes: the number of demes
    # N: population size of metapopulation

    NDeme=round(Int,N/nDemes) #initial value as integer

    # N_vec=fill(NDeme, nDemes) #vector vector containing NDeme corresponding to each deme
    # T: number of generation
    # twoNmu, twoNs & twoNm: the population scaled by the total population.
    # Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)
    # nu: 𝝂 is the degree of freedom of the student t distribution deciding dispersal of individuals

    ## Population Parameters
    mu = twoNmu/2/NDeme/nDemes # mutation coefficient
     #twoNmu: the meta-population scale mutation cases
    #println("mig  ",mig)
     #each deme instead of all demes
    #we kept twoNDememig constant in input so mig change with NDeme
    #same meta-population size and more deme=>more migrations=>more spatial structure


    # s0 = twoNs/2/NDeme/nDemes #selection coefficient, s_ to avoid name crash with s in Julia
    # mig = twoNmig/2/NDeme/nDemes # dispersal rate


    ## Preallocating dataframes
    migCell=Array{Array{Float64,1},1}(undef,nDemes) # stores how many migrants per deme
                               #one column is enough
    migLR=Array{Array{Float64,2},1}(undef,nDemes) #left right migration amount
    mutantCount = spzeros(nDemes,T+1) # counts number of de novo mutants generated
     #per deme for each time
    migrantCount = spzeros(nDemes,T+1) # counts number of migrants generated per deme
     #for each time
    ## Creating the initial frequency array, with each array representing one deme
    # D = cell(1,nDemes);
    D=Array{Array{Float64,2},1}(undef,nDemes)#deme cell array - each with their own frequency array
    X_D=spzeros(nDemes,T+1) #the matrix for the total frequency of all mutants
    X_D[:,1].=0 #no mutant in the first generation

    # X_D_deter= spzeros(1,T+1) #total mutant frequency predicted by the semi-Deterministic algorithm
    # X_D_deter[:,1].=0 #no predicted mutant in the first generation

    x = spzeros(1,T+1)
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
        #println("xwt  ",xwt)
        X = x[1:end-1,t-1] # mutant frequencies at t-1
        #println("X_size  ",size(X))

        # xx = X + xwt*s_*X # new mutant frequencies after selection coefficient
         #simplified to ss=xwt*s_ #yes while rounding issue cause it "not equal"
        # ss=s_[1:end-1] .- s_'*x[:,t-1] #ss: adjustment parameter for probability
        #K-1*1            1*K   K*1
        #K-1*1                1(constant now)
        ss=xwt*s0 #selection modification
        #println("ss ",ss, " ", typeof(ss))
        #println("X", size(X))

        xx=X+ss*X
            #constant*vector
            #println("xx",xx)


        # if isempty(xx)
        #     prinln("xx is empty")
        # end

           #s_=[fill(s0,K-1),0] #0 for wild type

        #in case ss<0 => xx<0
         #but xwt*s_ will never be negative?=>is this still necessary??
        xx[xx.<0].=0
        xxK = 1-sum(xx) #Frequency of wild type (wt is the Kth allele)
        # #println("xxK",xxK)

        if xxK < 0 # in case wild type frequency is negative
            xxK = 0 # set wt frequency to 0
            xx = xx/sum(xx) #resets all mutant frequencies to sum to 1
        end

        ## Genetic drift （Gaussian approximation of multinomial sampling）
        #println("p for drift", [xx;xxK])
        if Gauss == 1 #true
            z=gaussian_multinomial(NDeme,[xx;xxK]) #z in frequency
            n_=round.(Int,z*NDeme)#n in individuals, in integer individuals

            # if sum(n_)!=NDeme
            #     #println("rounding in converting frequency to individuals", sum(n_))
            # end
            #within range (-2,2) so not taking any treatment (2020/02/27)

            # #println("gaussian_sum(z)  ",sum(z))#yes,1.0
            # #println("z",z)
        elseif Gauss == 0 #false
            n_ = rand(Multinomial(Int(NDeme),[xx;xxK])) #Multinomial random numbers
            # #println("multinomial_sum_z  ",sum(z))
        end

        x[:,t] = n_ #n outputed in individuals in the column of this generation
        n=x #represent matrix with only using columnn in individuals
        D[i]=n
        #println("n_ after drift", n_)

        #a zoom/correction of the frequency needed or not??


        # Mutation
        # #println("NDeme  ",NDeme," mu ",mu," xwt ",xwt)
        nwt = n[end,t] #wild type individuals in t generation
        m = pois_rand(mu*nwt) #number of de novo mutants generated
        # #println("mutation",m)
        mutantCount[i,t] = m #stores number of mutants generated into the array
         #i deme, t generation

        #selection, genetic drift then mutation???

        if m != 0 # with >= mutant alleles generated
            if  m > nwt
                m=nwt
            end

            xN = spzeros(m,T+1) #rows of new mutant alleles
            for f in 1:nDemes #add xN to all demes
                x = D[f]
                x=vcat(xN,x) #new mutant piled on top
                x=sparse(x)#transform to sparse matrix to be stored back to D
                if f==i #frequency of the deme that mutant emerged as 1/NDeme
                    x[1:m,t] .=1 # 1 individual arising per mutant seeding into above array
                    x[end,t] = x[end,t] - m
                end
                D[f]=x
            end
        end

    end # for i in 1:nDemes


    ##Migration
     #former: (from t-1 generation =>t)
     #🎫trying: from t to t (migration following the current status of selection, genetic drift and mutation)
     #println("Migration########generation $t")
     migIn=[] #so the arrray counting migration in will renew each generation
                #but persist in this generation to next nDemes

    nHaplotypes=size(D[1])[1] #row number as haplotype number
     migIn=spzeros(nHaplotypes,nDemes) #to be filled by jumps in student t distributions

    ## Determining the number of migrants leaving the demes
    for i in 1:nDemes
        n_i = D[i] #extract from cell array, multi-step=>better computational efficiency
        n_B = n_i[:,t] #frequency column in t generation, i-th deme
        # #println("length(n_B)", length(n_B)) #should be a vector?
        MIG = mig*n_B #expected number of migrants individuals #column vec
        #println("mig ",mig," n_B ",maximum(n_B)," ",minimum(n_B))
        migLV = pois_rand.(MIG) #actual number of migrants
        #println("migLV  ", migLV)

        # for g in 1:nHaplotypes # poisson may generate migration cases > total individuals
        #                         # while having a relatively large dispersal rate
        #
        #     if migLV[g] > n_B[g] # looks for negative values and sets to zero
        #         migLV[g] = Int(n_B[g])
        #     end
        # end

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

                    # #println("i   ",i)
                    # println("check while functional or not##")
                    # println("jump[j]  ", jump[j])
                    # println("to left _can't large than ",i-1)
                    # println("to right   ",nDemes-i)
                    col=Int(i+jump[j])

                    # println("migIn  ",size(migIn))
                    # println("migIn[id_h,col]_before ",migIn[id_h,col])
                    migIn[id_h,col] += 1 #add one individual to the
                    # println("migIn[id_h,col]_after ",migIn[id_h,col])

                        #corresponding haplotype(id_h) and nDemes(i) column

                end #for j in 1:length(jump)

            end #for id_h

        # println("migLV before subtract  ", unique(migLV))
        n_lv = n_B - migLV #after leaving; both column vector
        n_i[:,t] = n_lv # storing population after migrants leave in
        #println("migration_leaving_size(n_lv)",size(n_lv))
        D[i] = n_i # stores "frequency" array back in to deme cell array

        migCell[i] = migLV # leaving array! not left-right array (migLR)
        # println("migLV  ", migLV)
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

        # println("n_B before +migIn",n_B)
        # if unique(migIn[:,i])!=[0]
            # println("migIn_unique", unique(migIn[:,i]))
            # println("migIn", migIn[:,i])
        # end

        n_aft=n_B + migIn[:,i]
        # println("n_aft after +migIn",n_aft)

        x_aft = n_aft/sum(n_aft)
        D[i][:,t]=x_aft #store vector of individuals back to D cell array
        #println("x_aft_end ", x_aft)
        #println("D[i][:,t] ", D[i][:,t])
    end

  # non-spatial
     # ##mean frequency among all demes
     #  #a. semi-deterministic prediction of the total mutant frequency
     #  new_Xd=((s0*0+mu)*exp(1)^((s0+mu)*t)-mu*(1-0))/
     #  ((s0*0+mu)*exp(1)^((s0+mu)*t)+s0*(1-0)) #0 here as the total frequency of pre-existig mutants
     #  X_D_deter[t]=new_Xd

  end #for t = 2:T+1


  #b. total frequency of all mutants
  for i in 1:nDemes
      X_D[i,:]=sum(D[i][1:end-1,:], dims=1) #store total frequency of mutants
  end
  X_D=mean(X_D,dims=1) #now the mean value of all the demes

 # #c. for each haplotype---non-spatial
  # D_alldemes=D[1]
  # for i in 2:nDemes
  #     D_alldemes=D[i]+D_alldemes
  # end
  # D_alldemes=D_alldemes/nDemes

  return D, X_D, migrantCount#, X_D_deter, D_alldemes
 end #function Gaussian_1d_localmig()
