using PoissonRandom
using LazyArrays
using SparseArrays
using Multibreak
using Base.Threads
using Plots
include("2_2gaussian_multinomial.jl")

#scaling without traking the migrant individuals
#constant population size

# function gaussian_1d_localmig(NDeme,twoNmu,twoNs,twoNm,T,nDemes, Gauss)
function gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)#N

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

    # nDemes = 10 # deme number
    # twoNmu = 10;
    # s0 = 0.05 # selection coefficient

    ## Population Parameters
    mu = twoNmu/2/NDeme/nDemes # mutation coefficient
    # s0 = twoNs/2/NDeme/nDemes #selection coefficient, s_ to avoid name crash with s in Julia
    # mig = twoNm/2/NDeme/nDemes # dispersal rate
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
    # D=Array{Array{Float64,2},1}(undef,nDemes)#deme cell array - each with their own frequency array
    D=Array{SparseMatrixCSC{Float64},1}(undef,nDemes)

    x = spzeros(1,T+1)
    x[end,1] = 1 #Initial frequency array of the wt (no mutants present yet)

    for i in 1:nDemes # populate deme cell array with initial wt-frequency of 1 @ T= 0
        D[i] = x
    end

 for t in 2:T+1 #Generations #T+1
    # println("generation######################",t)
     for i in 1:nDemes
        # println("deme################", i)

        ##Selection
        ## Filling up the frequency array
        x=D[i]
        ##mutation, selection, genetic drift
        xwt = x[end,t-1] #wild type frequency at t-1
        # println("xwt",xwt)
        X = x[1:end-1,t-1] # mutant frequencies at t-1
        # println("X_size",size(X))

        # xx = X + xwt*s_*X # new mutant frequencies after selection coefficient
         #simplified to ss=xwt*s_ #yes while rounding issue cause it "not equal"
        # ss=s_[1:end-1] .- s_'*x[:,t-1] #ss: adjustment parameter for probability
        #K-1*1            1*K   K*1
        #K-1*1                1(constant now)
        ss=xwt*s0 #selection modification
        # println("ss ",ss, " ", typeof(ss))
        # println("X", size(X))

        xx=X+ss*X
            #constant*vector
            # println("xx",xx)


        # if isempty(xx)
        #     prinln("xx is empty")
        # end

           #s_=[fill(s0,K-1),0] #0 for wild type

        #in case ss<0 => xx<0
         #but xwt*s_ will never be negative?=>is this still necessary??
        xx[xx.<0].=0
        xxK = 1-sum(xx) #Frequency of wild type (wt is the Kth allele)
        # println("xxK",xxK)

        if xxK < 0 # in case wild type frequency is negative
            xxK = 0 # set wt frequency to 0
            xx = xx/sum(xx) #resets all mutant frequencies to sum to 1
        end

        ## Genetic drift （Gaussian approximation of multinomial sampling）
        if Gauss == 1 #true
            z=gaussian_multinomial(NDeme,[xx;xxK]) #z in frequency
            n_=round.(Int,z*NDeme)#n in individuals, in integer individuals

            # if sum(n_)!=NDeme
            #     println("rounding in converting frequency to individuals", sum(n_))
            # end
            #within range (-2,2) so not taking any treatment (2020/02/27)

            # println("gaussian_sum(z)  ",sum(z))#yes,1.0
            # println("z",z)
        elseif Gauss == 0 #false
            n_ = rand(Multinomial(Int(NDeme),[xx;xxK])) #Multinomial random numbers
            # println("multinomial_sum_z  ",sum(z))
        end

        x[:,t] = n_ #n outputed in individuals in the column of this generation
        n=x #represent matrix with only using columnn in individuals
        D[i]=n
        # println("z_after drift", z)

        #a zoom/correction of the frequency needed or not??


        # Mutation
        # println("NDeme  ",NDeme," mu ",mu," xwt ",xwt)
        nwt = n[end,t] #wild type individuals in t generation
        m = pois_rand(mu*nwt) #number of de novo mutants generated
        # println("mutation",m)
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
                x=sparse(x)
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
     #trying: from t to t (migration following the current status of selection, genetic drift and mutation)
    ## Determining the number of migrants leaving the demes
    for i in 1:nDemes
        n_i = D[i] #extract from cell array, multi-step=>better computational efficiency
        n_B = n_i[:,t] #frequency column in t generation, i-th deme
        # println("length(n_B)", length(n_B)) #should be a vector?
        MIG = mig*n_B #expected number of migrants individuals #column vec
        migLV = pois_rand.(MIG) #actual number of migrants
        #size: (allele_amount,)
        # println("size(migCell[i])",size(migCell[i]))

                          #Before leaving in migration
                          #column vector
        for g in 1:length(n_B) # poisson may generate>there is=>negative frequency
            #"migCell[i] = poissrnd(MIG)' #actual number of migrants"--line111
            if migLV[g] > n_B[g] # looks for negative values and sets to zero
                # println("migLV[g]  ", migLV[g])
                # println("n_B[g]  ", n_B[g])
                migLV[g] = Int(n_B[g]) #both float vector
            end
        end

        n_lv = n_B - migLV #after leaving; both column vector
        n_i[:,t] = n_lv # storing population after migrants leave in
        # println("migration_leaving_size(n_lv)",size(n_lv))
        D[i] = n_i # stores "frequency" array back in to deme cell array

        migCell[i] = migLV # leaving array! not left-right array (migLR)
        # println("migLV  ", migLV)
        migrantCount[i,t] = sum(migLV) #sum of migrants per deme
    end

        ## Determining if migrants move left or right
    for i in 1:nDemes
        c = size(migCell[i])[1] #c=column number, registers size of migrant
        #array to preallocate size of migLR_i array
        migLR_i = Int.(spzeros(c,2)) #column 1 amount going left, 2 going right
        tot_mig = migCell[i] #total: left+right
        # println("tot_mig ", " max " maximum(tot_mig)," type ",typeof(tot_mig))
        tot_mig = Int.(tot_mig)
        if i == 1 #first deme
            # println("first_size",size(tot_mig))
            migLR_i[:,2] = tot_mig # left end, only toward right
            migLR[i] = migLR_i
        elseif i == nDemes #last deme
            # println("last_size", size(tot_mig))
            migLR_i[:,1] = tot_mig # right end, only toward left
            migLR[i] = migLR_i
        else # middle demes
            # migLR_i(1,:) = binornd(tot_mig,0.5);
            # println("middle_size", size(tot_mig))
            # println("mig")
            # println("middle_size",size(tot_mig))
            migLR_i[:,1] = rand.(Binomial.(tot_mig,0.5))#toward left
                            #two elemetwise dot is needed
            migLR_i[:,2] = tot_mig - migLR_i[:,1] #toward right
            migLR[i] = migLR_i
        end
    end
#migLR:    column1           column2
#leaving:  toward left       toward right
#entering: fr Right          fr Left

    ## Accounting for migrants entering in the "frequency" array, D
    for i in 1:nDemes
        if i == 1 # for first deme
            n_i = D[i] # n before accounting for entering migrants
            n_B = n_i[:,t]
            n_frR = migLR[i+1] # n entering array !ROW VECTOR
            n_frR = n_frR[:,1] #only migration from right, columnn vector
            #column vector
            # println("size(n_B)",size(n_B))
            # println("size(n_frR)",size(n_frR))

            n_aft = n_B + n_frR
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
            elseif i == nDemes # for last deme
            n_i = D[i]
            n_B = n_i[:,t]
            n_frL = migLR[i-1] # n entering array !ROW VECTOR
            n_frL = n_frL[:,2] #only migration from left, columnn vector
            # println("size(n_B)",size(n_B))
            # println("size(n_frL)",size(n_frL))
            n_aft = n_B + n_frL # adds entering array into existing array
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
        else # for the middle demes
            n_i = D[i] # n before entering, COLUMN VECTOR
            n_B = n_i[:,t]
            # println("migLR_size",size(migLR))
            n_frR = migLR[i+1] # n entering array from right deme !ROW VECTOR
            n_frR = n_frR[:,1] # n entering (only left migration therefore                #row 1)!COLUMN VECTOR
            n_frL = migLR[i-1] # n entering array from left deme !ROW VECTOR
            n_frL = n_frL[:,2] # n entering (only right migration therefore
                #row 2)!COLUMN VECTOR

                 # println("size(n_B)",size(n_B),
                 #         "size(n_frL)",size(n_frL),
                 #         "size(n_frR))",size(n_frR))

            n_aft = n_B + n_frR + n_frL # adds entering array into existing array
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
        end
    end

    #convert individual to frequency
    for i in 1:nDemes
       #not tracking migration individuals
        n_B=D[i][:,t]#access matrix in individuals
        # println("n_B",n_B)
        x_aft=n_B/sum(n_B) #transform into frequency matrix
        D[i][:,t]=x_aft #store frequency matrix back to D cell array
        # println("x_after_ind=>freq",D[i][:,t])
     end


  end #for t = 2:T+1

  # Setting all frequencies to sum up to 1 after migration (convert
   #from population into frequency)
  for i in 1:nDemes
      x_i = D[i] # access frequency matrix of i-th deme
      for t in 1:T+1
          x_B = x_i[:,t]

          if sum(x_B)!=1
              # println("before scale at end, ", sum(x_B))
              # println("x_B", x_B)
              x_B=round.(x_B, digits=10)
              #could be due to rounding in line 266
              x_aft = x_B/sum(x_B)
              # println("after scale at end, ", sum(x_aft))

              x_i[:,t]=x_aft #store frequency vector back
              # println("x_after_final scaling",x_aft)

          end
      end
      D[i] = x_i # stores back into cell array
  end

  return D
 end #function Gaussian_1d_localmig()
