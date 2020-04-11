using Statistics
using Multibreak
include("3_1_2notrack_cons.jl")

#1d local migration model
#sampling of num of origins along the time points
#constant population size

function num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
    Nsamp = 1000 #total sample size of meta-population
    t_fix_vec=[] #matrix
    eta_fix_vec=[]
    fix_M=0 #how many cases reaching fixation in replicates

    @multibreak for m in 1:M ##Replicates
     #@multibreak used for line40, 43(more than one break or continue)
        # println("replicate#######",m) #Showing the progress of the simulation

        D = gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes,Gauss)
         for t in 1:T+1 #generation
            for i in 1:nDemes
                 n_i=D[i]
                 xwt=n_i[end,t]
                 if xwt==0
                     if i!=nDemes
                         continue
                     elseif i==nDemes
                         # println("reached fixation#####  t  ",t)
                         fix_M+=1 #+1 for counting fixation cases when reachign fixation
                         push!(t_fix_vec,t)

                         #the number of origins at fixation time point
                         n_mu=size(D[1])[1]-1 #total mutant haplotype amount
                         x=zeros(n_mu,1)
                         for i in 1:nDemes
                             x_=D[i][1:end-1,t]
                             x = x + x_
                         end
                         eta_fix=count(a->a>0,x)
                         push!(eta_fix_vec, eta_fix)

                         break; break; continue # to start next replication
                     end
                 else
                     break; continue #break nDemes loop and skip this generation
                 end
                 # println(" m ", m, " t ", t," i ", i )
                  #the multibreak above functions if all those not printed out
            end #for i
            # println("m ", m, " t ", t)
        end #for t
        # println("m ", m)

     end #for m

     t_avn=mean(t_fix_vec)
     t_se=std(t_fix_vec)/sqrt(M) #standard deviation

     eta_avn=mean(eta_fix_vec)
     eta_se=std(eta_fix_vec)/sqrt(M)

     return t_avn, t_se, eta_avn, eta_se, fix_M
end #function num_origin_1d_fix()
