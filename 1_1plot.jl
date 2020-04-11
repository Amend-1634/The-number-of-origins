#testing whether it function appropriately by plotting
#@time to assess the time needed


pwd()
cd("D://Julia/func")
include("1_1allele_frequency.jl")
@time x1,X1,X_deter1,K1=allele_frequency(1,1e6,1,1/50,2000)
                       #allele_frequency(K,N,twoNmu,s0,T)

@time x10,X10,X_deter10,K10=allele_frequency(1,1e6,10,1/50,2000)

#mu=1
using Plots.PlotMeasures, Plots#for mm
 @time begin
 T=size(x1)[2]-1
 t=collect(0:T)
 p1=plot(t, vec(X1),color=:black,line=:solid,yaxis=:log10,ylims=(1e-6,1),label="sum of mutants",dpi=300,margin=6mm)
 plot!(t, X_deter1,color=:red)
 plot!(t,x1[end,:],color=:black,line=:dash)
 for i in 1:size(x1)[1]-1
   plot!(t,x1[i,:],legend=false)
 end
 n_s=1e3#sample size
 f = Plots.font("DejaVu Sans", 5)
 hline!([1/n_s],label="1/n_s")
 annotate!(2050, 1/n_s, text("1/n_s",f, :left))
 N=1e6;s=1/50
 hline!([1/(2N*s)])#critical line for gene to survive
 annotate!(2050, 1/(2N*s), text("1/(2N*s)",f, :left))
 #plot!(1/(2N*s))
 xlabel!("Generation")
 ylabel!("Frequency")
 end

#mu=10
using Plots.PlotMeasures
 @time begin
 T=size(x10)[2]-1
 t=collect(0:T)
 p10=plot(t, vec(X10),color=:black,line=:solid,yaxis=:log10,ylims=(1e-6,1),label="sum of mutants",dpi=300,margin=6mm)
 plot!(t, X_deter10,color=:red)
 plot!(t,x10[end,:],color=:black,line=:dash)
 for i in 1:size(x10)[1]-1
   plot!(t,x10[i,:],legend=false)
 end
 n_s=1e3#sample size
 hline!([1/n_s],label="1/n_s")
 f = Plots.font("DejaVu Sans", 5)
 annotate!(2050, 1/n_s, text("1/n_s",f, :left))
 N=1e6;s=1/50
 hline!([1/(2N*s)])#critical line for gene to survive
 annotate!(2050, 1/(2N*s), text("1/(2N*s)",f, :left))
 #plot!(1/(2N*s))
 xlabel!("Generation")
 ylabel!("Frequency")
 end

@time plot(p1, p10, layout = 2,title=["twoNmu=1" "twoNmu=10"],legend = false)#plot with mutation coeffient 1 or 10
 #plot!(size=(800,400))#cannot set the backend at one time
savefig("D:\\research_project\\Photo\\mu_1_10")

#Task2: Deterministic growth model line--prediction line in plotting
 #of total frequency of all mutant alleles
T=2000
s0=0.05
mu=10/2/1e6

#computational efficiency
#a.
X_deter=Float64[]
 X=Float64[]
 push!(X_deter,0)
 push!(X,0)

 @time begin
  for i in 2:T+1
   new=((s0*X[1]+mu)*exp(1)^((s0+mu)*i)-mu*(1-X[1]))/
    ((s0*X[1]+mu)*exp(1)^((s0+mu)*i)+s*(1-X[1]))
    push!(X_deter,new)
   end
  end
  plot(layout=2)
  plot!(X_deter,yaxis=:log10,ylims=(1e-6,1),subplot=1)

#b.tanh form #???
X_deter=Float64[]
 X=Float64[]
 push!(X_deter,0)
 push!(X,0)

 @time begin
  for i in 2:T+1
   γ=(s0+mu)/2
   i_star= 1/γ * tanh((s0-mu-2s0*X[1])/(s0+mu))^-1 #i_star as t*
   new=(s0-mu)/(2s0) + γ/s0 * tanh(i-i_star)
   push!(X_deter,new)
   end
  end
 plot!(X_deter,yaxis=:log10,ylims=(1e-6,1),subplot=2)

#tanh more computationally intense and weird plot


#multibreak=> establishment (x[]>1/2Ns) means
 #there are cases like this=>estalishment not equal to survive till mutant fixation
 #x>1/ns; could it be lower than 1/2Ns?

#Task3: figure2: 3 series
 #Semideterministic theory: 0~t_k time window over which mutants can be generated
 #that could contribute to a sample at time T.



 #Task4: using Profile-optimizing the codes
 # @profiler allele_frequency(1,1e6,10,1/50,2000) #commented tobe readable in debian

 simu_seq=Float64[]
  for k in 1:size(x1)[2]#each generation
     simu=count(a->a!=0,x1[1:end-1,k]) #established allele--simulation
     push!(simu_seq,simu)
   end#for k in generations
   plot(simu_seq)
