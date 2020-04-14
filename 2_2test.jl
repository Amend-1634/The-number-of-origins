cd("D:\\Julia\\func")#pathway to your working directory
#fake frequecy vector generated
function simu_p(n)
    x=randn(n)
    x=x/sum(x)
    return x
 end
 xx=simu_p(1000)
 N=1e6
include("2_2gaussian_multinomial2.jl") #former version
 #transformed into positive semidefinite matrix
  #but this no longer work any more
   #as MvNormal() require positive definite and Hermitian matrix
@time z=gaussian_multinomial(N,xx)


include("2_2gaussian_multinomial.jl")
 #amended--transformed into positive definite
  #but line82 cannot filter the D(eigenvalue matrix composed by complex number)
   #by integer
@time z=gaussian_multinomial(N,xx)

sum(z)
println(xx)
plot(z) #same output as Matlab code
