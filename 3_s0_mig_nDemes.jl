#s0_vec,series: mig_vec nDemes; 2subplot: N
function onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0_vec,mig_vec,M,T,Gauss)
    len_mig=length(mig_vec)
    len_nDemes=length(nDemes_vec)
    len_s0=length(s0_vec)

    #y: time to fixation
    t_avn_mtx=zeros(len_s0,len_nDemes*len_mig)
    t_se_mtx=zeros(len_s0,len_nDemes*len_mig)
    eta_avn_mtx=zeros(len_s0,len_nDemes*len_mig)
    eta_se_mtx=zeros(len_s0,len_nDemes*len_mig)
    fix_M_mtx=zeros(len_s0,len_nDemes*len_mig)

    for i in 1:len_nDemes
        for j in 1:len_mig
            if len_nDemes==1
                for k in 1:len_s0
                    println("nDemes=$i, mig=$j, s0=$k")
                    t_avn, t_se, eta_avn, eta_se, fix_M=
                    num_origin_1d_fix(N,nDemes_vec[i],twoNmu,s0_vec[k],mig_vec[j],M,T,Gauss)
                    t_avn_mtx[k,j]=t_avn
                    t_se_mtx[k,j]=t_se
                    eta_avn_mtx[k,j]=eta_avn
                    eta_se_mtx[k,j]=eta_se
                    fix_M_mtx[k,j]=fix_M
                end
            else
                for k in 1:len_s0
                    t_avn, t_se, eta_avn, eta_se, fix_M=
                    num_origin_1d_fix(N,nDemes_vec[i],twoNmu,s0_vec[k],mig_vec[j],M,T,Gauss)
                    t_avn_mtx[k,j+2]=t_avn
                    t_se_mtx[k,j+2]=t_se
                    eta_avn_mtx[k,j+2]=eta_avn
                    eta_se_mtx[k,j+2]=eta_se
                    fix_M_mtx[k,j+2]=fix_M
                end
            end#if
        end#for j
    end#for i

    legend=Array{String,2}(undef,1,len_nDemes*len_mig)
        for i in 1:len_nDemes
            for j in 1:len_mig
                legend[k,j] = "nDemes = $(nDemes_vec[i]), s0_vec = $(s0_vec[j])"
            end #for j
        end #for i
    fix_M_mtx=fix_M_mtx/M #fixation probability
    return legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx
end
