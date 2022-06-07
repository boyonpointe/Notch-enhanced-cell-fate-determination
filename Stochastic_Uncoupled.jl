using DifferentialEquations
using Random
using NPZ
#================== simple morphogen  gradient ==================#
function stoc_uncoupled_F!(du,u,p,t)
    N,D_M,alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiA,xiB,tON = p
    M = u[1:N]
    A = u[(N+1):(2*N)]
    B = u[(2*N+1):(3*N)]
    for i in 1:N
        if i == 1
            du[i] = alphaM + D_M*(M[2] - M[1]) - tauM*M[i]
        elseif i == N
            du[i] = D_M*(M[N-1] - M[N]) - tauM*M[i]
        else
            du[i] = D_M*(M[i-1] + M[i+1] -2*M[i]) - tauM*M[i]
        end
        #----------------------------------------------------#
        if (t > tON)
            du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) - A[i]/tauA
            du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) - B[i]/tauB
        end
        #----------------------------------------------------#
    end
end

#========================================================================#
function stoc_uncoupled_G!(du,u,p,t)
    N,D_M,alphaM,tauM,K1,K2,K3,K4,alphA,alphaB,tauA,tauB,h,xiM,xiA,xiB,tON = p
    du[1:N]           = xiM * u[1:N]
    du[(N+1):(2*N)]   = xiA * u[(N+1):(2*N)]
    du[(2*N+1):(3*N)] = xiB * u[(2*N+1):(3*N)]
end
#==========================================================================#
function trial_run()
    N = 50 ;
    D_M = 50.0 ;
    alphaM = 500*N ;
    tauM   = 1.0 ;
    K1,K2,K3,K4 = 50,200,100,200 ;
    h = 2 ;

    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;
    xiM,xiA,xiB = 1.0,0.1,0.1 ;

    tmax = 60.0 ;
    tON  = 30.0 ;
    p = N,D_M,alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiA,xiB,tON ;

    u0 = zeros(3*N) ;
    tspan = (0.0,tmax) ;
    prob  = SDEProblem(stoc_uncoupled_F!,stoc_uncoupled_G!,u0,tspan,p) ;
    sol   = solve(prob,SRIW1()) ;
    println(size(sol))
    npzwrite("/home/chandrashekar/Documents/TooBigForDropbox/Boundary/tmp_uncoupled_t.npy",convert(Array,sol.t)) ;
    npzwrite("/home/chandrashekar/Documents/TooBigForDropbox/Boundary/tmp_uncoupled.npy",convert(Array,sol)) ;
end
#=======================================================================#

function sample_ensemble(iters)
    N = 50 ;
    D_M = 50.0 ;
    alphaM = 500*N ;
    tauM   = 1.0 ;
    K1,K2,K3,K4 = 50,200,100,200 ;

    alphaA,alphaB,tauA,tauB,h = 2000,2000,0.25,0.25,2 ;
    xiM,xiA,xiB = 1.0,0.1,0.1 ;
    tmax = 60.0 ;
    tON  = 30.0 ;
    p = N,D_M,alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiA,xiB,tON ;
    GeneExp = zeros(2,iters,N) ;
    ExpData = zeros(6,iters,N) ;
    for i = 1:iters
        println(i)
        u0 = zeros(3*N) ;
        tspan = (0.0,tmax) ;
        prob  = SDEProblem(stoc_uncoupled_F!,stoc_uncoupled_G!,u0,tspan,p) ;

        sol  = solve(prob,SRIW1(),force_dtmin=true) ;

        idx = sum(sol.t .< 31) ;
        A = sol[51:100,idx:end] ;
        B = sol[101:150,idx:end] ;

        ExpData[1,i,:] = A[:,end] ;
        ExpData[2,i,:] = B[:,end] ;
        ExpData[3,i,:] = mean(A,dims=2) ;
        ExpData[4,i,:] = mean(B,dims=2) ;
        ExpData[5,i,:] = var(A,dims=2) ;
        ExpData[6,i,:] = var(B,dims=2) ;

        GeneExp[1,i,:] = sol[N+1:2*N,end] ;
        GeneExp[2,i,:] = sol[2*N+1:3*N,end] ;
    end
    #npzwrite("/home/chandrashekar/Documents/TooBigForDropbox/Boundary/Uncoupled_ensemble_10_1_1_niters_"*string(iters)*".npy",GeneExp) ;
    return ExpData ;
end
#========================================================================#
#trial_run()

#sample_ensemble(10000);

ED = sample_ensemble(500) ;
#npzwrite("Data/tauS/timeseries_spec_uncoupled_niters_"*string(iters)*".npy",ED) ;
