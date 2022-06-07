using DifferentialEquations
using Random
using NPZ
using Plots
using Statistics
#=================================================================#
function stoc_coupled_F!(du,u,p,t)
    N,D_M,alphaM,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,ctype,xiM,xiA,xiB,xiR,xiL,xiS = p

    M = u[1:N]
    A = u[(N+1):(2*N)]
    B = u[(2*N+1):(3*N)]
    R = u[(3*N+1):(4*N)]
    L = u[(4*N+1):(5*N)]
    S = u[(5*N+1):(6*N)]

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
            if ctype == 0
                du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) - A[i]/tauA
                du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) - B[i]/tauB
            elseif ctype == 1
                du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) * ((Q^g)/(Q^g + S[i]^g)) - A[i]/tauA
                du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) - B[i]/tauB
            elseif ctype == 2
                du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) - A[i]/tauA
                du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) * ((Q^g)/(Q^g + S[i]^g)) - B[i]/tauB
            elseif ctype == 3
                du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) +  delta*((S[i]^g)/(Q^g + S[i]^g)) - A[i]/tauA
                du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) - B[i]/tauB
            elseif ctype == 4
                du[1*N + i] = alphaA * ((M[i]^h)/(K1^h + M[i]^h)) * ((K3^h)/(K3^h + B[i]^h)) - A[i]/tauA
                du[2*N + i] = alphaB * ((M[i]^h)/(K2^h + M[i]^h)) * ((K4^h)/(K4^h + A[i]^h)) +  delta*((S[i]^g)/(Q^g + S[i]^g)) - B[i]/tauB
            end
            #----------------------------------------------------------------#
            if i == 1
                L_trans = L[2]
                R_trans = R[2]
            elseif i == N
                L_trans = L[N-1]
                R_trans = R[N-1]
            else
                L_trans = L[i-1] + L[i+1]
                R_trans = R[i-1] + R[i+1]
            end
            #----------------------------------------------------------------#
            du[3*N + i] = betaR0 + betaR * ((A[i]^g)/(J^g + A[i]^g)) * ((B[i]^g)/(J^g + B[i]^g)) - R[i]/tauR - k_trans*R[i]*L_trans
            du[4*N + i] = betaL0 + betaL * ((W^g)/(W^g + S[i]^g)) - L[i]/tauL - k_trans*L[i]*R_trans
            du[5*N + i] = k_trans*R[i]*L_trans - S[i]/tauS
        end
    end
end
#============================================================================#
function stoc_coupled_G!(du,u,p,t)
    N,D_M,alphaM,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,ctype,xiM,xiA,xiB,xiR,xiL,xiS = p
    du[1:N]           = xiM * u[1:N]
    du[(N+1):(2*N)]   = xiA * u[(N+1):(2*N)]
    du[(2*N+1):(3*N)] = xiB * u[(2*N+1):(3*N)]
    du[(3*N+1):(4*N)] = xiR * u[(3*N+1):(4*N)]
    du[(4*N+1):(5*N)] = xiL * u[(4*N+1):(5*N)]
    du[(5*N+1):(6*N)] = xiS * u[(5*N+1):(6*N)]
end
#===============================================================================#
function trial_run()
    N = 50
    D_M,alphaM,tauM = 50.0,500*N,1.0
    K1,K2,K3,K4 = 50,200,100,200

    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;


    betaR0,betaL0,betaR,betaL,k_trans = 20.0,20.0,200.0,200.0,1.0
    tauR,tauL,tauS = 0.2,0.2,2.0

    ktau_const = 2.0 ;
    #S_half_max = (k_trans*0.5*betaR*tauR * 0.5*0.5*betaL*tauL * tauS) ;
    S_half_max = (ktau_const*0.5*betaR*tauR * 0.5*0.5*betaL*tauL) ;
    println("S_half_max = "*string(S_half_max)) ;


    W = 5*S_half_max #explore different Ws
    h,g = 2,4
    tON = 30.0

    tmax = 60.0

    #----------------------------#
    ctype = 2
    Q,J = 0.5*S_half_max,44.67 #bad : Q = 1.5
    delta = 0.0 * alphaA
    #----------------------------#
    #ctype = 3
    #Q,J = 0.5 * S_half_max, 96.78
    #delta = 0.0 * alphaA
    #----------------------------#

    xiM,xiA,xiB,xiR,xiL,xiS = 1.0,0.1,0.1,0.01,0.01,0.01 ;


    p = N,D_M,alphaM,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,ctype,xiM,xiA,xiB,xiR,xiL,xiS ;

    u0 = zeros(6*N) ; # it cant be zeros -  the R,L should start at max
    u0[(3*N+1):(4*N)] .= betaR0 * tauR ;
    u0[(4*N+1):(5*N)] .= (betaL0 + betaL) * tauL ;

    tspan = (0.0,tmax) ;
    prob  = SDEProblem(stoc_coupled_F!,stoc_coupled_G!,u0,tspan,p) ;
    sol   = solve(prob,SRIW1(),force_dtmin=true) ;
    # println(size(sol))
    # npzwrite("/home/chandrashekar/Documents/TooBigForDropbox/Boundary/tmp_coupled_shortlivedS_t.npy",convert(Array,sol.t)) ;
    # npzwrite("/home/chandrashekar/Documents/TooBigForDropbox/Boundary/tmp_coupled_shortlivedS.npy",convert(Array,sol)) ;
    return sol ;
end
#================================================================================#
#trial_run()
function sample_ensemble(tauS,iters,ctype,Q,J,delta)
    N = 50
    D_M,alphaM,tauM = 50.0,500*N,1.0
    K1,K2,K3,K4 = 50,200,100,200

    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;


    #betaR0,betaL0,betaR,betaL,k_trans = 20.0,20.0,200.0,200.0,1.0

    betaR0,betaL0,betaR,betaL = 20.0,20.0,200.0,200.0 ;

    tauR,tauL = 0.2,0.2

    ktau_const = 2.0 ;

    k_trans = ktau_const/tauS ;


    S_half_max = (ktau_const*0.5*betaR*tauR * 0.5*0.5*betaL*tauL) ;

    #S_half_max = (k_trans*0.5*betaR*tauR * 0.5*0.5*betaL*tauL * tauS)

    W = 5*S_half_max #explore different Ws = 5*S_half_max is quite large(=2000)
    h,g = 2,4
    tON = 30.0

    tmax = 60.0

    Q = Q * S_half_max ;
    println(Q) ;

    delta = delta * alphaA ;


    xiM,xiA,xiB,xiR,xiL,xiS = 1.0,0.1,0.1,0.01,0.01,0.01 ;
    #xiM,xiA,xiB,xiR,xiL,xiS = 1.0,0.1,0.1,0.1,0.1,0.1 ;

    p = N,D_M,alphaM,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,ctype,xiM,xiA,xiB,xiR,xiL,xiS ;

    GeneExp = zeros(2,iters,N)

    for i = 1:iters
        println(i)
        u0 = zeros(6*N) ; # it cant be zeros -  the R,L should start at max
        u0[(3*N+1):(4*N)] .= betaR0 * tauR ;
        u0[(4*N+1):(5*N)] .= (betaL0 + betaL) * tauL ;
        tspan = (0.0,tmax) ;
        prob  = SDEProblem(stoc_coupled_F!,stoc_coupled_G!,u0,tspan,p) ;
        sol  = solve(prob,SRIW1(),force_dtmin=true) ;
        GeneExp[1,i,:] = sol[N+1:2*N,end] ;
        GeneExp[2,i,:] = sol[2*N+1:3*N,end] ;
    end
    #npzwrite("/home/kuyyamudi/Documents/TooBigForDropbox/Data/Coupled_Ensemble_niters_"*string(iters)*"tauS_"*string(tauS)*".npy",GeneExp);
    npzwrite("Data/tauS/Coupled_Ensemble_niters_"*string(iters)*"tauS_"*string(tauS)*".npy",GeneExp);
    return GeneExp
end
#===========================================================================================#
function sample_ensemble_timeseries(tauS,iters,ctype,Q,J,delta)
    N = 50
    D_M,alphaM,tauM = 50.0,500*N,1.0
    K1,K2,K3,K4 = 50,200,100,200

    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;


    #betaR0,betaL0,betaR,betaL,k_trans = 20.0,20.0,200.0,200.0,1.0

    betaR0,betaL0,betaR,betaL = 20.0,20.0,200.0,200.0 ;

    tauR,tauL = 0.2,0.2

    ktau_const = 2.0 ;

    k_trans = ktau_const/tauS ;


    S_half_max = (ktau_const*0.5*betaR*tauR * 0.5*0.5*betaL*tauL) ;

    #S_half_max = (k_trans*0.5*betaR*tauR * 0.5*0.5*betaL*tauL * tauS)

    W = 5*S_half_max #explore different Ws = 5*S_half_max is quite large(=2000)
    h,g = 2,4
    tON = 30.0

    tmax = 60.0

    Q = Q * S_half_max ;
    println(Q) ;

    delta = delta * alphaA ;


    xiM,xiA,xiB,xiR,xiL,xiS = 1.0,0.1,0.1,0.01,0.01,0.01 ;
    #xiM,xiA,xiB,xiR,xiL,xiS = 1.0,0.1,0.1,0.1,0.1,0.1 ;

    p = N,D_M,alphaM,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,ctype,xiM,xiA,xiB,xiR,xiL,xiS ;

    ExpData = zeros(6,iters,N) #A at tmax, B at tmax, meanA, meanB, varA, varB

    for i = 1:iters
        u0 = zeros(6*N) ; # it cant be zeros -  the R,L should start at max
        u0[(3*N+1):(4*N)] .= betaR0 * tauR ;
        u0[(4*N+1):(5*N)] .= (betaL0 + betaL) * tauL ;
        tspan = (0.0,tmax) ;
        prob  = SDEProblem(stoc_coupled_F!,stoc_coupled_G!,u0,tspan,p) ;
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
    end
    #npzwrite("/home/kuyyamudi/Documents/TooBigForDropbox/Data/Coupled_Ensemble_niters_"*string(iters)*"tauS_"*string(tauS)*".npy",GeneExp);
    #npzwrite("Data/tauS/Coupled_Ensemble_ts_spec_niters_"*string(iters)*"tauS_"*string(tauS)*".npy",ExpData);
    return ExpData ;
end


#===========================================================================================#
function Explore_hill(ctype,delta,iters)
    N = 50 ;
    Q = LinRange(0.5,2.0,20) ;
    J = LinRange(1,100,15) ;
    OD = zeros(20,15,2,iters,N) ;
    for i = 1:20
        for j = 1:15
            println(i,":",j)
            G = sample_ensemble(iters,ctype,Q[i],J[j],delta) ;
            OD[i,j,:,:,:] = G ;
        end
    end
    npzwrite("Data/Stochastic_Q_J_ctype_2.npy",OD) ;
end
#=========================================================================================#
function Explore_J(ctype,delta,iters,qid)
    N = 50 ;
    QVals = LinRange(0.5,2.0,40) ;
    Q = QVals[qid] ;
    J = LinRange(1,100,30) ;
    OD = zeros(30,2,iters,N) ;
    for i = 1:30
        println(i)
        G = sample_ensemble(iters,ctype,Q,J[i],delta) ;
        OD[i,:,:,:] = G ;
    end
    ofilename = "Data/Exp_J_ctype_"*string(ctype)*"_qid_"*string(qid)*".npy" ;
    npzwrite(ofilename,OD)
end

#==========================================================================================#


iters = 500 ;
N = 50 ;
#tauS_Vals = LinRange(0.75,1.25,30) ;
tauS_Vals = LinRange(0.5,2.5,100) ;
OD = zeros(length(tauS_Vals),6,iters,N) ;

for i = 1:length(tauS_Vals)
    tauS = tauS_Vals[i] ;
    println(i) ;
    println(tauS) ;

    ctype,Q,J,delta = 2,0.5,44.67,0.5 ; #tauAB = 0.25
    specs = sample_ensemble_timeseries(tauS,iters,ctype,Q,J,delta) ;
    OD[i,:,:,:] = specs ;
end


