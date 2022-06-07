using DifferentialEquations
using Random
using NPZ
using Plots
using Statistics
#=================================================================#

function stoc_uncoupled_F!(du,u,p,t)
    alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiAB,tON = p ;
    M1,A1,B1 = u ;
    du[1] = alphaM - M1/tauM ;
    if t > tON
        du[2] = alphaA * ((M1^h)/(K1^h + M1^h)) * ((K3^h)/(K3^h + B1^h)) - A1/tauA ;
        du[3] = alphaB * ((M1^h)/(K2^h + M1^h)) * ((K4^h)/(K4^h + A1^h)) - B1/tauB ;
    else
        du[2] = 0.0 ;
        du[3] = 0.0 ;
    end
end

function stoc_uncoupled_G!(du,u,p,t)
    alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiAB,tON = p ;
    du[1] = xiM * u[1] ;
    du[2] = xiAB * u[2] ;
    du[3] = xiAB * u[3] ;
end

#=================================================================#
#=
Should consider noise in all terms and not just M
=#
function stoc_coupled_two_cells_F!(du,u,p,t)
    alphaM1,alphaM2,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,xiM,xiAB,xiS = p ;
    M1,A1,B1,R1,L1,S1,M2,A2,B2,R2,L2,S2 = u ;

    du[1] = alphaM1 - M1/tauM ;
    du[2] = alphaA * ((M1^h)/(K1^h + M1^h)) * ((K3^h)/(K3^h + B1^h)) - A1/tauA ;
    du[3] = alphaB * ((M1^h)/(K2^h + M1^h)) * ((K4^h)/(K4^h + A1^h)) * ((Q^g)/(Q^g + S1^g)) - B1/tauB ;
    du[4] = betaR0 + betaR * ((A1^g)/(J^g + A1^g)) * ((B1^g)/(J^g + B1^g)) - R1/tauR - k_trans*R1*L2 ;
    du[5] = betaL0 + betaL * ((W^g)/(W^g + S1^g)) - L1/tauL - k_trans*L1*R2
    du[6] = k_trans*R1*L2 - S1/tauS ;

    du[7]  = alphaM2 - M2/tauM ;
    du[8]  = alphaA * ((M2^h)/(K1^h + M2^h)) * ((K3^h)/(K3^h + B2^h)) - A2/tauA ;
    du[9]  = alphaB * ((M2^h)/(K2^h + M2^h)) * ((K4^h)/(K4^h + A2^h)) * ((Q^g)/(Q^g + S2^g)) - B2/tauB ;
    du[10] = betaR0 + betaR * ((A2^g)/(J^g + A2^g)) * ((B2^g)/(J^g + B2^g)) - R2/tauR - k_trans*R2*L1 ;
    du[11] = betaL0 + betaL * ((W^g)/(W^g + S2^g)) - L2/tauL - k_trans*L2*R1
    du[12] = k_trans*R2*L1 - S2/tauS ;
end

#===================================================================================#
function stoc_coupled_two_cells_G!(du,u,p,t)
    alphaM1,alphaM2,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,xiM,xiAB,xiS = p ;
    du[1] = xiM * u[1] ;
    du[7] = xiM * u[7] ;

    du[2] = xiAB * u[2] ;
    du[3] = xiAB * u[3] ;
    du[8] = xiAB * u[8] ;
    du[9] = xiAB * u[9] ;

    du[4] = xiS * u[4] ;
    du[5] = xiS * u[5] ;
    du[6] = xiS * u[6] ;
    du[10] = xiS * u[10] ;
    du[11] = xiS * u[11] ;
    du[12] = xiS * u[12] ;
end
#===================================================================================#
#ktau_const = 2.0 ;

function trial_run_uncoupled(alphaM)
    xiM = 0.1 ;
    xiAB = 0.1 ;
    tauM = 1.0 ;
    K1,K2,K3,K4 = 50,200,100,200 ;
    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;
    h = 2 ;
    tON = 30.0 ;
    tmax = 60.0 ;
    #println(alphaM) ;
    p = alphaM,tauM,K1,K2,K3,K4,alphaA,alphaB,tauA,tauB,h,xiM,xiAB,tON ;
    u0 = zeros(3) ;
    u0[1] = alphaM*tauM ;
    tspan = (0.0,tmax) ;
    prob  = SDEProblem(stoc_uncoupled_F!,stoc_uncoupled_G!,u0,tspan,p) ;
    sol   = solve(prob,SRIW1(),force_dtmin=true) ;

    idx = sum(sol.t .< 30) ;
    sz = size(sol)[2]-idx + 1 ;
    # #-----------------------------------------------#
    Output = zeros(4,sz) ;
    Output[1,:] = sol.t[idx:end] ; #time
    Output[2,:] = sol[2,idx:end] ; #A1
    Output[3,:] = sol[3,idx:end] ; #B1
    Output[4,:] = sol[1,idx:end] ; #M1
    return Output ;
end
#==================================================================================#

function trial_run_coupled(alphaM1,alphaM2,tauS)
    xiM = 1.0 ;
    xiAB = 0.1 ;
    xiS = 0.01 ;

    tauM = 1.0 ;
    K1,K2,K3,K4 = 50,200,100,200 ;

    alphaA,alphaB,tauA,tauB = 2000,2000,0.25,0.25 ;

    betaR0,betaL0,betaR,betaL = 20.0,20.0,200.0,200.0 ;
    tauR,tauL = 0.2,0.2 ;
    ktau_const = 2.0 ;
    S_half_max = (ktau_const*0.5*betaR*tauR * 0.5*0.5*betaL*tauL) ;

    k_trans = ktau_const/tauS ;

    h,g = 2,4 ;
    tON = 30.0 ;
    tmax = 40.0 ;

    Q,J = 0.5*S_half_max,44.67
    W = 5*S_half_max #explore different Ws
    #println(Q) ;
    delta = 0.0 * alphaA

    p = alphaM1,alphaM2,tauM,K1,K2,K3,K4,Q,J,W,alphaA,alphaB,betaR0,betaL0,betaR,betaL,k_trans,tauA,tauB,tauR,tauL,tauS,h,g,delta,tON,xiM,xiAB,xiS ;
    #println(k_trans) ;

    u0 = zeros(12) ;
    u0[1] = alphaM1 * tauM ;
    u0[7] = alphaM2 * tauM ;
    tspan = (0.0,tmax) ;
    prob  = SDEProblem(stoc_coupled_two_cells_F!,stoc_coupled_two_cells_G!,u0,tspan,p) ;
    sol   = solve(prob,SRIW1(),force_dtmin=true) ;

    idx = sum(sol.t .< 10) ;
    sz = size(sol)[2]-idx + 1 ;
    # #-----------------------------------------------#
    Output = zeros(5,sz) ;
    Output[1,:] = sol.t[idx:end] ; #time
    Output[2,:] = sol[2,idx:end] ; #A1
    Output[3,:] = sol[3,idx:end] ; #B1
    Output[4,:] = sol[8,idx:end] ; #A2
    Output[5,:] = sol[9,idx:end] ; #B2
    #Output[6,:] = sol[1,idx:end] ; #M1
    #Output[7,:] = sol[7,idx:end] ; #M2
    #----------------------------------------------#
    #npzwrite("Data/TwoCell_Timeseries_tauS_"*string(tauS)*".npy",Output) ;
    return Output ;
end
#==================================================================#
function explore_tauS_and_M(iters)
    alphaM_Vals = npzread("Data/alphaM_values_50_cell.npy") ;
    tauS_Vals = LinRange(0.5,2.5,100) ;
    ncells = 30 ;

    Specs = zeros(length(tauS_Vals),ncells,iters,4) ;#means and variance
    for i = 1:length(tauS_Vals)
        println(i) ;
        tauS = tauS_Vals[i] ;
        for j = 1:ncells
            alphaM1 = alphaM_Vals[j] ;
            alphaM2 = alphaM_Vals[j+1] ;
            for q = 1:iters
                Out = trial_run_coupled(alphaM1,alphaM2,tauS) ;
                A,B = Out[2,:],Out[3,:] ;
                Specs[i,j,q,:] = [mean(A),mean(B),var(A),var(B)] ;
            end
        end
    end
    npzwrite("Data/tauS/two_cell_coupled_tauS_M_niters_"*string(iters)*".npy",Specs) ;
end
#=================================================================#
function explore_M_uncoupled(iters)
    alphaM_Vals = npzread("Data/alphaM_values_50_cell.npy") ;
    ncells = 30 ;
    Specs = zeros(ncells,iters,4) ;#means and variance
    for j = 1:ncells
        println(j) ;
        alphaM = alphaM_Vals[j] ;
        for q = 1:iters
            Out = trial_run_uncoupled(alphaM) ;
            A,B = Out[2,:],Out[3,:] ;
            Specs[j,q,:] = [mean(A),mean(B),var(A),var(B)] ;
        end
    end
    npzwrite("Data/tauS/specs_uncoupled_M_niters_"*string(iters)*".npy",Specs) ;
end
#==================================================================#
#explore_tauS_and_M(300) ;
explore_M_uncoupled(300) ;


