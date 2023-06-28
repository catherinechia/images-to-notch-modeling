##############################################
#Author: Catherine Chia
#Contact: catherine.mmc@gmail.com

##Main references
#https://docs.sciml.ai/GlobalSensitivity/stable/tutorials/juliacon21/

##Refer to in-text comments for respective references
##############################################
# Launch REPL with 
# $ julia --threads 4
# Threads.nthreads()
] generate sens_env
] activate sens_env

# Jump to the righty folder
# if !occursin("sysbio", pwd())
#     cd("src/sysbio")
# ends

##############################################
# Load model
include("0_model_Chia2023_scaled.jl");

#Pycall
using PyCall
using Conda
import Pandas

#For DDE
using DifferentialEquations
using RecursiveArrayTools # for VectorOfArray
using Plots
using DataFrames
using DataFramesMeta

#Global sensitivity analysis: Periods, Amplitude, Phase
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using CairoMakie
using DifferentialEquations
using RecursiveArrayTools # for VectorOfArray

#Call Python package - pyboat
# Installing pyBOAT in its own environment - Once
# Conda.add("pyboat")
# Conda.add("scipy")

# Include core script
@pyinclude("../anaconda3/envs/i2n/lib/python3.10/site-packages/pyboat/api.py")
@pyinclude("../anaconda3/envs/i2n/lib/python3.10/site-packages/pyboat/core.py")

##############################################
#Objective functions

fn_pyboat = function (l_signal_sol1, l_signal_sol2, measurement_type)
    # --- set basic parameters and initialize the Analyzer---
    duration_per_timepoint = (120 + 800 + 120) / 147 # 7.074829931972789 min
    dt = round(duration_per_timepoint, digits=3) # the sampling interval, 7.075 min
    lowest_period = 15 #30.0    #min. half-life and fluorescent maturation time (15mins)
    highest_period = 400.0 #240.0  #min. 4hrs x 60min/hr from Ultradian range 2-4 hours
    period_stepsize = 1.0   #1.05 used in dataprocessingV7   
    periods = collect(lowest_period:period_stepsize:highest_period)# period range
    wAn_sol1 = py"WAnalyzer"(periods, dt, time_unit_label="min")
    wAn_sol2 = py"WAnalyzer"(periods, dt, time_unit_label="min")
    ridge_threshold = 0.0 #min, assume no background noise

    #Wavelet transforms - Compute spectrum
    ar_spectrum_sol1 = wAn_sol1.compute_spectrum(l_signal_sol1, do_plot="False") #pyboat API uses self dt, periods
    ar_spectrum_sol2 = wAn_sol2.compute_spectrum(l_signal_sol2, do_plot="False") 

    #Ridge evaluation and readout
    pydf_norm_ridge_sol1 = wAn_sol1.get_maxRidge(power_thresh=ridge_threshold)
    pydf_norm_ridge_sol2 = wAn_sol2.get_maxRidge(power_thresh=ridge_threshold)

    #Convert python dataframe to julia dataframe
    df_norm_ridge_sol1 = pydf_norm_ridge_sol1 |> Pandas.DataFrame |> DataFrames.DataFrame
    df_norm_ridge_sol2 = pydf_norm_ridge_sol2 |> Pandas.DataFrame |> DataFrames.DataFrame
    
    if measurement_type=="periods"
        l_periods_sol1 = df_norm_ridge_sol1.periods
        l_periods_sol2 = df_norm_ridge_sol2.periods
        res_sol1 = median(l_periods_sol1)
        res_sol2 = median(l_periods_sol2)
    elseif measurement_type=="amplitude"
        l_amplitude_sol1 = df_norm_ridge_sol1.amplitude
        l_amplitude_sol2 = df_norm_ridge_sol2.amplitude
        res_sol1 = median(l_amplitude_sol1)
        res_sol2 = median(l_amplitude_sol2)
    elseif measurement_type=="phase"
        l_phase_sol1 = df_norm_ridge_sol1.phase 
        l_phase_sol2 = df_norm_ridge_sol2.phase 
        res_sol1 = median(l_phase_sol1)
        res_sol2 = median(l_phase_sol2)
    else

    end
    f1_res = [res_sol1, res_sol2] 
    
    return f1_res
end

##Specific objective functions
f1_periods = function (p)
    prob1 = remake(prob;p=p)
    sol = solve(prob1, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)
    f1_res = fn_pyboat(sol[1,:], sol[2,:], "periods")
    return f1_res
end

f1_amplitude = function (p)
    prob1 = remake(prob;p=p)
    sol = solve(prob1, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)
    f1_res = fn_pyboat(sol[1,:], sol[2,:], "amplitude")
    return f1_res
end

f1_phase = function (p)
    prob1 = remake(prob;p=p)
    sol = solve(prob1, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)
    f1_res = fn_pyboat(sol[1,:], sol[2,:], "phase")
    return f1_res
end

##############################################
#Extract sol period and amplitude
sol_init_control = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_control, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T1000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T1000, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T5000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T10000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T10000, tstops=tstops_DAPT_ON, saveat=tstep)

#Set basic parameters and initialize the Analyzer
duration_per_timepoint = 7.074829931972789 #min 
dt = 7.0 # the sampling interval, min
lowest_period = 15 #30.0    #min. half-life and fluorescent maturation time (15mins)
highest_period = 400.0 #240.0  #min. 4hrs x 60min/hr from Ultradian range 2-4 hours
period_stepsize = 1.0     
periods = collect(lowest_period:period_stepsize:highest_period)# period range
wAn_sol2_control = py"WAnalyzer"(periods, dt, time_unit_label="min")
wAn_sol2_T1000 = py"WAnalyzer"(periods, dt, time_unit_label="min")
wAn_sol2_T5000 = py"WAnalyzer"(periods, dt, time_unit_label="min")
wAn_sol2_T10000 = py"WAnalyzer"(periods, dt, time_unit_label="min")
ridge_threshold = 0.0 #min, assume no background noise

#Normalize to first time point
sol_fc_control = stack(sol_init_control.u) ./ sol_init_control.u[1]
sol_fc_T1000 = stack(sol_init_T1000.u) ./ sol_init_T1000.u[1]
sol_fc_T5000 = stack(sol_init_T5000.u) ./ sol_init_T5000.u[1]
sol_fc_T10000 = stack(sol_init_T10000.u) ./ sol_init_T10000.u[1]
l_signal_sol2_control = sol_fc_control[2,:]
l_signal_sol2_T1000 = sol_fc_T1000[2,:]
l_signal_sol2_T5000 = sol_fc_T5000[2,:]
l_signal_sol2_T10000 = sol_fc_T10000[2,:]



###########################################
##Sobol Global Sensitivity Analysis
# Define DDE problem
l_u0 = [3.76, 143.0]
prob = DDEProblem(chia2023_scaled_model, l_u0, h, tspan, p; constant_lags = l_lags)

#Test solve
sol_init_control = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_control, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T5000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)


#Global sensitivity analysis
bounds = [[0.001,0.07],[0.001,0.07],[1.0,31.0],[2.0,11.0],[50.0,300.0],[0.1,7.0],[0.1, 7.0]]

sobol_sens = gsa(f1_periods, Sobol(), bounds, samples = 200)
sobol_sens = gsa(f1_amplitude, Sobol(), bounds, samples = 200)
sobol_sens_phase = gsa(f1_phase, Sobol(), bounds, samples = 200)

fig = Figure(resolution = (600, 400))
CairoMakie.barplot(fig[1,1], [1,2,3,4,5], sobol_sens_phase.S1[1, :], color = :darkblue, axis = (xticksvisible = false, xticklabelsvisible = false, title = "Hes1 mRNA", ylabel = "First order"))
CairoMakie.barplot(fig[2,1], [1,2,3,4,5], sobol_sens_phase.ST[1, :], color = :darkblue, axis = (xticksvisible = false, xticklabelsvisible = false, ylabel = "Total order"))
CairoMakie.barplot(fig[1,2], [1,2,3,4,5], sobol_sens_phase.S1[2, :], color = :green, axis = (xticksvisible = false, xticklabelsvisible = false, title = "Hes1 Protein", ylabel = "First order"))
CairoMakie.barplot(fig[2,2], [1,2,3,4,5], sobol_sens_phase.ST[2, :], color = :green, axis = (xticksvisible = false, xticklabelsvisible = false, ylabel = "Total order"))
fig





