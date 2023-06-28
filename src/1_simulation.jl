##############################################
#Author: Catherine Chia
#Contact: catherine.mmc@gmail.com

##Refer to in-text comments for respective references
##############################################

# Jump to the righty folder
# if !occursin("sysbio", pwd())
#     cd("src/sysbio")
# ends

# Load model
include("0_model_Chia2023_scaled.jl");


#Libraries, packages
using DifferentialEquations
using Plots


#Callback functions 
#https://discourse.julialang.org/t/solving-ode-parameters-using-experimental-data-with-control-inputs/66614/20
#handling discrete callback tstops
#https://docs.sciml.ai/DiffEqDocs/latest/features/callback_functions/#SciMLBase.DiscreteCallback
#https://github.com/SciML/DiffEqCallbacks.jl/issues/44

time_start  = 0.0   #min
time_end    = 1024.0  #min
DAPT_start  = [140.0]   #min
DAPT_stop   = [700.0]   #min
tstops_DAPT_ON = [DAPT_start; DAPT_stop]
tspan = (time_start[1], time_end[1])
tstep=7 #1 step = 7 mins

##########################
# Define DDE problem
prob = DDEProblem(chia2023_scaled_model, l_u0, h, tspan, p; constant_lags = l_lags)

##########################
#Initial value of DAPT 
const DAPT_const=[0.0]


#ON OFF conditions
function condition_DAPTON(u, t, integrator)
    t in DAPT_start
end
function condition_DAPTOFF(u, t, integrator)
    t in DAPT_stop
end

#Use @show for debugging https://github.com/SciML/DifferentialEquations.jl/issues/805
function affect_control!(integrator)
    DAPT_const[1] = 0.0 # uM
    @show integrator.t, DAPT_const[1]
end

function affect_treat1!(integrator)
    DAPT_const[1] = 1.0 # uM
    @show integrator.t, DAPT_const[1]
end

function affect_treat2!(integrator)
    DAPT_const[1] = 5.0 # uM
    @show integrator.t, DAPT_const[1]
end

function affect_treat3!(integrator)
    DAPT_const[1] = 10.0
    @show integrator.t, DAPT_const[1]
end


#https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/
#save_positions=(true,true): Boolean tuple for whether to save before and after the affect!. 
cb_DAPTON_control = DiscreteCallback(condition_DAPTON, affect_control!, save_positions=(false, true))
cb_DAPTON_T1000 = DiscreteCallback(condition_DAPTON, affect_treat1!, save_positions=(false, true))
cb_DAPTON_T5000 = DiscreteCallback(condition_DAPTON, affect_treat2!, save_positions=(false, true))
cb_DAPTON_T10000 = DiscreteCallback(condition_DAPTON, affect_treat3!, save_positions=(false, true))
cb_DAPTOFF = DiscreteCallback(condition_DAPTOFF, affect_control!, save_positions=(false, true))

cbs_control = CallbackSet(cb_DAPTON_control, cb_DAPTOFF)
cbs_T1000 = CallbackSet(cb_DAPTON_T1000, cb_DAPTOFF)
cbs_T5000 = CallbackSet(cb_DAPTON_T5000, cb_DAPTOFF)
cbs_T10000 = CallbackSet(cb_DAPTON_T10000, cb_DAPTOFF)


# Solving DDE problem with initial parameters
sol_init_control = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_control, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T1000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T1000, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T5000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T5000, tstops=tstops_DAPT_ON, saveat=tstep)
sol_init_T10000 = solve(prob, MethodOfSteps(Tsit5()), callback=cbs_T10000, tstops=tstops_DAPT_ON, saveat=tstep)


