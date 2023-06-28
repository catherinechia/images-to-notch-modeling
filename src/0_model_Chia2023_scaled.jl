##############################################
#Author: Catherine Chia
#Contact: catherine.mmc@gmail.com
#Refer to in-text comments for respective references

#Added DAPT stimulus
#Monk N. A. (2003). Oscillatory expression of Hes1, p53, and NF-kappaB driven by transcriptional time delays. Current biology : CB, 13(16), 1409–1413. https://doi.org/10.1016/s0960-9822(03)00494-9

##############################################
#Libraries, packages
#None

# Define delayed differential equations
function chia2023_scaled_model(du, u, h, p, t)
    #Parameters
    mu_m, mu_p, tau, n, rep_rescaled = p
    
    #Delay
    p_rescaled_hist = h(p, t - tau)[2]

    #Drug-specific constants
    V_max = 0.8 #unitless   
    K_m = 1.0   #concentration

    #DDES
    D = 1 - (V_max * DAPT_const[1]/ (K_m + DAPT_const[1]))
    g = 1/(1 + (p_rescaled_hist/rep_rescaled)^n)
    du[1] = D*g - (mu_m * u[1]) 
    du[2] = u[1] - (mu_p * u[2])
end


#History function 
h(p, t) = ones(2)

#Initial set of parameters and values
#Scaled
d_params_scaled = Dict(
    :mu_m => 0.03,           #1 (1/min),  degradation rate
    :mu_p => 0.03,           #2 (1/min),  degradation rate  
    :tau => 18.5,            #3 (min),    delay
    :n => 5.0,               #4 (unitless) Hill coefficient
    :rep => 100.0,           #5 (unitless)
)  

p = [d_params_scaled[:mu_m], d_params_scaled[:mu_p], d_params_scaled[:tau], d_params_scaled[:n], d_params_scaled[:rep]]    
    
l_lags = [d_params_scaled[:tau]]
                

##Initial conditions
d_initial = Dict(
    :M0 => 3.0,                # expression level
    :P0 => 100.0               # expression level
)

l_u0 = [d_initial[:M0], d_initial[:P0]]



