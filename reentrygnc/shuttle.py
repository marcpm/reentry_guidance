import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from math import sin, cos, pi, sqrt, atan, exp, tan
# from numpy import sin, cos, pi
from . import coesa
import collections
from .base import EARTH_R, G_EARTH
from .base import (lbf2N, ft2m, m2ft, m2mi, deg2rad,
                    rad2deg, lbfsqf2Nsqm, Pa2lbfsqf, 
                    Btusqft2Wsqm, Wsqm2Btusqft)
from .atmos import (get_g, get_mu, get_pressure, get_rho,
                    get_temperature, get_TempPressRhoMu)

from .flight_control import OPTIMAL_AOA, OPTIMAL_BANK, NASA_AOA, NASA_BANK

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# plotting config
plt.style.use('classic')
font = {'family' : 'monospace',
        'weight' : 'regular',
        'size'   : 7}

plt.rc('font', **font) 
plt.rc('legend',fontsize=7)

def bank_control_gells(LD_max,phi,psi):
    return 0.5 * exp(LD_max/5)*atan(cos(phi)/tan(psi))


def bank_control_chern(u, phi, psi, theta, theta_f):
    return atan(((1-u)/u) * (  (cos(phi)*sin(theta_f-theta)) / (cos(theta_f-theta)*sin(psi)
                                                - cos(psi)*sin(phi)*sin(theta_f-theta))    ))

                                                
def not_bank_chern(theta, u, phi, psi, theta_f):
        sol, _, success, msg = fsolve(bank_control_chern, 1.0*pi/180, args=[u, phi, psi, theta_f])
        if not success:
            raise RuntimeError(f"Error solving for bank angle: fsolve says: {msg}")
        return sol
    
def alpha_control(x, method="hiu"):

    if method == "hiu":
        V = x
        if V/250.0 >= 10:
            alpha = 45 
        else:
            alpha = 45 - 0.612*(V/250.0 -10)**2

    elif method == "hao":
        mach = x
        if mach > 10:
            alpha = 20
        else:
            alpha = 20 - 0.45*(mach-10)

    return (alpha) # deg

def alpha_opt_control(t):
    ts, alphas = zip(*OPTIMAL_AOA)
    interpolator = interp1d(ts, alphas)

    return interpolator(t) # returns deg for CL_CD NASA default profile

def bank_opt_control(t):
    ts, banks = zip(*OPTIMAL_BANK)
    interpolator = interp1d(ts, banks)

    return deg2rad(interpolator(t)) # returns rads

def alpha_control_NASA(V):
    Vs, alphas = zip(*NASA_AOA)
    interpolator = interp1d(Vs, alphas)

    return interpolator(V)   # returns deg for CL_CD NASA default profile


def bank_control_NASA(V):
    Vs, banks = zip(*NASA_BANK)
    interpolator = interp1d(Vs, banks)

    return deg2rad(interpolator(V)) # returns rads


def get_q_wingLE(alpha):
    
    q0 =  1.0672181
    q1 = -0.19213774e-1
    q2 =  0.21286289e-3
    q3 = -0.10117249e-5

    return q0 + q1*alpha + q2*alpha**2 + q3*alpha**3 

def get_CL_CD(alpha):
    # takes alpha in deg
    cl0 = -0.20704
    cl1 =  0.029244
    cd0 =  0.07854
    cd1 = -0.61592e-2
    cd2 =  0.621408e-3

    CL = cl0 + cl1*alpha
    CD = cd0 + cd1*alpha + cd2*alpha**2

    return CL, CD

def shuttle_odes_hicks(t, x, rho, g, spacecraft):
    # Intro Astrodynamic Reentry Book (Colonel)

    # shuttle specs todo into spacecraft class
    LD_max = spacecraft.LD_max    # LD_max = ??? 2?
    S = spacecraft.A_ref 
    m = spacecraft.m   # mass

    # ode variable aliases
    V = x[0]
    gamma = x[1]
    psi = x[2]
    r = x[3]
    theta = x[4]
    phi = x[5]

    h = r - EARTH_R
    g_h = g(h)
    temperature, pressure, rho, mu = get_TempPressRhoMu(h)
    mach = V / sqrt(1.4 * 287.1 * temperature)
    u = V**2/(g_h * h ) # energy

    ######## control laws  ##########
    
    if  "chern" in spacecraft.control_algorithm :
        theta_f = deg2rad(28.5971) # Kennedy SPace Shuttle Landing latitude
        sigma = bank_control_chern(u, phi, psi, theta, theta_f)
        AOA_control = spacecraft.control_algorithm.split("-")[-1]
        alpha = alpha_control(V, method=AOA_control)

    elif spacecraft.control_algorithm == "gear-optim":
        sigma = bank_opt_control(t)
        alpha = alpha_opt_control(t)

    elif  "gells" in spacecraft.control_algorithm:
        sigma = bank_control_gells(LD_max, phi,psi)  # bank/roll angle
        AOA_control = spacecraft.control_algorithm.split("-")[-1]
        alpha = alpha_control(V, method=AOA_control)

    elif spacecraft.control_algorithm == "nasa-guidance":
        sigma = bank_control_NASA(V)
        alpha = alpha_control_NASA(V)

    # try:

    # except:
    #     print("V: ", V)
    #     raise TypeError("Control Laws Broke")
    ###################################

    CL, CD = get_CL_CD(alpha)
    L = 0.5 * rho * CL * S * V**2
    D = 0.5 * rho * CD * S * V**2

    q = 0.5 * rho * V**2

    # intro to Astrodynamic Reentry book eqs
    Vdot = -D/m - g_h * sin(gamma)
    gammadot = 1/V * (L/m * cos(sigma) - g_h*cos(gamma) + V**2/r * cos(gamma))
    psidot = 1/V * ( (L/m * sin(sigma)/cos(gamma))  -    (V**2/r * cos(gamma)*cos(psi)*tan(phi))      )
    rdot = V * sin(gamma)
    thetadot = (V*cos(gamma)*cos(psi))/(r*cos(phi))
    phidot = V*cos(gamma)*sin(psi)/r

    return [Vdot, gammadot, psidot,rdot, thetadot, phidot]



def shuttle_odes_betts(t, x, rho, g, spacecraft):
    # Optimal Control Non Linear. Book equations

    # shuttle specs todo into spacecraft class
    LD_max = spacecraft.LD_max    # LD_max = ??? 2?
    S = spacecraft.A_ref 
    m = spacecraft.m   # mass

    # ode variable aliases
    V = x[0]
    gamma = x[1]
    psi = x[2]
    h = x[3]
    theta = x[4]
    phi = x[5]

    r = h + EARTH_R
    g_h = g(h)
    temperature, pressure, rho, mu = get_TempPressRhoMu(h)
    mach = V / sqrt(1.4 * 287.1 * temperature)
    u = V**2/(g_h * h ) # energy

    # ######## control laws  ##########
    # theta_f = deg2rad(28.5971)
    # # sigma = bank_control_gells(LD_max, phi,psi)  # bank/roll angle
    # sigma = bank_control_chern(u, phi, psi, theta, theta_f)
    # alpha = alpha_control(V, method="hiu")
    # sigma = bank_profile_control(t)
    # alpha = alpha_profile_control(t)
    # try:
    #     sigma = bank_control_NASA(V)
    #     alpha = alpha_control_NASA(V)
    # except:
    #     print("V: ", V)
    #     raise TypeError("Control Laws Broke")
    ###################################

    ######## control laws  ##########
    
    if  "chern" in spacecraft.control_algorithm :
        theta_f = deg2rad(28.5971) # Kennedy SPace Shuttle Landing latitude
        sigma = bank_control_chern(u, phi, psi, theta, theta_f)
        AOA_control = spacecraft.control_algorithm.split("-")[-1]
        alpha = alpha_control(V, method=AOA_control)
        
    elif spacecraft.control_algorithm == "gear-optim":
        sigma = bank_opt_control(t)
        alpha = alpha_opt_control(t)

    elif  "gells" in spacecraft.control_algorithm:
        sigma = bank_control_gells(LD_max, phi,psi)  # bank/roll angle
        AOA_control = spacecraft.control_algorithm.split("-")[-1]
        alpha = alpha_control(V, method=AOA_control)

    elif spacecraft.control_algorithm == "nasa-guidance":
        sigma = bank_control_NASA(V)
        alpha = alpha_control_NASA(V)

    # try:

    # except:
    #     print("V: ", V)
    #     raise TypeError("Control Laws Broke")
    ###################################




    CL, CD = get_CL_CD(alpha)
    L = 0.5 * rho * CL * S * V**2
    D = 0.5 * rho * CD * S * V**2

    q = 0.5 * rho * V**2

    # intro to Astrodynamic Reentry book eqs
    Vdot = -D/m - g_h * sin(gamma)
    gammadot = L/(m*V) * cos(sigma) + cos(gamma) * (V/r - g_h/V )
    psidot = (L* sin(sigma))/(m*V*cos(gamma)) + (V*cos(gamma)*sin(psi)*sin(theta))/(r*cos(theta))
    hdot = V * sin(gamma)
    thetadot = (V*cos(gamma)*cos(psi))/(r)
    phidot = V*cos(gamma)*sin(psi)/(cos(theta)*r)

    return [Vdot, gammadot, psidot,hdot, thetadot, phidot]


    
def run_shuttle_simulation( 
                     V_0=23_000.0, gamma_0=1.0, psi_0=90.0, altitude=250_000.0, theta_0=0.0,phi_0=0.0,
                     time_elapsed = 2000.0, 
                     spacecraft=False, input_units="imperial", plot_units="imperial", solver="RK45", eqs="betts", plotting=False):
    


    gamma_0 = deg2rad(gamma_0)  # set gamma  deg2rad
    psi_0 = deg2rad(psi_0)
    theta_0 = deg2rad(theta_0)
    phi_0 = deg2rad(phi_0)

    # convert altitude ft 2 m
    altitude = ft2m(altitude) if  input_units=="imperial" else altitude
    r_0 = altitude + EARTH_R
    
    time_span = [0, time_elapsed] # seconds to integrate over

    V_0 = ft2m(V_0) if  input_units == "imperial" else V_0

    initial_conditions = [V_0, gamma_0, psi_0, r_0, theta_0, phi_0] # initial conditions 
 
    time_spans_dense = np.linspace(*time_span, num=200) # sample time for later interpolation of solutions

    if plotting:
        fig, axes = plt.subplots(4,3, figsize=(17,22),  ) # create plot figures

    if eqs == "hicks":
        result = solve_ivp(shuttle_odes_hicks, t_span=time_span, 
                        y0=initial_conditions, args=[get_rho, get_g, spacecraft], 
                        dense_output=True, method=solver, atol=1e-6, rtol=1e-3)

        sol = result.sol(time_spans_dense)
        # dump invalid solutions due unstiff solver sampling , integrating over negative altitudes
        time_spans_dense = time_spans_dense[ np.where(sol[3,:]>EARTH_R)[0]]
        sol = sol[:, np.where(sol[3,:]>EARTH_R)[0]]

        v_sol = sol[0]
        gamma_sol = sol[1]
        psi_sol = sol[2]
        r_sol = sol[3]
        theta_sol = sol[4]
        phi_sol = sol[5]

        h_sol = r_sol - EARTH_R

    elif eqs == "betts":
        initial_conditions[3] -= EARTH_R
        result = solve_ivp(shuttle_odes_betts, t_span=time_span, 
                        y0=initial_conditions, args=[get_rho, get_g, spacecraft], 
                        dense_output=True, method=solver, atol=1e-6, rtol=1e-3)

        sol = result.sol(time_spans_dense)
        # dump invalid solutions due unstiff solver sampling , integrating over negative altitudes
        time_spans_dense = time_spans_dense[ np.where(sol[3,:]>0)[0]]
        sol = sol[:, np.where(sol[3,:]>0)[0]]
        time_sol = time_spans_dense
        v_sol = sol[0]
        gamma_sol = sol[1]
        psi_sol = sol[2]
        h_sol = sol[3]
        theta_sol = sol[4]
        phi_sol = sol[5]

        
        _altitudes  = m2ft(h_sol)/1e3 if plot_units == "imperial" else h_sol/1e3
        
        deccel_sol  = np.gradient(v_sol) / np.gradient(time_spans_dense) / -9.81 # in terms of g
        # atmospheric conditions
        temperatures, pressures, rhos, mus = get_TempPressRhoMu(h_sol) # array of atmospheric conditions at the dense altitude sampling

        L_ref = spacecraft.L_ref
        R_nose = spacecraft.R_nose

        mach_num = v_sol / np.sqrt(1.4 * 287 * temperatures)
        reynolds_num = rhos * v_sol * L_ref / mus
        dynamic_press = 0.5 * rhos * v_sol**2
        stagnation_heat = -1.83e-4 * np.sqrt(rhos/R_nose) * v_sol**3
        stagnation_pressure = rhos * v_sol**2 # assuming Cp=2 for a reentry capsule
        dynamic_energy = v_sol * dynamic_press
        R  = 8.3144626
        # K_i = 0.1235
        # theta = 3055.5556
        # gamma_heat = 1 + ((1.4-1) / (1+(1.4-1)*((theta/temperatures)**2 * (np.exp(theta/temperatures)/(np.exp(theta/temperatures)-1)**2))) )   # calorically imperfect gas gamma 
        # c_p  = 1004.5  * (1 + (gamma_heat-1)/gamma_heat * (theta/temperatures)**2 * np.exp(theta/temperatures)/(np.exp(theta/temperatures)-1)**2)  
        # h_wall = c_p * temperatures + v_sol**2 / 2
        # stagnation_enthalpy = -1*stagnation_heat*np.sqrt(R_nose/(stagnation_pressure/101325)/K_i)  + h_wall
        gamma_heat = v_sol
        c_p = R*1.4 / (1.4-1)
        stagnation_enthalpy = c_p * temperatures + v_sol**2 / 2
       
    if not plotting:
        return sol, time_spans_dense

    if plot_units == "imperial":
        dynamic_press = Pa2lbfsqf(dynamic_press)
        v_sol = m2ft(sol[0]) 
        gamma_sol = rad2deg(sol[1])
        range_sol = m2mi(sol[3])
        stagnation_heat = Wsqm2Btusqft(stagnation_heat)
        stagnation_pressure = Pa2lbfsqf(stagnation_pressure)
        
        axes[0,0].plot(_altitudes, v_sol)
        axes[0,1].plot(_altitudes, deccel_sol,   )
        axes[0,2].plot(_altitudes, dynamic_press,  )
    
        axes[1,0].plot(_altitudes, mach_num,  )
        axes[1,1].plot(_altitudes, reynolds_num,  )
        axes[1,2].plot(_altitudes, stagnation_pressure,  ) ####
        
        axes[2,0].plot(_altitudes, stagnation_enthalpy,  )
        axes[2,1].plot(_altitudes, stagnation_heat,  )
        axes[2,2].plot(_altitudes, time_sol,   )
        
        axes[3,0].plot(_altitudes, v_sol*0,  )
        axes[3,1].plot(_altitudes, dynamic_energy,  )
        axes[3,2].plot(_altitudes, gamma_heat)
    else:
        axes[0,0].plot(v_sol,_altitudes,   )
        axes[0,1].plot(deccel_sol,_altitudes,   )
        axes[0,2].plot(dynamic_press,_altitudes,  )
    
        axes[1,0].plot(mach_num,_altitudes,  )
        axes[1,1].plot(reynolds_num,_altitudes,  )
        axes[1,2].plot(stagnation_pressure,_altitudes,  ) ####
        
        axes[2,0].plot(stagnation_enthalpy,_altitudes,  )
        axes[2,1].plot(stagnation_heat,_altitudes,  )
        axes[2,2].plot(time_sol,_altitudes,   )
        
        axes[3,0].plot(v_sol*0,_altitudes,  )
        axes[3,1].plot(dynamic_energy,_altitudes,  )
        axes[3,2].plot(gamma_heat,_altitudes)    

    if plot_units == "imperial":
        axes[0,0].set_ylabel("Velocity \n$[ft/s]$")
        axes[0,2].set_ylabel("Dynamic Pressure \n$[lbf / ft^{2}]$")
        axes[0,1].set_ylabel("Deceleration $[g]$")
        axes[1,0].set_ylabel("Mach Number")
        axes[1,1].set_ylabel("Reynolds Number")
        axes[1,2].set_ylabel("Stagnation Point \nPressure \n$[lbf/ ft^{2}]$")
        axes[2,0].set_ylabel("Stagnation Point \nEnthalpy \n$[Btu/lbm]$")
        axes[2,1].set_ylabel("Stagnation Point \nHeat Transfer \n$[Btu/s - ft^{2}]$")
        axes[2,2].set_ylabel("Entry Time $[s]$")
        axes[3,0].set_ylabel("Range $[mi]$")
        axes[3,1].set_ylabel("Dynamic Energy\n $[Btu/s-ft^{2}]$")
    else:
        axes[0,0].set_xlabel("Velocity \n$[m/s]$")
        axes[0,2].set_xlabel("Dynamic Pressure \n$[Pa]$")
        axes[0,1].set_xlabel("Deceleration $[g]$")
        axes[1,0].set_xlabel("Mach Number")
        axes[1,1].set_xlabel("Reynolds Number")
        axes[1,2].set_xlabel("Stagnation Point \nPressure \n$[Pa]$")
        axes[2,0].set_xlabel("Stagnation Point \nEnthalpy \n$[J/kg]$")
        axes[2,1].set_xlabel("Stagnation Point \nHeat Transfer \n$[W/m^{2}]$")
        axes[2,2].set_xlabel("Entry Time $[s]$")
        axes[3,0].set_xlabel("Range $[m]$")
        axes[3,1].set_xlabel("Dynamic Energy\n $[W/m^{2}]$")


    for idx,ax in enumerate(axes.reshape(-1)):
        if plot_units == "imperial":
            ax.set_xlabel("Altitude $[10^3 \;ft]$")
            ax.set_xticks(np.arange(0, 300, step=50))
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        else:
            ax.set_ylabel("Altitude $[km]$")
            # ax.xaxis.set_major_locator(plt.MaxNLocator(5))

        
        ax.grid(which="both")
        ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0);
        # ax.legend()

    fig.delaxes(axes[3,2])
    fig.delaxes(axes[3,0])

    # fig.tight_layout()
    plt.subplots_adjust( wspace= 1.05, hspace=0.55, )

    plt.show()

            
    return sol, time_spans_dense
