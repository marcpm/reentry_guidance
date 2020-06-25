# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from reentrygnc.ballistic import ballistic_odes, run_ballistic_simulation
from reentrygnc.atmos import get_rho, get_g
from reentrygnc.base import  EARTH_R, lbfsqf2Nsqm, Pa2lbfsqf, Shuttle, ft2m, rad2deg
from reentrygnc.shuttle import bank_opt_control, alpha_opt_control, bank_control_NASA, alpha_control_NASA
from reentrygnc.shuttle import run_shuttle_simulation
from scipy.integrate import solve_ivp
from math import pi
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from reentrygnc.coesa import table

import numpy as np


# %%
plt.style.use('Solarize_Light2')
font = {'family' : 'monospace',
        'weight' : 'regular',
        'size'   : 14}

plt.rc('font', **font) 
plt.rc('legend',fontsize=10)


# %%

def shuttle_example(equations="colonel", control_algorithm = "chern-hiu", altitude = 250_000,
                    V_0 = 25_600.0, gamma_0s = [-1.0], psi_0 = 90.0, theta_0 = 0.0, phi_0 = 0.0,         
                    time_lapse = 2000 , L_star=2.0):

    m = 182_986 # lbf
    A_ref = 2690 # ft^2 m^2
    L_ref = 107.5 # ft
    R_nose = 1.0  # ft
    LD_max = L_star

    spacecraft = Shuttle(m, A_ref, R_nose, L_ref, LD_max, 
                                parachute=False, control_algorithm=control_algorithm,       
                                imperial_units=True)

    # altitude = 250_000
    # V_0 = 25_600.0
    # gamma_0s = [-1.0]
    # psi_0 = 90.0

    # theta_0 = 0.0
    # phi_0 = 0.0
    # time_lapse = 1300 # 1900 # optim control / 1200 NASA control
    
    sol, t = run_shuttle_simulation(
                     V_0=V_0, gamma_0s=gamma_0s,  
                      psi_0=psi_0, altitude=altitude, theta_0=theta_0,phi_0=phi_0,
                    time_elapsed=time_lapse,
                    spacecraft=spacecraft,
                     input_units="imperial", solver="RK45", eqs=equations)

    return sol,t       


# %%



# %%



# %%



# %%



# %%



# %%
fig, axes = plt.subplots(1,2, figsize=(14,6))
for control_law in ["chern-hiu", "chern-hao", "gells-hao", "gells-hiu", "gear-optim",]: # "nasa-guidance"
    sol, t = shuttle_example(equations="hicks", control_algorithm=control_law, psi_0=)
    v_sol = sol[0]
    gamma_sol = sol[1]
    psi_sol = sol[2]
    r_sol = sol[3]
    theta_sol = sol[4]
    phi_sol = sol[5]
    h_sol = r_sol - EARTH_R
    axes[0].plot(t, h_sol/1e3, label=control_law)
    axes[1].plot(t, v_sol, label=control_law)

axes[0].set_xlabel("Time [s]"); axes[1].set_xlabel("Time [s]")
axes[0].set_ylabel("Altitude [km]")    ; axes[1].set_ylabel("Velocity [m/s]")    
axes[0].set_title("Hicks DAEs"); axes[1].set_title("Hicks DAEs")
axes[0].legend(); axes[1].legend()
plt.show()


# %%
fig, axes = plt.subplots(1,2, figsize=(14,6))
for control_law in ["chern-hiu", "chern-hao", "gells-hao", "gells-hiu", "gear-optim",]: # "nasa-guidance"
    sol, t = shuttle_example(equations="betts", control_algorithm=control_law, )
    v_sol = sol[0]
    gamma_sol = sol[1]
    psi_sol = sol[2]
    h_sol = sol[3]
    theta_sol = sol[4]
    phi_sol = sol[5]
    # h_sol = r_sol - EARTH_R
    axes[0].plot(t, h_sol/1e3, label=control_law)
    axes[1].plot(t, v_sol, label=control_law)

axes[0].set_xlabel("Time [s]"); axes[1].set_xlabel("Time [s]")
axes[0].set_ylabel("Altitude [km]")    ; axes[1].set_ylabel("Velocity [m/s]")    
axes[0].set_title("T.Betts DAEs"); axes[1].set_title("T.Betts DAEs")
# axes.legend()
plt.show()


# %%
from reentrygnc.shuttle import bank_opt_control, alpha_opt_control, bank_control_NASA, alpha_control_NASA, get_q_wingLE

fig, axes = plt.subplots(2,4, figsize=(18,6))
for control_law in ["gear-optim" ,"nasa-guidance" ]: # "nasa-guidance"
    
    timelapse = 1200 if control_law == "nasa-guidance" else 1410
    sol, t = shuttle_example(equations="betts", control_algorithm=control_law, gamma_0s=[-1.0], psi_0=90, altitude=260_000, time_lapse=timelapse)  #time_lapse=1410 for betts optim
    v_sol = sol[0]
    gamma_sol = sol[1]
    psi_sol = sol[2]
    h_sol = sol[3]
    theta_sol = sol[4]
    phi_sol = sol[5]
    # h_sol = r_sol - EARTH_R
    axes[0,0].plot(t, h_sol/1e3, label=control_law)
    axes[0,1].plot(t, v_sol, label=control_law)
    axes[0,2].plot(phi_sol*180/pi, theta_sol*180/pi, label=control_law)
    axes[0,3].plot(t, gamma_sol*180/pi, label=control_law)
    axes[1,0].plot(t, theta_sol*180/pi, label=control_law)
    axes[1,1].plot(t, phi_sol*180/pi, label=control_law)
    if control_law == "gear-optim":
        axes[1,2].plot(t, alpha_opt_control(t), label=control_law)
        axes[1,3].plot(t, bank_opt_control(t)*180/pi, label=control_law)
    elif control_law == "nasa-guidance":
        axes[1,2].plot(t, alpha_control_NASA(v_sol), label=control_law)
        axes[1,3].plot(t, bank_control_NASA(v_sol)*180/pi, label=control_law)
    

axes[0,0].set_xlabel("Time [s]"); 
axes[0,0].set_ylabel("Altitude [km]")    
# axes[0,0].set_title("T. Betts DAEs");
axes[0,1].set_xlabel("Time [s]"); axes[0,1].set_ylabel("Velocity [m/s]")
axes[0,2].set_xlabel("Longitude $\phi$ [deg]"); axes[0,2].set_ylabel(r"Latitude $\theta$ [deg]")
axes[0,3].set_xlabel("Time [s]"); axes[0,3].set_ylabel("$\gamma$ [deg]")
axes[1,0].set_xlabel("Time [s]"); axes[1,0].set_ylabel(r"Latitude $\theta$ [deg]")
axes[1,1].set_xlabel("Time [s]"); axes[1,1].set_ylabel("Longitude $\phi$ [deg]")
axes[1,2].set_xlabel("Time [s]"); axes[1,2].set_ylabel(r"AOA $\alpha$ [deg]")
axes[1,3].set_xlabel("Time [s]"); axes[1,3].set_ylabel("Bank angle $\sigma$ [deg]")
# axes[1,2].set_xlabel("Time [s]"); axes[1,2].set_ylabel("Time [s]")



for ax in axes.reshape(-1):
    ax.legend(["gear-optim" , "nasa-guidance"])

plt.tight_layout()
plt.show()


# %%
fig, axes = plt.subplots(1,2, figsize=(14,6))
for control_law in ["chern-hiu", "chern-hao", "gells-hao", "gells-hiu", ]: # "nasa-guidance"
    sol, t = shuttle_example(equations="betts", control_algorithm=control_law, time_lapse=2200)
    v_sol = sol[0]
    gamma_sol = sol[1]
    psi_sol = sol[2]
    h_sol = sol[3]
    theta_sol = sol[4]
    phi_sol = sol[5]
  
    axes[0].plot(t, h_sol/1e3, label=control_law)
    axes[1].plot(t, v_sol, label=control_law)

axes[0].set_xlabel("Time [s]"); axes[1].set_xlabel("Time [s]")
axes[0].set_ylabel("Altitude [km]")    ; axes[1].set_ylabel("Velocity [m/s]")    
axes[0].set_title("Hicks DAEs"); axes[1].set_title("T. Betts DAEs")
axes[0].legend(); axes[1].legend()
plt.show()


# %%



# %%
fig, axes = plt.subplots(1,2, figsize=(14,6))
for control_law  in [ "gells-hiu"]:
    for L_star in np.linspace(0, 4, 5):
        sol, t = shuttle_example(equations="betts", control_algorithm=control_law, L_star=L_star)
        v_sol = sol[0]
        gamma_sol = sol[1]
        psi_sol = sol[2]
        h_sol = sol[3]
        theta_sol = sol[4]
        phi_sol = sol[5]
        # h_sol = r_sol - EARTH_R
        axes[0].plot(t, h_sol/1e3, label=control_law)
        axes[1].scatter(L_star, theta_sol[-1]*180/pi, label=control_law)

axes[0].set_xlabel("Time [s]"); axes[1].set_xlabel("Max Lift-Drag Ratio, $L_{*}$")
axes[0].set_ylabel("Altitude [km]")    ; axes[1].set_ylabel(r"Final Latitude $\theta$ [deg]")    
axes[0].set_title("T.Betts DAEs"); axes[1].set_title("T.Betts DAEs")
# axes.legend()
plt.show()


# %%



# %%
sol, t = shuttle_example(equations="optimalbook", control_algorithm="chern-hiu")
v_sol = sol[0]
gamma_sol = sol[1]
psi_sol = sol[2]
h_sol = sol[3]
theta_sol = sol[4]
phi_sol = sol[5]

plt.plot(t, h_sol/1e3)
plt.xlabel("Time [s]")
plt.ylabel("Altitude [km]")
plt.title("Optimal Book equations $\phi$ decoupled")


# %%



# %%



# %%
plt.plot(t, v_sol)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")


# %%
plt.plot(t, gamma_sol/pi*180)
plt.xlabel("Time [s]")
plt.ylabel("$\gamma$ [deg]")


# %%
plt.plot(t, psi_sol/pi*180)
plt.xlabel("Time [s]")
plt.ylabel("$\psi$ [deg]")


# %%
plt.plot(t, phi_sol/pi*180)
plt.xlabel("Time [s]")
plt.ylabel("$\phi$ [deg]")


# %%
plt.plot(t, theta_sol/pi*180)
plt.xlabel("Time [s]")
plt.ylabel(r"$\theta$ [deg]")


# %%
import plotly.graph_objects as go

fig = go.Figure(data=go.Scatter3d(x=theta_sol/pi*180, y=phi_sol/pi*180, z=h_sol,
    marker=dict(
        size=4,
        color=h_sol,
        colorscale='Viridis',
    ),line=dict(
        color='darkblue',
        width=2
    )
))
fig.show()


# %%
import plotly.express as px

fig = px.line_geo(lat=28.5971-theta_sol/pi*180, lon=-80.6833-phi_sol/pi*180, )
# fig.update_geos(fitbounds="locations")
# fig.update_layout(height=300, margin={"r":0,"t":0,"l":0,"b":0})
# fig.update_geos(projection_type="orthographic")
fig.add_scattergeo(lat=[28.5971], lon=[-80.6837], text="Kennedy Space Center Shuttle Landing Facility")
fig.update_geos(projection_type="orthographic",
    resolution=50,
    showcoastlines=True, coastlinecolor="RebeccaPurple",
    showland=True, landcolor="LightGreen",
    showocean=True, oceancolor="LightBlue",
    showlakes=False, lakecolor="Blue",
    showrivers=False, rivercolor="Blue"
)
fig.update_layout(height=300, margin={"r":0,"t":0,"l":0,"b":0})

fig.show()


# %%



# %%
from reentrygnc.shuttle import get_q_wingLE
from reentrygnc.base import m2ft
plt.plot(t, get_q_wingLE(alpha_control_NASA((v_sol))) * 100)
t_ = np.linspace(0, 1600)
plt.plot(t_, get_q_wingLE(alpha_opt_control((t_))) *100 )


# %%



# %%
import plotly.express as px
import plotly.graph_objects as go


thetasol = []
phisol = []
hsol = []
vsol =[]
for control_law in ["gear-optim" ,"nasa-guidance" ]: # "nasa-guidance"
    
    timelapse = 1200 if control_law == "nasa-guidance" else 1410
    sol, t = shuttle_example(equations="betts", control_algorithm=control_law, gamma_0s=[-1.0], psi_0=90, altitude=260_000, time_lapse=timelapse)  #time_lapse=1410 for betts optim
    v_sol = sol[0]
    gamma_sol = sol[1]
    psi_sol = sol[2]
    h_sol = sol[3]
    theta_sol = sol[4]
    phi_sol = sol[5]
    vsol.append(v_sol)
    thetasol.append(theta_sol)
    phisol.append(phi_sol)
    hsol.append(h_sol)


    # fig = px.line_geo(lat=28.5971-theta_sol/pi*180, lon=-80.6833-phi_sol/pi*180, )

# fig.update_geos(fitbounds="locations")
# fig.update_layout(height=300, margin={"r":0,"t":0,"l":0,"b":0})
# fig.update_geos(projection_type="orthographic")

fig = go.Figure()
  
# lat=28.5971-thetasol[0]/pi*180, lon=-80.6833-phisol[0]/pi*180, 


fig.add_trace(go.Scattergeo(lat=28.5971-thetasol[0]/pi*180, lon=-80.6833-phisol[0]/pi*180, name="Optimal Control Non-Linear Programming - Gear (1987) Trajectory", mode="lines", line=dict(width = 3,)), )
fig.add_trace(go.Scattergeo(lat=28.5971-thetasol[1]/pi*180, lon=-80.6833-phisol[1]/pi*180, name="NASA JSC - Harpold (1979) Trajectory", mode="lines", line=dict(width = 3,)), )

fig.add_trace(go.Scattergeo(lat=[28.5971], lon=[-80.6833], name="Kennedy Space Center Shuttle Landing Facility", mode="markers", marker=dict(size=9, color='rgb(255, 153, 0)') ))



fig.update_geos(projection_type="orthographic",
resolution=50,
showcoastlines=True, coastlinecolor="RebeccaPurple",
showland=True, landcolor="LightGreen",
showocean=True, oceancolor="LightBlue",
showlakes=False, lakecolor="Blue",
showrivers=False, rivercolor="Blue"
)
fig.update_layout(height=500, width=800, margin={"r":0,"t":0,"l":0,"b":0}, legend_orientation="h")

fig.show()


# %%
fig = go.Figure()

def gen_color_bar(line_trace):
    """
    Generates a trace which shows a colorbar based on a line plot.

    Relevant issue: https://github.com/plotly/plotly.py/issues/1085
    """
    return go.Scatter3d(
        x=line_trace.x, y=line_trace.y, z=line_trace.z,
        mode="markers",
        marker=go.scatter3d.Marker(
            color=line_trace.line.color,
            colorscale=line_trace.line["colorscale"], # Waiting on https://github.com/plotly/plotly.py/issues/1087
            showscale=line_trace.line.showscale,
            opacity=0.00000000000001 # Make invisible, visible=False disables color bar
        ),
        hoverinfo="none",
        showlegend=False
    )

a = go.Scatter3d(x=thetasol[0]/pi*180, y=phisol[0]/pi*180, z=hsol[0],
    marker=dict(
        size=4,
        color=np.gradient(vsol[0]),
        colorscale='Viridis',
    ),line=dict(
        color='darkblue',
        width=2, showscale=True
    )
)
fig.add_trace(go.Scatter3d(x=thetasol[0]/pi*180, y=phisol[0]/pi*180, z=hsol[0],
    marker=dict(
        size=4,
        color=np.gradient(vsol[0]),
        colorscale='Viridis',
    ),line=dict(
        color='darkblue',
        width=2, showscale=True
    )
))
b = go.Scatter3d(x=thetasol[1]/pi*180, y=phisol[1]/pi*180, z=hsol[1],
    marker=dict(
        size=4,
        color= np.gradient(vsol[1]),
        colorscale='Viridis',
    ),line=dict(
        color='darkblue',
        width=2, 
    )
)
fig.add_trace(go.Scatter3d(x=thetasol[1]/pi*180, y=phisol[1]/pi*180, z=hsol[1],
    marker=dict(
        size=4,
        color= np.gradient(vsol[1]),
        colorscale='Viridis',
    ),line=dict(
        color='darkblue',
        width=2, 
    )
))


fig.add_trace(go.Scatter3d(x=[thetasol[1][-1]/pi*180], y=[phisol[1][-1]/pi*180], z=[hsol[1][-1]],
    marker=dict(
        size=4,
    ),line=dict(
        color='red',
        width=20, 
    )
))




fig.update_layout(showlegend=True, showscale=True, scene = dict(

                    xaxis_title='Latitude [deg]',
                    yaxis_title='Longitude [deg]',
                    zaxis_title='Altitude [m]'),
                    width=700,
                    margin=dict(r=20, b=10, l=10, t=10))
# colorbar_trace = gen_color_bar(a)


# go.FigureWidget(data=[a,colorbar_trace])




fig.show()


# %%
trace=go.Scatter3d(x=rad2deg(thetasol[1]),y=rad2deg(phisol[1]),z=hsol[1],
                   mode='markers', 
                   marker=dict(size=3,color=vsol[1],colorscale='Viridis',colorbar=dict(title="Velocity [m/s]"), showscale=True),)
data=[trace]

layout=go.Layout(scene = dict(
               xaxis_title='Latitude [deg]',
               yaxis_title='Longitude [deg]',
               zaxis_title='Altitude [m]'),
               width=700, height=700,
               margin=dict(r=40, b=40, l=40, t=60))
fig=go.Figure(data=data,layout=layout, )
fig.update_layout(title="NASA JSC - Harpold (1979) Re-entry to TAEM Trajectory",)
fig.show()


# %%
trace=go.Scatter3d(x=rad2deg(thetasol[0]),y=rad2deg(phisol[0]),z=hsol[0],
                   mode='markers', 
                   marker=dict(size=3,color=vsol[0],colorscale='Viridis',colorbar=dict(title="Velocity [m/s]"), showscale=True),)
data=[trace]

layout=go.Layout(scene = dict(
               xaxis_title='Latitude [deg]',
               yaxis_title='Longitude [deg]',
               zaxis_title='Altitude [m]'),
               width=700, height=700,
               margin=dict(r=40, b=40, l=40, t=60))
fig=go.Figure(data=data,layout=layout, )
fig.update_layout(title="Optimal Control Non-Linear Prog. - Gear (1987) Re-entry to TAEM Trajectory",)
fig.show()


# %%
trace=go.Scatter3d(x=rad2deg(thetasol[0]),y=rad2deg(phisol[0]),z=hsol[0],
                   mode='markers', 
                   marker=dict(size=3,color=np.gradient(vsol[0])/-9.80655,colorscale='Jet',colorbar=dict(title="Deceleration [g]"), showscale=True),)
data=[trace]

layout=go.Layout(scene = dict(
               xaxis_title='Latitude [deg]',
               yaxis_title='Longitude [deg]',
               zaxis_title='Altitude [m]'),
               width=700, height=700,
               margin=dict(r=40, b=40, l=40, t=60))
fig=go.Figure(data=data,layout=layout, )
fig.update_layout(title="Optimal Control Non-Linear Prog. - Gear (1987) Re-entry to TAEM Trajectory",)
fig.show()


# %%
trace=go.Scatter3d(x=rad2deg(thetasol[1]),y=rad2deg(phisol[1]),z=hsol[0],
                   mode='markers', 
                   marker=dict(size=3,color=np.gradient(vsol[1])/-9.80655,colorscale='Jet',colorbar=dict(title="Deceleration [g]"), showscale=True),)
data=[trace]

layout=go.Layout(scene = dict(
               xaxis_title='Latitude [deg]',
               yaxis_title='Longitude [deg]',
               zaxis_title='Altitude [m]'),
               width=700, height=700,
               margin=dict(r=40, b=40, l=40, t=60))
fig=go.Figure(data=data,layout=layout, )
fig.update_layout(title="NASA JSC - Harpold (1979) Re-entry to TAEM Trajectory",)
fig.show()


# %%
1+1


# %%


