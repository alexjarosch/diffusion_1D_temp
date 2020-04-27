import numpy as np
import matplotlib.pyplot as plt
import copy
import tkinter as tk

"""A 1D implementation of thermal conductivity through soil and air"""

__author__ = "Alexander H. Jarosch (research@alexj.at)"
__date__ = "2020-01-29 -- 2020-04-16"
__copyright__ = "Copyright (C) 2020 ThetaFrame Solutions"
__license__ = "GNU GPL Version 3"
__version__ = "4.2"

"""
This script is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this script.  If not, see <http://www.gnu.org/licenses/>.
"""

# time constant as we calculate in SI
sec_in_year = 365.25 * 24 * 3600

window = tk.Tk()
window.title("Cave Heat Conduction Model with Snow")
window.geometry('700x600')

# input for depth
lbl_dc = tk.Label(window, text="Depth to Simulate [m]")
lbl_dc.grid(column=0, row=0, sticky="W")
txt_dc = tk.Entry(window, width=10)
txt_dc.insert(tk.END, '100.0')
txt_dc.grid(column=1, row=0)

lbl_zres = tk.Label(window, text="Vertical Resolution N")
lbl_zres.grid(column=0, row=1, sticky="W")
txt_zres = tk.Entry(window, width=10)
txt_zres.insert(tk.END, '51')
txt_zres.grid(column=1, row=1)

# input for jan tmp
lbl_jaT = tk.Label(window, text="January Mean TMP [C]")
lbl_jaT.grid(column=0, row=2, sticky="W")
txt_jaT = tk.Entry(window, width=10)
txt_jaT.insert(tk.END, '-12.0')
txt_jaT.grid(column=1, row=2)

# input for jul tmp
lbl_juT = tk.Label(window, text="July Mean TMP [C]")
lbl_juT.grid(column=0, row=3, sticky="W")
txt_juT = tk.Entry(window, width=10)
txt_juT.insert(tk.END, '8.0')
txt_juT.grid(column=1, row=3)

# input for snow delta_tmp
lbl_sdT = tk.Label(window, text="Snow delta TMP [C]")
lbl_sdT.grid(column=0, row=4, sticky="W")
txt_sdT = tk.Entry(window, width=10)
txt_sdT.insert(tk.END, '0.0')
txt_sdT.grid(column=1, row=4)

# input for soil temp
lbl_iTg = tk.Label(window, text="Initial Ground TMP [C]")
lbl_iTg.grid(column=0, row=5, sticky="W")
txt_iTg = tk.Entry(window, width=10)
txt_iTg.insert(tk.END, '-2.0')
txt_iTg.grid(column=1, row=5)

# input for cave temp
lbl_gHf = tk.Label(window, text="Ground Heat Flux [W/m^2]")
lbl_gHf.grid(column=0, row=6, sticky="W")
txt_gHf = tk.Entry(window, width=10)
txt_gHf.insert(tk.END, '0.05')
txt_gHf.grid(column=1, row=6)

# input for thermal diffusivity
lbl_alphaG = tk.Label(window, text="Thermal Diffusivity Ground [m^2 / s]")
lbl_alphaG.grid(column=0, row=7, sticky="W")
txt_alphaG = tk.Entry(window, width=10)
txt_alphaG.insert(tk.END, '1.0e-6')
txt_alphaG.grid(column=1, row=7)

# input for thermal diffusivity
lbl_kG = tk.Label(window, text="Thermal Conductivity Ground [W / (m K)]")
lbl_kG.grid(column=0, row=8, sticky="W")
txt_kG = tk.Entry(window, width=10)
txt_kG.insert(tk.END, '2.0')
txt_kG.grid(column=1, row=8)

# input for model runtime
lbl_rtime = tk.Label(window, text="Model Runtime [years]")
lbl_rtime.grid(column=0, row=9, sticky="W")
txt_rtime = tk.Entry(window, width=10)
txt_rtime.insert(tk.END, '400')
txt_rtime.grid(column=1, row=9)

# input for plot interval
lbl_intP = tk.Label(window, text="Interval to Plot [years]")
lbl_intP.grid(column=0, row=10, sticky="W")
txt_intP = tk.Entry(window, width=10)
txt_intP.insert(tk.END, '50')
txt_intP.grid(column=1, row=10)

# check if you want to see the points in the line
pltpoints = tk.IntVar()
cb_prtP = tk.Checkbutton(window, text='Plot model points', variable=pltpoints, onvalue=1, offvalue=0)
cb_prtP.grid(column=0, row=11, sticky="W")

# check if you want to see the points in the line
pltintT = tk.IntVar()
cb_prtiT = tk.Checkbutton(window, text='Plot initial TMP profile', variable=pltintT, onvalue=1, offvalue=0)
cb_prtiT.grid(column=0, row=12, sticky="W")

# check if you want to see the points in the line
pltLeg = tk.IntVar()
cb_prtLeg = tk.Checkbutton(window, text='Plot Legend', variable=pltLeg, onvalue=1, offvalue=0)
cb_prtLeg.grid(column=0, row=13, sticky="W")

# save solution
svSol = tk.IntVar()
cb_svSol = tk.Checkbutton(window, text='Save Final Solution', variable=svSol, onvalue=1, offvalue=0)
cb_svSol.grid(column=0, row=14, sticky="W")
txt_svSol = tk.Entry(window, width=20)
txt_svSol.insert(tk.END, 'solution.txt')
txt_svSol.grid(column=1, row=14)

# read solution
rdSol = tk.IntVar()
cb_rdSol = tk.Checkbutton(window, text='Read Initial TMP Profile', variable=rdSol, onvalue=1, offvalue=0)
cb_rdSol.grid(column=0, row=15, sticky="W")
txt_rdSol = tk.Entry(window, width=20)
txt_rdSol.insert(tk.END, 'solution_in.txt')
txt_rdSol.grid(column=1, row=15)


# forcing with annual cycle
def T_surf(t):
    jan_mean_tmp = float(txt_jaT.get())      # deg C
    jul_mean_tmp = float(txt_juT.get())      # deg C
    TMP_amp = (np.abs(jan_mean_tmp) + np.abs(jul_mean_tmp)) / 2
    TMP_now = jan_mean_tmp + TMP_amp * (1 + np.cos(2 * np.pi * 1/sec_in_year * t - np.pi))
    return TMP_now


def T_snow(t):
    snow_amp = float(txt_sdT.get())             # deg C
    TMP_snow = snow_amp/2 * (1 + np.cos(2 * np.pi * 2/sec_in_year * t))
    # block out the summer period, which is between 3 and 9 (months)
    # TMP_snow[np.logical_and(np.mod(t, sec_in_year) >= 7889400.0, np.mod(t, sec_in_year) <= 23668200.0)] = 0.0
    if np.mod(t, sec_in_year) >= 7889400.0 and np.mod(t, sec_in_year) <= 23668200.0:
        TMP_snow = 0.0

    return TMP_snow


def updateS():
    dt_range = sec_in_year/1500
    t_range = np.arange(0, 2*sec_in_year+dt_range, dt_range)
    TMP_snow_plot = np.zeros(len(t_range))
    for i in range(len(t_range)):
        TMP_snow_plot[i] = T_snow(t_range[i])

    # do tmp info update
    amean_T_a = np.mean(T_surf(t_range))
    amean_T_eff = np.mean(T_surf(t_range) + TMP_snow_plot)
    amean_T_txt = "%0.2f C" % amean_T_a
    amean_Te_txt = "%0.2f C" % amean_T_eff
    lbl_meanT2.configure(text=amean_T_txt)
    lbl_meanT2.update()
    lbl_meanT2e.configure(text=amean_Te_txt)
    lbl_meanT2e.update()

    # do vert res update
    [zga, dz_is] = np.linspace(0, float(txt_dc.get()), int(txt_zres.get()), retstep=True)
    dz_txt = "%0.4f m" % dz_is
    lbl_dz2.configure(text=dz_txt)
    lbl_dz2.update()


def updateSP():
    dt_range = sec_in_year/1500
    t_range = np.arange(0, 2*sec_in_year+dt_range, dt_range)
    TMP_snow_plot = np.zeros(len(t_range))
    for i in range(len(t_range)):
        TMP_snow_plot[i] = T_snow(t_range[i])

    # do tmp info update
    amean_T_a = np.mean(T_surf(t_range))
    amean_T_eff = np.mean(T_surf(t_range) + TMP_snow_plot)
    amean_T_txt = "%0.2f C" % amean_T_a
    amean_Te_txt = "%0.2f C" % amean_T_eff
    lbl_meanT2.configure(text=amean_T_txt)
    lbl_meanT2.update()
    lbl_meanT2e.configure(text=amean_Te_txt)
    lbl_meanT2e.update()

    # do vert res update
    [zga, dz_is] = np.linspace(0, float(txt_dc.get()), int(txt_zres.get()), retstep=True)
    dz_txt = "%0.4f m" % dz_is
    lbl_dz2.configure(text=dz_txt)
    lbl_dz2.update()

    plt.figure()
    plt.title("Surface TMP forcing")
    plt.plot(t_range/sec_in_year*12, T_surf(t_range), 'b', label="air TMP")
    plt.plot(t_range/sec_in_year*12, TMP_snow_plot, 'r', label="snow deltaTMP")
    plt.plot(t_range/sec_in_year*12, T_surf(t_range) + TMP_snow_plot, 'g', label="effective TMP")
    plt.xlabel("Time in months")
    plt.ylabel("TMP in C")
    plt.legend()

    plt.show()


def clicked():

    updateS()

    """ PARAMETER SPACE """

    # thermal diffusivities in SI units
    alpha_ground = float(txt_alphaG.get())      # m^2 / s
    k_ground = float(txt_kG.get())         # m^2 / s

    # geometry
    depth_to_cave = float(txt_dc.get())    # m
    Nz_ground = int(txt_zres.get())        # m

    # initial temperatures
    t_init_ground = float(txt_iTg.get())      # deg C
    # ground heat Flux
    g_heat_flux = float(txt_gHf.get())        # W/m2

    # time domain
    t_save = int(txt_intP.get())       # years
    t_total = int(txt_rtime.get())      # years
    t_solve = np.arange(0, t_total+t_save, t_save)
    Nt = len(t_solve)

    """ Model Code """

    # space domain
    [z1, dz1] = np.linspace(0, depth_to_cave, Nz_ground, retstep=True)
    # add a point to do the bottom BC which we don't plot
    [z, dz] = np.linspace(0, depth_to_cave+2*dz1, Nz_ground+2, retstep=True)
    Nz = len(z)
    z_total = z[-1]

    # diffusivity
    D = alpha_ground

    # calculate the TMP gradient for the lower boundary
    deltaT_BC = (g_heat_flux * dz)/k_ground
    deltaT_grad = g_heat_flux/k_ground

    T_data = np.zeros((Nt+1, Nz))

    if rdSol.get() == 1:
        T_init_read = np.loadtxt(txt_rdSol.get(), delimiter=',')
        T_data[0, :] = T_init_read
    else:
        # set initial TMP conditions
        T_data[0, :] = t_init_ground + deltaT_grad*z

    # lower boundary condition
    # T_data[0, -2] = t_init_ground + deltaT_BC

    # define FDM index array for space
    # make i index array
    index_i = np.arange(1, Nz-1)
    # make i plus array
    index_i_plus = np.arange(2, Nz)
    # make i minus array
    index_i_minus = np.arange(0, Nz-2)

    t_save_sec = t_save * sec_in_year

    # set initial profile
    T_data[0, 0] = T_surf(0) + T_snow(0)

    # overall time in seconds
    t = 0
    # time stepping look
    for n in range(Nt-1):
        t_stab = 0
        u = copy.deepcopy(T_data[n, :])
        while t_stab < t_save_sec:
            dt_stab = 0.5 * dz**2 / np.max(D)
            dt_use = min(dt_stab, t_save_sec - t_stab)
            t_stab = t_stab + dt_use
            t = t + dt_use

            space_diff = u[index_i_plus] - 2*u[index_i] + u[index_i_minus]
            # boundary conditions and time step
            u[index_i] = u[index_i] + D*(dt_use/dz**2)*space_diff
            # upper BC
            u[0] = T_surf(t) + T_snow(t)
            # lower BC
            # u[Nz-1] = u[Nz-2]
            u[Nz-1] = u[Nz-2]
            u[Nz-2] = u[Nz-1] + deltaT_BC

            time_info = "%0.4f years" % (t/sec_in_year)
            lbl_time2.configure(text=time_info)
            lbl_time2.update()

        T_data[n+1, :] = copy.deepcopy(u)

    # save solution

    if svSol.get() == 1:
        np.savetxt(txt_svSol.get(), T_data[-2,:], delimiter=',')


    # Visuals

    if pltpoints.get() == 0:
        mstring = ''
    elif pltpoints.get() == 1:
        mstring = '.'

    if pltintT.get() == 0:
        pltrange = range(1, len(t_solve))
    elif pltintT.get() == 1:
        pltrange = range(0, len(t_solve))

    plt.figure()
    plt.title("TMP with depth")
    for pidx in pltrange:
        plt.plot(T_data[pidx, :-2], z[:-2], marker=mstring, label=("Year %d" % t_solve[pidx]))

    plt.plot([0, 0], [0, z_total], '--k')
    plt.ylim(z_total, 0)
    if pltLeg.get() == 1:
        plt.legend()

    plt.xlabel("TMP in C")
    plt.ylabel("Depth in m")

    plt.show()


lbl_info = tk.Label(window, text="Model Setup Info:")
lbl_info.grid(column=1, row=16)

lbl_dz = tk.Label(window, text="Vertical Resolution dz:")
lbl_dz.grid(column=0, row=17, sticky="W")
lbl_dz2 = tk.Label(window, text="%0.4f m" % (2))
lbl_dz2.grid(column=1, row=17)

amean_T = np.mean([float(txt_jaT.get()), float(txt_juT.get())])
lbl_meanT = tk.Label(window, text="Annual Mean Air TMP:")
lbl_meanT.grid(column=0, row=18, sticky="W")
lbl_meanT2 = tk.Label(window, text="%0.2f C" % amean_T)
lbl_meanT2.grid(column=1, row=18)

lbl_meanTe = tk.Label(window, text="Annual Mean Effective TMP:")
lbl_meanTe.grid(column=0, row=19, sticky="W")
lbl_meanT2e = tk.Label(window, text="%0.2f C" % amean_T)
lbl_meanT2e.grid(column=1, row=19)

btn = tk.Button(window, text="Update Setup / Plot Input", command=updateSP)
btn.grid(column=1, row=20)

lbl_time = tk.Label(window, text="Model Time:")
lbl_time.grid(column=0, row=21, sticky="W")
lbl_time2 = tk.Label(window, text="%0.4f years" % (0))
lbl_time2.grid(column=1, row=21)

btn = tk.Button(window, text="Update Setup and Run Model", command=clicked)
btn.grid(column=1, row=22)

lbl_space = tk.Label(window, text="")
lbl_space.grid(column=0, row=23)

lbl_space = tk.Label(window, text="")
lbl_space.grid(column=0, row=24)

lbl_CR = tk.Label(window, text="Ver. 4.2 - Copyright (C) 2020 ThetaFrame Solutions, GPLv3")
lbl_CR.grid(column=1, row=25)

window.mainloop()
