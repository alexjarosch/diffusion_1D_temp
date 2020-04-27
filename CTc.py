import numpy as np
import matplotlib.pyplot as plt
import copy
import tkinter as tk

"""A 1D implementation of thermal conductivity through soil and air"""

__author__ = "Alexander H. Jarosch (research@alexj.at)"
__date__ = "2020-01-29 --"
__copyright__ = "Copyright (C) 2020 Alexander H. Jarosch"
__license__ = "GNU GPL Version 3"
__version__ = "0.5"

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
window.title("Cave Heat Conduction Model")
window.geometry('600x500')

# input for depth
lbl_dc = tk.Label(window, text="Depth to Cave [m]")
lbl_dc.grid(column=0, row=0, sticky="W")
txt_dc = tk.Entry(window, width=10)
txt_dc.insert(tk.END, '60.0')
txt_dc.grid(column=1, row=0)

# input for cave
lbl_hc = tk.Label(window, text="Height of Cave [m]")
lbl_hc.grid(column=0, row=1, sticky="W")
txt_hc = tk.Entry(window, width=10)
txt_hc.insert(tk.END, '6.0')
txt_hc.grid(column=1, row=1)

lbl_zres = tk.Label(window, text="Vertical Resolution N")
lbl_zres.grid(column=0, row=2, sticky="W")
txt_zres = tk.Entry(window, width=10)
txt_zres.insert(tk.END, '21')
txt_zres.grid(column=1, row=2)

# input for jan tmp
lbl_jaT = tk.Label(window, text="January Mean TMP [C]")
lbl_jaT.grid(column=0, row=3, sticky="W")
txt_jaT = tk.Entry(window, width=10)
txt_jaT.insert(tk.END, '-4.0')
txt_jaT.grid(column=1, row=3)

# input for jul tmp
lbl_juT = tk.Label(window, text="July Mean TMP [C]")
lbl_juT.grid(column=0, row=4, sticky="W")
txt_juT = tk.Entry(window, width=10)
txt_juT.insert(tk.END, '18.0')
txt_juT.grid(column=1, row=4)

# input for soil temp
lbl_iTg = tk.Label(window, text="Initial Ground TMP [C]")
lbl_iTg.grid(column=0, row=5, sticky="W")
txt_iTg = tk.Entry(window, width=10)
txt_iTg.insert(tk.END, '0')
txt_iTg.grid(column=1, row=5)

# input for cave temp
lbl_iTc = tk.Label(window, text="Initial Cave TMP [C]")
lbl_iTc.grid(column=0, row=6, sticky="W")
txt_iTc = tk.Entry(window, width=10)
txt_iTc.insert(tk.END, '-2')
txt_iTc.grid(column=1, row=6)

# input for thermal diffusivity
lbl_alphaG = tk.Label(window, text="Thermal Diffusivity Ground [m^2 / s]")
lbl_alphaG.grid(column=0, row=7, sticky="W")
txt_alphaG = tk.Entry(window, width=10)
txt_alphaG.insert(tk.END, ' 3.0e-6')
txt_alphaG.grid(column=1, row=7)

# input for thermal diffusivity
lbl_alphaA = tk.Label(window, text="Thermal Diffusivity Air [m^2 / s]")
lbl_alphaA.grid(column=0, row=8, sticky="W")
txt_alphaA = tk.Entry(window, width=10)
txt_alphaA.insert(tk.END, '1.96e-5')
txt_alphaA.grid(column=1, row=8)

# input for model runtime
lbl_rtime = tk.Label(window, text="Model Runtime [years]")
lbl_rtime.grid(column=0, row=9, sticky="W")
txt_rtime = tk.Entry(window, width=10)
txt_rtime.insert(tk.END, '10')
txt_rtime.grid(column=1, row=9)

# input for plot interval
lbl_intP = tk.Label(window, text="Interval to Plot [years]")
lbl_intP.grid(column=0, row=10, sticky="W")
txt_intP = tk.Entry(window, width=10)
txt_intP.insert(tk.END, '1')
txt_intP.grid(column=1, row=10)


# forcing with annual cycle
def T_surf(t):
    jan_mean_tmp = float(txt_jaT.get())      # deg C
    jul_mean_tmp = float(txt_juT.get())      # deg C
    TMP_amp = (np.abs(jan_mean_tmp) + np.abs(jul_mean_tmp)) / 2
    TMP_now = jan_mean_tmp + TMP_amp * (1 + np.cos(2 * np.pi * 1/sec_in_year * t - np.pi))
    return TMP_now


def updateS():
    # do tmp info update
    amean_T_a = np.mean([float(txt_jaT.get()), float(txt_juT.get())])
    amean_T_txt = "%0.4f C" % amean_T_a
    lbl_meanT2.configure(text=amean_T_txt)
    lbl_meanT2.update()

    # do vert res update
    [zga, dz_is] = np.linspace(0, float(txt_dc.get()), int(txt_zres.get()), retstep=True)
    z_cave_is = float(txt_dc.get()) + np.arange(dz_is, float(txt_hc.get())+dz_is, dz_is)
    Nz_cave = len(z_cave_is)
    dz_txt = "%0.4f m" % dz_is
    Nz_cave_txt = "%d" % Nz_cave
    lbl_dz2.configure(text=dz_txt)
    lbl_dz2.update()
    lbl_caveP2.configure(text=Nz_cave_txt)
    lbl_caveP2.update()


def clicked():

    updateS()

    """ PARAMETER SPACE """

    # thermal diffusivities in SI units
    alpha_ground = float(txt_alphaG.get())      # m^2 / s
    alpha_air = float(txt_alphaA.get())         # m^2 / s

    # geometry
    depth_to_cave = float(txt_dc.get())    # m
    height_of_cave = float(txt_hc.get())   # m
    Nz_ground = int(txt_zres.get())        # m

    # initial temperatures
    t_init_ground = float(txt_iTg.get())      # deg C
    t_init_caveair = float(txt_iTc.get())     # deg C

    # time domain
    t_save = int(txt_intP.get())       # years
    t_total = int(txt_rtime.get())      # years
    t_solve = np.arange(0, t_total+t_save, t_save)
    Nt = len(t_solve)

    """ Model Code """

    # space domain
    [z_ground, dz] = np.linspace(0, depth_to_cave, Nz_ground, retstep=True)
    z_cave = depth_to_cave + np.arange(dz, height_of_cave+dz, dz)
    z = np.hstack([z_ground, z_cave])
    Nz = len(z)
    z_total = z[-1]

    # diffusivity
    # find the soil index
    s_idx = np.where(z <= depth_to_cave)
    c_idx = np.where(z > depth_to_cave)
    D = np.zeros(len(z))
    D[s_idx] = alpha_ground
    D[c_idx] = alpha_air

    T_data = np.zeros((Nt+1, Nz))

    # set initial TMP conditions
    T_data[:, s_idx] = t_init_ground
    T_data[:, c_idx] = t_init_caveair

    # define FDM index array for space
    # make i index array
    index_i = np.arange(1, Nz-1)
    # make i plus array
    index_i_plus = np.arange(2, Nz)
    # make i minus array
    index_i_minus = np.arange(0, Nz-2)

    t_save_sec = t_save * sec_in_year

    # set initial profile
    T_data[0, 0] = T_surf(0)

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
            u[index_i] = u[index_i] + D[index_i]*(dt_use/dz**2)*space_diff
            u[0] = T_surf(t)
            u[Nz-1] = u[Nz-2]

            time_info = "%0.4f years" % (t/sec_in_year)
            lbl_time2.configure(text=time_info)
            lbl_time2.update()

        T_data[n+1, :] = copy.deepcopy(u)

    """ Visuals """

    plt.figure()
    plt.title("Surface TMP forcing")
    dt_range = sec_in_year/48
    t_range = np.arange(0, sec_in_year+dt_range, dt_range)
    plt.plot(t_range/sec_in_year, T_surf(t_range))
    plt.xlabel("Time in years")
    plt.ylabel("TMP in C")

    plt.figure()
    plt.title("TMP with depth")
    for pidx in range(len(t_solve)):
        plt.plot(T_data[pidx, :], z, marker='.', label=("Year %d" % t_solve[pidx]))

    plt.plot(plt.xlim(), [depth_to_cave, depth_to_cave], '--k', label="cave roof")
    plt.ylim(z_total, 0)
    plt.legend()
    plt.xlabel("TMP in C")
    plt.ylabel("Depth in m")

    plt.show()


lbl_info = tk.Label(window, text="Model Setup Info:")
lbl_info.grid(column=1, row=11)

lbl_dz = tk.Label(window, text="Vertical Resolution dz:")
lbl_dz.grid(column=0, row=12, sticky="W")
lbl_dz2 = tk.Label(window, text="%0.4f m" % (3))
lbl_dz2.grid(column=1, row=12)

amean_T = np.mean([float(txt_jaT.get()), float(txt_juT.get())])
lbl_meanT = tk.Label(window, text="Annual Mean TMP:")
lbl_meanT.grid(column=0, row=13, sticky="W")
lbl_meanT2 = tk.Label(window, text="%0.4f C" % amean_T)
lbl_meanT2.grid(column=1, row=13)

lbl_caveP = tk.Label(window, text="Points in Cave:")
lbl_caveP.grid(column=0, row=14, sticky="W")
lbl_caveP2 = tk.Label(window, text="%d" % 2)
lbl_caveP2.grid(column=1, row=14)

btn = tk.Button(window, text="Update Setup", command=updateS)
btn.grid(column=1, row=15)

btn = tk.Button(window, text="Update Setup and Run Model", command=clicked)
btn.grid(column=1, row=17)

lbl_time = tk.Label(window, text="Model Time:")
lbl_time.grid(column=0, row=16, sticky="W")
lbl_time2 = tk.Label(window, text="%0.4f years" % (0))
lbl_time2.grid(column=1, row=16)

lbl_space = tk.Label(window, text="")
lbl_space.grid(column=0, row=18)

lbl_space = tk.Label(window, text="")
lbl_space.grid(column=0, row=19)

lbl_CR = tk.Label(window, text="Ver. 1.0 - Copyright (C) 2020 Alexander H. Jarosch, GPLv3")
lbl_CR.grid(column=1, row=20)

window.mainloop()
