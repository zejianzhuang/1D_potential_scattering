

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


my_setting = {
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,     # use inline math for ticks  
        "text.latex.preamble": "\n".join([
            r"\usepackage[T1]{fontenc}",
            r"\usepackage{newtxtext}",
            r"\usepackage{newtxmath}"]),
        "axes.labelsize": 16.0,
        "axes.titlesize": 16.0,
        "xtick.labelsize": 14.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "ytick.labelsize": 14.0,
        "legend.fontsize": 14.0,
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True
    }
plt.rcParams.update(my_setting)


df0 = pd.read_csv("./output/denominator_of_T_with_E_less_than_0.csv")

fig0, ax0 = plt.subplots()

ax0.plot(df0.iloc[:, 0], df0.iloc[:, 1], label=r"$T^{-1}$")
ax0.set_xlabel(r"$E$ [eV]")
ax0.set_ylabel(r"Denominator of $T$")
ax0.legend()
fig0.savefig("../figure/deno_T.png")
#fig0.savefig("./figure/deno_T.pdf")

df1 = pd.read_csv("./output/pole_structure.csv")
r = df1.iloc[:, 1]
bs = df1.iloc[:, 0]
r = np.array(list(map(complex, r) ) )
bs = np.array(list(map(complex, bs) ) )
bs_in_k = np.array(list(map(complex, df1.iloc[:, 2]) ) ) * 1e-3
r_in_k = np.array(list(map(complex, df1.iloc[:, 3]) ) ) * 1e-3


color = ["#A80326", "#FDB96B"]
fig1, ax1 = plt.subplots(1, 2, figsize=(9, 3.5), layout='constrained')
ax1[0].scatter([r.real, r.real], [r.imag, -r.imag], c=color[1], marker='^', label="Resonace")
ax1[0].scatter(bs.real, bs.imag, c=color[0], marker='s', label="Bound state")
ax1[0].set_xlabel(r"Real part of $E$ [eV]")
ax1[0].set_ylabel(r"Image part of $E$ [eV]")
ax1[0].set_title("Pole in complex energy plane")

#fig1.savefig("./figure/pole_stru.pdf")

# Plot a position of a pole in k plane
ax1[1].scatter(bs_in_k.real, bs_in_k.imag, c=color[0], marker='s', label="Bound state")
ax1[1].scatter([r_in_k.real, -r_in_k.real], [r_in_k.imag, r_in_k.imag], c=color[1], marker='^', label="Resonance")
ax1[1].set_xlabel(r"Real part of $k$ [KeV]")
ax1[1].set_ylabel(r"Image part of $k$ [KeV]")
ax1[1].set_title(r"Pole in complex $k$ plane")

l, label = ax1[0].get_legend_handles_labels()
fig1.legend(l, label, loc="center right")
fig1.savefig("../figure/pole_stru.png")
#fig1.savefig("./figure/pole_stru.pdf")


# Diagram for trans pro
df2 = pd.read_csv("./output/trans_pro.csv")
fig2, ax2 = plt.subplots()
ax2.plot(df2.iloc[:, 0], df2.iloc[:, 1], label=r"$|T|^2$")
ax2.set_xlabel(r"$E$ [eV]")
ax2.set_ylabel("Transmission probability")
ax2.legend()
fig2.savefig("../figure/trans_pro.png")
#fig2.savefig("./figure/trans_pro.pdf")