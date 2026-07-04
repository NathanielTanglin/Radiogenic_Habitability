import matplotlib.pyplot as plt
import numpy
import os
from lib import *
from paths import path

# "eunits" means "Earth units"
# "si" means "standard international units"

G_SI = 6.674 * 1e-11 # SI.
EPSILON_0_SI = 8.82 * 1e-12 # C^2 kg^-1 m^-3 s^2 (SI)
MU_0_SI = 1.256637 * 1e-6 # SI.
EARTH_MAGNETIC_MOMENT_SI = 7.94 * 1e22 # A m^2 (SI). From https://modern-physics.org/earths-dipole-moment/
AU_TO_METER = 149597870691.0; # SI.
EARTH_RADIUS_SI = 6.371 * 1e6 # Meters (SI). 

def dipole_moment_to_B_field_si(dipole_moment_si, radius_si):
    return MU_0_SI*dipole_moment_si/(4*np.pi*np.power(radius_si, 3.0)) # SI unites

# F_AP is a unitless fraction.
def F_AP(stellar_B_si, planet_B_si, stellar_radius_si, planet_semi_major_axis_si):
    beta = (stellar_B_si/planet_B_si)*np.power(stellar_radius_si/planet_semi_major_axis_si, 3.0)
    beta[beta > 1.0] = 1.0

    F_AP = 1.0 - np.power(1.0 - (3.0*np.power(beta, 1.0/3.0))/(2.0+beta), 0.5)

    return F_AP

solar_dynamics = read_vplanet(path("data", "radioisotope_planet_B_F_AP", "sol.sun.forward"))
time_Gyr = solar_dynamics["Time"]/1e9 # Converts from yr to Gyr.
mask = (time_Gyr < 6)
time_Gyr = time_Gyr[mask]

stellar_B_si = 0.0005 # In Tesla
stellar_radius_si = solar_dynamics["Radius"][mask] * EARTH_RADIUS_SI # Converts from eunits to si.

hot_planet_dynamics       = read_vplanet(path("data", 'radioisotope_planet_B_F_AP', 'sol.earth_hot.forward'))
cool_planet_dynamics      = read_vplanet(path("data", 'radioisotope_planet_B_F_AP', 'sol.earth_cool.forward'))
planet_dynamics           = read_vplanet(path("data", 'radioisotope_planet_B_F_AP', 'sol.earth.forward'))

hot_planet_mag_moment_si  = hot_planet_dynamics["MagMom"][mask] * EARTH_MAGNETIC_MOMENT_SI # Converts from eunits to si.
cool_planet_mag_moment_si = cool_planet_dynamics["MagMom"][mask] * EARTH_MAGNETIC_MOMENT_SI # Converts from eunits to si.
planet_mag_moment_si      = planet_dynamics["MagMom"][mask] * EARTH_MAGNETIC_MOMENT_SI # Converts from eunits to si.

hot_planet_B_si           = dipole_moment_to_B_field_si(hot_planet_mag_moment_si, EARTH_RADIUS_SI)
cool_planet_B_si          = dipole_moment_to_B_field_si(cool_planet_mag_moment_si, EARTH_RADIUS_SI)
planet_B_si               = dipole_moment_to_B_field_si(planet_mag_moment_si, EARTH_RADIUS_SI)

hot_planet_F_AP           = F_AP(stellar_B_si, hot_planet_B_si, stellar_radius_si, AU_TO_METER)
cool_planet_F_AP          = F_AP(stellar_B_si, cool_planet_B_si, stellar_radius_si, AU_TO_METER)
planet_F_AP               = F_AP(stellar_B_si, planet_B_si, stellar_radius_si, AU_TO_METER)

# Old incorrect equation.
# (4.0*np.pi*EPSILON_0*3.0*mag_mom*EARTH_MAGNETIC_MOMENT)/(2.0*np.power(planet_radius, 3.0))

# Formatting parameters.
label_fontsize = 14

colors = {
    "hot": (1.0, 0.0, 0.0),
    "fiducial": (0.5, 0.0, 0.5),
    "cool": (0.0, 0.0, 1.0),
}

labels = {
    "hot": r"$\times3$ primordial" + "\nEarth radioisotopes",
    "cool": r"$\times0.3$ primordial" + "\nEarth radioisotopes",
    "fiducial": r"Primordial Earth",
}

fig, ax = plt.subplots(2, 1, sharex = True)
fig.set_size_inches(5, 9)
fig.subplots_adjust(hspace = 0.1)

# Plots the data.
ax[0].plot(time_Gyr, hot_planet_B_si, color = colors["hot"], label = labels["hot"])
ax[0].plot(time_Gyr, planet_B_si, color = colors["fiducial"], label = labels["fiducial"])
ax[0].plot(time_Gyr, cool_planet_B_si, color = colors["cool"], label = labels["cool"])
ax[0].legend()

def order_of(num):
    try:
        return int(math.floor(math.log10(num)))
    except ValueError:
        print(num)

        return 0

ax[0].yaxis.set_major_formatter(tck.FuncFormatter(lambda data, pos: r"$" + str(round(data/(10**(order_of(data))), 2)) + r"\times 10^{" + str(order_of(data)) + r"}$"))

ax[1].plot(time_Gyr, hot_planet_F_AP, color = colors["hot"], label = labels["hot"])
ax[1].plot(time_Gyr, planet_F_AP, color = colors["fiducial"], label = labels["fiducial"])
ax[1].plot(time_Gyr, cool_planet_F_AP, color = colors["cool"], label = labels["cool"])
ax[1].legend()

# Axes labels.
ax[1].set_xlabel("Time [Gyr]", fontsize = label_fontsize)
ax[0].set_ylabel(r"$B_p$ $[\text{Tesla}]$", fontsize = label_fontsize)
ax[1].set_ylabel("Fraction of surface supporting outflow\n" + r"($F_{AP}$)", fontsize = label_fontsize)

save_path = path("plots", "planet_B_F_AP_evolution.png")
fig.savefig(save_path, bbox_inches = "tight", dpi = 300)
plt.close()