import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import os
import numpy as np
import pandas as pd
import math
import planet
from planet import Planet
import scipy.interpolate as sci

# Unused.
def sigfigs(val, figs):
    order = math.floor(math.log10(val))

    return round(val, -(order + 1) + figs)

def read_vplanet(filepath, output_options=None, usecols=None):
    if usecols == None:
        usecols = output_options

    if output_options == None:
        path_list = filepath.split(os.sep)

        directory = os.sep.join(path_list[:-1])

        body_file_name = path_list[-1].split('.')[-2] + '.in'

        with open(os.path.join(directory, body_file_name)) as body_file:
            body = body_file.read()

            output_options = body.split('saOutputOrder')[-1].replace('$', '').replace('\n', '').replace('-', '').strip()

            # Sequentially removes unneeded space.
            while output_options.find('  ') != -1:
                output_options = output_options.replace('  ', ' ')

        output_options = output_options.split(' ')

    return pd.read_csv(filepath, sep = ' ', names = output_options, index_col = False, usecols=usecols)

def solar_planet(planet_dynamics, stellar_dynamics, semi_major_axis):
    #semi_major_axis = 1.49598 * 1e11 # In meters.
    solar_mean_surface_B = 1e-4 # In Tesla. From https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html.
    solar_radius = 6.957e8 # Volumetric mean radius in meters. From https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html

    return Planet(solar_mean_surface_B, solar_radius, semi_major_axis, planet_dynamics, stellar_dynamics)

def M_dwarf_planet(planet_dynamics, stellar_dynamics, semi_major_axis):
    #semi_major_axis = 2.89308763e9 # In meters.
    
    #7.08658848e9 # In meters. # 1.49598 * 1e11 # In meters.
    M_dwarf_mean_surface_B = (1.5*1e3) * 1e-4 # In Tesla.
    R_EARTH = 6.371 * 1e6 # In SI base units.
    M_dwarf_radius = 12.3 * R_EARTH # 6.957e8 # Volumetric mean radius in meters. From https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html

    return Planet(stellar_mean_surface_B=M_dwarf_mean_surface_B, stellar_radius=M_dwarf_radius, planet_star_distance=semi_major_axis, planet_dynamics=planet_dynamics, stellar_dynamics=stellar_dynamics)

def solar_metric(planet_dynamics, stellar_dynamics, semi_major_axis):
    return solar_planet(planet_dynamics, stellar_dynamics, semi_major_axis).integrate()

def M_dwarf_metric(planet_dynamics, stellar_dynamics, semi_major_axis):
    return M_dwarf_planet(planet_dynamics, stellar_dynamics, semi_major_axis).integrate()

semi_major_axis = 1.49598 * 1e11 # In meters.
solar_mean_surface_B = 1/1e4 # In Tesla. From https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html.
solar_radius = 6.957e8 # Volumetric mean radius in meters. From https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
EARTH_FLUX = 1361.87527448 # In W/m^2

# Trappist-1 HZ limits.
trappist_1_flux_range = np.linspace(0.25, 1.5, 2) * EARTH_FLUX

# Central optimistic habitable zone. From https://www.astroexplorer.org/details/apjaa34dff1
trappist_1_central_HZ = 0.6 * EARTH_FLUX

kv_central_HZ = 0.6 * EARTH_FLUX

# Solar HZ limits.
sol_flux_range = np.linspace(0.3, 2.75, 2) * EARTH_FLUX

def get_a(flux, luminosity):
    return np.sqrt(luminosity / (4*np.pi*flux))

def get_a_limits(flux, stellar_dynamics):
    SOLAR_LUMINOSITY = (3.83*10**26) # SI base units.

    luminosity = np.array(stellar_dynamics['Luminosity']) * SOLAR_LUMINOSITY

    if type(flux) == float or type(flux) == int:
        return get_a(flux, luminosity)
    else:
        return [get_a(f, luminosity) for f in flux] # SI base units.

INNER = 1
OUTER = 0

# Bitwise flags
ATM = 1
WATER = 2
DESSIC_TIME = 4
PERCENT = 8

def build_contour_VPLanet(save_name = '', directory = 'Parameter_Sweep', planet_file_name = 'planet.in', mode = ATM):
    print('Plotting parameter sweep data from directory {}'.format(directory))

    x = np.zeros((30, 30))
    y = np.zeros((30, 30))
    z = np.zeros((30, 30))

    files = os.listdir(directory)

    unit_label = str()

    for (i, subdirectory) in enumerate(files):
        vpl_file = os.path.join(directory, subdirectory)

        terms = subdirectory.split('_')

        row = int(terms[1].replace('K', ''))
        column = int(terms[2].replace('ThU', ''))

        if i % 50 == 0:
            print('> Fetching data for runs {}-{}...'.format(i+1, i+50))

        planet_dynamics = read_vplanet(os.path.join(vpl_file, planet_file_name))

        initial_40K = planet_dynamics['40KNumCore'].iloc[0] / planet.EARTH_40K_CORE
        initial_232Th = planet_dynamics['232ThNumMan'].iloc[0] / planet.EARTH_232TH_MANTLE

        # Simple bitwise mapping to choose the mass type.
        mass_type = {
            WATER: 'SurfWaterMass',
            ATM: 'EnvelopeMass'
        }[(mode & WATER) | (mode & ATM)]
        
        mass = planet_dynamics[mass_type]

        metric = 0.0

        # DESSIC_TIME and PERCENT are mutually exclusive flags.
        assert (mode & DESSIC_TIME) | (mode & PERCENT) != (DESSIC_TIME | PERCENT)

        time = planet_dynamics['Time']

        if (mode & DESSIC_TIME):
            # Gets the first time snap shot where the mass is gone (equal to zero).
            metric = round(time[mass == 0.0].iloc[0] / 1e9, 2)
            unit_label = 'Time Until Ocean Evaporated [Gyr]' if mass_type == 'SurfWaterMass' else 'Time Until Atmosphere Lost [Gyr]'
        elif (mode & PERCENT):
            metric = 100.0 * (mass.iloc[0] - mass[time <= 6.025300e9].iloc[-1]) / mass.iloc[0]
            unit_label = '% Water Loss' if mass_type == 'SurfWaterMass' else '% Atmospheric Mass Loss'
        else: # Just mass loss.
            metric = mass.iloc[0] - mass[time <= 6.025300e9].iloc[-1]
            unit_label = 'Water Loss [TO]' if mass_type == 'SurfWaterMass' else 'Atmospheric Mass Loss [kg]' 

        # Row -> K, column -> Th, U
        x[row][column] = initial_40K
        y[row][column] = initial_232Th
        z[row][column] = metric

    print('Building plot...'.format(i+1))

    fig, axis = plt.subplots()
    fig.set_size_inches(6, 6)

    im = axis.contourf(x, y, z, levels = 72)
    axis.set_ylabel('$^{232}$Th, $^{235}$U, $^{238}$U initial mantle num\n[primordial Earth abundances]', fontsize = 12)
    axis.set_xlabel('$^{40}$K initial core num\n[primordial Earth abundances]', fontsize = 12)
    axis.xaxis.set_major_locator(tck.LinearLocator(10))
    axis.yaxis.set_major_locator(tck.LinearLocator(10))
    axis.xaxis.set_major_formatter(tck.FuncFormatter(lambda x,pos: '{}'.format(round(x, 2))))
    axis.yaxis.set_major_formatter(tck.FuncFormatter(lambda x,pos: '{}'.format(round(x, 2))))

    # Marks Earth on the contour.
    axis.scatter([1.0], [1.0], s = 30, marker = 'x', color = 'red')
    axis.text(0.9, 1.2, "Earth", fontsize = 12, color = 'red')

    cbar = fig.colorbar(mappable=im, label=unit_label, orientation = 'horizontal', location = 'top')
    cbar.ax.xaxis.set_major_formatter(tck.FuncFormatter(lambda x,pos: '{}'.format(round(x, 2))))

    interp = sci.RegularGridInterpolator((x.T[0], y[0]), z)((1, 1))

    cbar.add_lines([interp], colors = ['red'], linewidths = [2])

    if save_name != '':
        fig.savefig(save_name, dpi = 300, bbox_inches = 'tight')
    
    print('Plot built!')

def get_param_pos(x, y, p):
    X = x[int(len(x) * p)]
    Y = y[int(len(y) * p)]

    return (X, Y)

def plot_thermal_evolution(save_name, planet_dynamics):
    (fig, axes) = plt.subplots(3, 1, sharex = True)

    fig.set_size_inches(4, 10)
    fig.subplots_adjust(hspace=0)

    axes[0].plot(planet_dynamics['Time'], planet_dynamics['40KNumCore'], label = '$^{40}$K (core)')
    axes[0].plot(planet_dynamics['Time'], planet_dynamics['232ThNumMan'], label = '$^{232}$Th (mantle)')
    axes[0].plot(planet_dynamics['Time'], planet_dynamics['235UNumMan'], label = '$^{235}$U (mantle)')
    axes[0].plot(planet_dynamics['Time'], planet_dynamics['238UNumMan'], label = '$^{238}$U (mantle)')

    (x_min, x_max) = axes[0].get_xlim()
    (y_min, y_max) = axes[0].get_ylim()

    x_scale = x_max-x_min
    y_scale = y_max-y_min

    (x, y) = get_param_pos(planet_dynamics['Time'], planet_dynamics['40KNumCore'], 0.1)
    axes[0].text(x, y, '$^{40}$K (core)')

    (x, y) = get_param_pos(planet_dynamics['Time'], planet_dynamics['232ThNumMan'], 0.1)
    axes[0].text(x, y + y_scale*0.01, '$^{232}$Th (mantle)')

    (x, y) = get_param_pos(planet_dynamics['Time'], planet_dynamics['238UNumMan'], 0.2)
    axes[0].text(x + x_scale*0.05, y, '$^{238}$U (mantle)')

    (x, y) = get_param_pos(planet_dynamics['Time'], planet_dynamics['235UNumMan'], 0.2)
    axes[0].text(x + x_scale*0.01, y, '$^{235}$U (mantle)')

    axes[0].set_ylabel('Number of Atoms', fontsize = 12)
    axes[0].grid(True, alpha = 0.25)
    axes[0].xaxis.set_major_formatter(tck.FuncFormatter(lambda x, pos: x*1e-9))
    axes[0].set_yscale('log')

    magmom = planet_dynamics['MagMom']

    min_mag_moment_index = np.argmin(magmom)
    time_at_min = planet_dynamics['Time'][min_mag_moment_index]

    axes[1].plot(planet_dynamics['Time'], magmom, alpha = 0.7, color = 'purple')
    axes[1].xaxis.set_major_formatter(tck.FuncFormatter(lambda x, pos: x*1e-9))
    axes[1].set_ylabel('Magnetic Moment\n[Earth Magnetic Moments]', fontsize=12)
    axes[1].grid(True, alpha = 0.25)

    axes[1].text(4.1e9, 1.4, 'Core\nsolidifying', fontsize = 10, color = 'black')

    (vmin, vmax) = axes[1].get_ylim()
    
    axes[1].vlines([time_at_min], vmin, vmax, linestyles = '--', color = 'black', alpha = 0.75)

    percent_mass_loss = 100 * (planet_dynamics['EnvelopeMass'][0] - planet_dynamics['EnvelopeMass'])/planet_dynamics['EnvelopeMass'][0]

    axes[2].plot(planet_dynamics['Time'], percent_mass_loss)
    axes[2].set_ylabel('% Atmospheric Mass Loss', fontsize = 12)
    axes[2].grid(True, alpha = 0.25)
    axes[2].xaxis.set_major_formatter(tck.FuncFormatter(lambda x, pos: x*1e-9))
    axes[2].set_xlabel('Time [Gyr]', fontsize = 12)

    axes[1].set_ylim(vmin, vmax)

    if save_name != '':
        fig.savefig(save_name, dpi = 300, bbox_inches = 'tight')