import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

# Units in number of atoms.
EARTH_233TH_CORE = 1.35657e+40
EARTH_232TH_MANTLE = 6.07873e+41
EARTH_238U_CORE = 2.95900e+39
EARTH_238U_MANTLE = 2.86829e+41
EARTH_235U_CORE = 2.45149e+39
EARTH_235U_MANTLE = 9.89663e+40
EARTH_40K_MANTLE = 9.72110e+42
EARTH_40K_CORE = 9.10261e+42

class Planet:
    G = 6.674 * 1e-11 # SI base units.
    EPSILON_0 = 8.82 * 1e-12 # C^2 kg^-1 m^-3 s^2 (SI base units)

    def __init__(self,
                 stellar_mean_surface_B,
                 stellar_radius,
                 planet_star_distance,
                 planet_dynamics,
                 stellar_dynamics,
                 planet_radius = 6.371 * 1e6, # In meters (SI base units). 
                 planet_mass = 5.97 * 1e24, # In kg (SI base units).
                 eta = 0.1):
        self.planet_dynamics = planet_dynamics
        self.stellar_dynamics = stellar_dynamics

        SOLAR_LUMINOSITY = 3.828 * 1e26 # Watts (SI base units).
        XUV_flux = (stellar_dynamics['LXUVStellar'] * SOLAR_LUMINOSITY) / (4*np.pi*(planet_star_distance**2)) # SI base units.
        EARTH_MAGNETIC_MOMENT = 7.94 * 1e22 # A m^2 (SI base units). From https://modern-physics.org/earths-dipole-moment/
        self.planet_dipole_strength = planet_dynamics['MagMom'] * EARTH_MAGNETIC_MOMENT # SI base units.

        self.LIMITING_INDEX = len(XUV_flux)

        self.eta = eta
        self.planet_radius = planet_radius
        self.planet_mass = planet_mass
        self.planet_mean_surface_B = (4*np.pi*self.EPSILON_0*3*self.planet_dipole_strength[:self.LIMITING_INDEX]) / (2.0*planet_radius**3) # Tesla (SI base units).
        self.stellar_mean_surface_B = stellar_mean_surface_B
        self.stellar_radius = stellar_radius
        self.planet_star_distance = planet_star_distance
        self.time = stellar_dynamics['Time'][:self.LIMITING_INDEX]
        self.XUV_flux = XUV_flux
        self.eta = eta

    def dm_dt_energy_limited(self):
        return (self.eta*np.pi*(self.planet_radius**3)*self.XUV_flux) / (self.G*self.planet_mass)

    def Beta(self):
        # Prevents division by zero error outputs that are annoying and unnecessary.
        with np.errstate(divide = 'ignore'):
            return (self.stellar_mean_surface_B/self.planet_mean_surface_B) * np.power(self.stellar_radius/self.planet_star_distance, 3) # Unitless.

    def F_AP(self):
        Beta = self.Beta()

        F_AP = 1 - np.power(1 - (3*np.power(Beta, 1/3))/(2+Beta), 1/2) # Unitless.
        
        F_AP[Beta >= 1] = 1

        return F_AP
    
    def dm_dt_magnetic(self):
        return self.dm_dt_energy_limited() * self.F_AP()
    
    def integrate(self, magnetic = True):
        dm_dt = self.dm_dt_magnetic() if magnetic else self.dm_dt_energy_limited()

        max_time_index = len(dm_dt)
        mass_loss = np.zeros(max_time_index)

        for time_index in range(1, max_time_index):
            dt = (self.time[time_index] - self.time[time_index-1]) * 3.1536e7/1 # Convert to seconds.
            dm = dm_dt[time_index] * dt

            mass_loss[time_index] = mass_loss[time_index-1] + dm
        
        return mass_loss

    def output(self):
        df = pd.DataFrame()
        df['F_AP'] = self.F_AP()
        df['beta'] = self.Beta()
        df['planetB'] = self.planet_mean_surface_B
        
        semimajorAxis = np.array(self.LIMITING_INDEX)
        semimajorAxis.fill(self.planet_star_distance)

        stellarRadius = np.array(self.LIMITING_INDEX)
        stellarRadius.fill(self.stellar_radius)

        df['semimajorAxis'] = semimajorAxis
        df['stellarRadius'] = stellarRadius

        df.to_csv('planet_output.csv', index=False)
    
    def plot_thermal_evolution(self):
        (fig, axes) = plt.subplots(1, 3)

        fig.set_size_inches(14.5, 3.75)
        fig.subplots_adjust(wspace=0.225)

        MAX_TIME = self.time.iloc[-1]

        axes[0].plot(self.time, self.planet_mean_surface_B[:len(self.time)], alpha = 0.7)
        axes[0].set_xlim(0.0, MAX_TIME)
        axes[0].xaxis.set_major_formatter(tck.FuncFormatter(lambda x,pos: x*1e-9))
        axes[0].set_ylabel('Mean surface B-field [Tesla]', fontsize=12)
        axes[0].set_xlabel('Time [Gyr]')
        axes[0].grid(True, alpha = 0.25)

        #axes[0].legend(loc = 'best', fontsize = 10)

        axes[1].plot(self.time, self.planet_dynamics['40KNumCore'][:len(self.time)], label = '$^{40}$K (core)')
        axes[1].plot(self.time, self.planet_dynamics['232ThNumMan'][:len(self.time)], label = '$^{232}$Th (mantle)')
        axes[1].plot(self.time, self.planet_dynamics['235UNumMan'][:len(self.time)], label = '$^{235}$U (mantle)')
        axes[1].plot(self.time, self.planet_dynamics['238UNumMan'][:len(self.time)], label = '$^{238}$U (mantle)')
        axes[1].legend(ncol = 1, fontsize = 10)

        axes[1].set_ylabel('Number of atoms', fontsize = 12)
        axes[1].set_xlim(0.0, MAX_TIME)
        axes[1].grid(True, alpha = 0.25)
        axes[1].xaxis.set_major_formatter(tck.FuncFormatter(lambda x, pos: x*1e-9))
        axes[1].set_xlabel('Time [Gyr]')
        axes[1].set_yscale('log')

        axes[2].plot(self.time, self.planet_dynamics['TCore'][:len(self.time)], label = 'Core temp')
        axes[2].plot(self.time, self.planet_dynamics['TLMan'][:len(self.time)], label = 'Mantle temp')
        axes[2].plot(self.time, self.planet_dynamics['TCMB'][:len(self.time)], label = 'CMB temp')

        axes[2].legend(loc='best', fontsize = 10)
        axes[2].set_xscale('log')
        axes[2].set_ylabel('Temperature [K]', fontsize = 12)
        axes[2].set_xlim(0.0, MAX_TIME)
        axes[2].grid(True, alpha = 0.25)
        axes[2].xaxis.set_major_formatter(tck.FuncFormatter(lambda x, pos: x*1e-9))
        axes[2].set_xlabel('Time [Gyr]')

        initial_K = self.planet_dynamics['40KNumCore'].iloc[0] / EARTH_40K_CORE
        initial_ThU = self.planet_dynamics['232ThNumMan'].iloc[0] / EARTH_232TH_MANTLE

        fig.suptitle('Thermal Intererior Evolution of Earth-Like Planet (Initial K = {_40K} [Primordial Earth]; Initial Th, U = {_ThU} [Primordial Earth])'.format(_40K=round(initial_K, 2), _ThU=round(initial_ThU, 2)), fontsize = 15)
