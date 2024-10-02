'''Johnson and Brandis Dragonfly Proposed Entry Trajectory into Titan '''
'''Johnston CO, West TK, Brandis AM. Features of Afterbody Radiative Heating for
Titan Entry. In: AIAA Aviation 2019 Forum; 2019. .'''

import numpy as np
from scipy.interpolate import CubicSpline

CpN2_300 = 1039  # J/kg-K
CpCH4_300 = 2220  # J/kg-K
Cp_mix = (0.98 * 1039) + (0.02 * 2220)

'''Find R for N2 98% and CH4 2%'''
''' Methane % at 247 km?'''
gamma = 1.396  # previous tests
R_N2 = 296.8  # J
R_CH4 = 518.3  # J
Rmix = (0.98 * R_N2) + (0.02 * R_CH4)
print(f'R mixture is {Rmix}')
MX2 = 10
Area_Dfly = (np.pi * Dflight ** 2) / 4

Titan_interface = 1270  # km
Ve = 7330  # m/s
q_y_max = 247  # km
v_y_max = 5770  # m/s
q_max = 291  # W/cm^2

y_g_max = 219  # km
v_g_max = 4550  # m/s

yo = 40000  # 40 km Scale height

Dflight = 4.5  # m
Rnose = 1.309  # m
D_X2 = 0.09  # m 90 mm
mass = 2255  # kg entry mass pre ablation


def pressure(rho, Rmix, T):
    pressure = rho * Rmix * T
    return pressure


def enthalpy(v_inf, Cp_mix, T_inf):
    enthalpy = (v_inf ** 2 / 2) + (Cp_mix * T_inf)
    return enthalpy


# density dictionary
rho_flight = 9.02e-4


def rho_lab(rho_inf, D_flight, D_x2):
    rho_x2 = (rho_inf * D_flight) / D_x2
    return rho_x2


# Mach 10 nozzle in X2
def T_X2(ho, MX2, gamma, Rmix, Cp_mix):
    T_X2 = (ho / ((MX2 ** 2 * gamma * Rmix / 2) + Cp_mix))
    return T_X2


def V_X2(MX2, gamma, Rmix, T_X2):
    V_X2 = MX2 * np.sqrt(gamma * Rmix * T_X2)
    return V_X2


def sound_speed(gamma, R, Temperature):
    return np.sqrt(gamma * R * Temperature)


def Mach(velocity, sound_speed):
    return velocity / sound_speed


def normal_shock_relations(M1, gamma):
    """Calculate downstream Mach number M2 in a normal shock."""
    M2 = np.sqrt((1 + 0.5 * (gamma - 1) * M1 ** 2) / (gamma * M1 ** 2 - 0.5 * (gamma - 1)))
    return M2


def shock_P2_P1(M1, gamma):
    """ Calculate the pressure ratio across the shock (p2/p1)."""
    pressure_ratio = 1 + 2 * gamma / (gamma + 1) * (M1 ** 2 - 1)
    return pressure_ratio


def shock_density_ratio(M1, gamma):
    """ Calculate the density ratio across the shock (ρ2/ρ1)."""
    density_ratio = ((gamma + 1) * M1 ** 2) / (2 + (gamma - 1) * M1 ** 2)
    return density_ratio


def T2_T1_ratio(M1, gamma):
    """     Calculate the temperature ratio across the shock (T2/T1). """
    left = 1 + ((2 * gamma) / (gamma + 1)) * (M1 ** 2 - 1)
    right = ((2 + (gamma - 1) * M1 ** 2) / ((gamma + 1) * M1 ** 2))
    T2_T1 = left * right
    return T2_T1


''' Dictionary of Table Values of Entry Trajectory'''
entry_dict = {
    "time": [160, 180, 190, 202, 211, 224, 235, 241, 245, 248, 254, 263, 280],  # List of time values
    "velocity": [7330, 7290, 7240, 7100, 6890, 6350, 5590, 5070, 4700, 4420, 3850, 3080, 1940],
    # List of velocity values
    "density": [4.4e-6, 1.94e-5, 3.98e-5, 8.74e-5, 1.58e-4, 3.56e-4, 6.64e-4, 9.02e-4, 1.09e-3, 1.24e-3, 1.57e-3,
                2.12e-3, 3.26e-3],  # List of density values
    "temperature": [162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 160, 158]  # List of temperature values
}

# New dictionary to store interpolated values
interpolated_dict = {
    "time": [],
    "velocity": [],
    "density": [],
    "temperature": []
}

# Define the new range of time for interpolation (you can adjust this as needed)
time_interp = np.linspace(min(entry_dict["time"]), max(entry_dict["time"]),
                          num=50)  # 50 points between the min and max time

# Perform cubic spline interpolation for each variable
cs_velocity = CubicSpline(entry_dict["time"], entry_dict["velocity"])
cs_density = CubicSpline(entry_dict["time"], entry_dict["density"])
cs_temperature = CubicSpline(entry_dict["time"], entry_dict["temperature"])

# Interpolate the values at the new time points
velocity_interp = cs_velocity(time_interp)
density_interp = cs_density(time_interp)
temperature_interp = cs_temperature(time_interp)

# Store the interpolated values in the new dictionary
interpolated_dict["time"] = time_interp
interpolated_dict["velocity"] = velocity_interp
interpolated_dict["density"] = density_interp
interpolated_dict["temperature"] = temperature_interp

# Example loop to perform calculations across all values
for i in range(len(interpolated_dict["time"])):
    time = interpolated_dict["time"][i]
    velocity = interpolated_dict["velocity"][i]
    density = interpolated_dict["density"][i]
    temperature = interpolated_dict["temperature"][i]

    titan_freestream_P = pressure(density, Rmix, temperature)
    print(f"When Titan velocity is {velocity:.0f} m/s, Titan Freestream Pressure is {titan_freestream_P:.2f} Pa.")
    enthalpy_titan = enthalpy(velocity, Cp_mix, temperature)
    print(f"When Titan velocity is {velocity:.0f} m/s, Titan Entry Enthalpy is {enthalpy_titan / 1e6:.2f} MJ/Kg-K.")
    rho_X2 = rho_lab(density, Dflight, D_X2)
    print(f"When Titan velocity is {velocity:.0f} m/s, X2 Freestream density is {density:.6f} kg/m^3.")
    Temp_X2 = T_X2(enthalpy_titan, MX2, gamma, Rmix, Cp_mix)
    print(f"When Titan velocity is {velocity:.0f} m/s, X2 Freestream Temperature is {Temp_X2:.2f} K.")
    Velocity_X2 = V_X2(MX2, gamma, Rmix, Temp_X2)
    print(f"When Titan velocity is {velocity:.0f} m/s, X2 Freestream Velocity is {Velocity_X2:.2f} m/s.")
    Pressure_X2 = pressure(rho_X2, Rmix, Temp_X2)
    print(f"When Titan velocity is {velocity:.0f} m/s, X2 Freestream Pressure is {Pressure_X2:.2f} Pa.")
    print()
    Titan_soundspeed = sound_speed(gamma, Rmix, temperature)
    print(f"When Titan velocity is {velocity:.0f} m/s, Titan local sound speed is {Titan_soundspeed:.2f} m/s.")
    Titan_Mach = Mach(velocity, Titan_soundspeed)
    print(f"When Titan velocity is {velocity:.0f} m/s, Titan local Mach number is {Titan_Mach:.2f}.")

    Mach_behind = normal_shock_relations(Titan_Mach, gamma)
    print(f"When Titan velocity is {velocity:.0f} m/s, Mach number behind Titan bow shock is {Mach_behind:.2f}.")

    Pressure_PostDfly_shock = shock_P2_P1(Titan_Mach, gamma)
    print(
        f"When Titan velocity is {velocity:.0f} m/s, Pressure behind Titan bow shock is {(Pressure_PostDfly_shock * titan_freestream_P) / 1000:.2f}.kPa ")

    Density_PostDfly_shock = shock_density_ratio(Titan_Mach, gamma)
    print(
        f"When Titan velocity is {velocity:.0f} m/s, Density behind Titan bow shock is {(Density_PostDfly_shock * density):.6f}.kg/m^3 ")

    Temperature_PostDfly_shock = T2_T1_ratio(Titan_Mach, gamma)
    print(
        f"When Titan velocity is {velocity:.0f} m/s, Temperature behind Titan bow shock is {(Temperature_PostDfly_shock * temperature):.2f} K.")

    print('---' * 20)

''' Ballistic Equations '''


def y_f_max(yo, Cd, rho, Area, mass, flight_angle):
    ymax_decel = yo * np.log((Cd * rho * Area * yo) / (mass * np.sin(np.radians(flight_angle))))
    return ymax_decel


alt_gmax = y_f_max(yo, 0.7, 5.3, Area_Dfly, mass, 90.0)

print(f" Altitude of g Max for Titan entry is {alt_gmax:.2f} m")

titan_surface_dens = 5.3  # 219 km gmax

V_gmax = 7.5 * (np.e ** -0.5)
print(f"velocity at g max is {V_gmax:.2f} m/s, they get 4.55 km/s")
#
# '''Normal Shock Relations'''
#
# def normal_shock_relations(M1, gamma):
#     """ Calculate downstream Mach number M2 in a normal shock."""
#     M2 = np.sqrt((1 + 0.5 * (gamma - 1) * M1**2) / (gamma * M1**2 - 0.5 * (gamma - 1)))
#     return M2
#
# def shock_P2_P1(M1, gamma):
#     """ Calculate the pressure ratio across the shock (p2/p1)."""
#     pressure_ratio = 1 + 2*gamma/(gamma + 1) * (M1**2 - 1)
#     return pressure_ratio
#
# u = shock_P2_P1(21.88, gamma)
# print('thing',u)
# post_pressure = u * 29.85
# print('post shock pressure',post_pressure)
#
# def shock_density_ratio(M1, gamma):
#     """ Calculate the density ratio across the shock (ρ2/ρ1)."""
#     density_ratio = ((gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2 )
#     return density_ratio
#
# def T2_T1_ratio(M1, gamma):
#     """     Calculate the temperature ratio across the shock (T2/T1). """
#     # Using previously defined functions to get p2/p1 and ρ2/ρ1
#     p_ratio = shock_P2_P1(M1, gamma)
#     rho_ratio = shock_density_ratio(M1, gamma)
#     temperature_ratio = p_ratio / rho_ratio
#     return temperature_ratio
#
#
#
#




