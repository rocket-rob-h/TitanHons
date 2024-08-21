# import numpy as np
#
# # Constants
# PI = np.pi
#
# # Function to determine the area of the disc
# def disc_area(radius):
#     """
#     Determines the area of the disc.
#
#     Parameters:
#     radius (float): The radius of the disc (m).
#
#     Returns:
#     float: The area of the disc (m^2).
#     """
#     return PI * radius**2
#
# # Function to determine the force at the stagnation point
# def stagnation_force(stagnation_pressure, area):
#     """
#     Determines the force at the stagnation point on the disc.
#
#     Parameters:
#     stagnation_pressure (float): The stagnation pressure acting on the disc (Pa).
#     area (float): The area of the disc (m^2).
#
#     Returns:
#     float: The force at the stagnation point (N).
#     """
#     return stagnation_pressure * area
#
# # Function to determine the bending moment
# def bending_moment(force, beam_length, angle_rad):
#     """
#     Determines the bending moment at the base of the beam.
#
#     Parameters:
#     force (float): The force applied to the disc (N).
#     beam_length (float): The length of the beam (m).
#     angle_rad (float): The angle of the beam with respect to the horizontal, in radians.
#
#     Returns:
#     float: The bending moment at the base of the beam (Nm).
#     """
#     moment_arm = beam_length * np.cos(angle_rad)
#     return force * moment_arm
#
# # Function to determine the section modulus for a rectangular beam
# def section_modulus_rectangular(beam_width, beam_thickness):
#     """
#     Determines the section modulus of a rectangular beam.
#
#     Parameters:
#     beam_width (float): The width of the beam (m).
#     beam_thickness (float): The thickness of the beam (m).
#
#     Returns:
#     float: The section modulus of the rectangular beam (m^3).
#     """
#     return (beam_width * beam_thickness**2) / 6
#
# # Function to determine the second moment of area for a rectangular beam
# def second_moment_area_rectangular(beam_width, beam_thickness):
#     """
#     Determines the second moment of area (area moment of inertia) for a solid rectangular beam.
#
#     Parameters:
#     beam_width (float): The width of the beam (m).
#     beam_thickness (float): The thickness of the beam (m).
#
#     Returns:
#     float: The second moment of area of the rectangular beam (m^4).
#     """
#     return (beam_width * beam_thickness**3) / 12
#
# # Function to determine the second moment of area for a hollow rectangular beam
# def second_moment_area_hollow_rectangular(external_width, external_height, thickness):
#     """
#     Determines the second moment of area (area moment of inertia) for a hollow rectangular beam.
#
#     Parameters:
#     external_width (float): The external width of the beam (m).
#     external_height (float): The external height of the beam (m).
#     thickness (float): The thickness of the walls of the beam (m).
#
#     Returns:
#     float: The second moment of area of the hollow rectangular beam (m^4).
#     """
#     internal_width = external_width - 2 * thickness
#     internal_height = external_height - 2 * thickness
#     return (external_width * external_height**3 - internal_width * internal_height**3) / 12
#
# # Function to determine the second moment of area for a circular rod
# def second_moment_area_circular(diameter):
#     """
#     Determines the second moment of area (area moment of inertia) for a solid circular rod.
#
#     Parameters:
#     diameter (float): The diameter of the circular rod (m).
#
#     Returns:
#     float: The second moment of area of the circular rod (m^4).
#     """
#     return (PI * diameter**4) / 64
#
# # Function to determine the second moment of area for a hollow circular rod
# def second_moment_area_hollow_circular(external_diameter, internal_diameter):
#     """
#     Determines the second moment of area (area moment of inertia) for a hollow circular rod.
#
#     Parameters:
#     external_diameter (float): The external diameter of the hollow circular rod (m).
#     internal_diameter (float): The internal diameter of the hollow circular rod (m).
#
#     Returns:
#     float: The second moment of area of the hollow circular rod (m^4).
#     """
#     return (PI / 64) * (external_diameter**4 - internal_diameter**4)
#
# # Function to determine the stress in the beam using section modulus
# def stress_from_section_modulus(moment, section_modulus):
#     """
#     Determines the maximum stress in the beam using the section modulus.
#
#     Parameters:
#     moment (float): The bending moment applied to the beam (Nm).
#     section_modulus (float): The section modulus of the beam (m^3).
#
#     Returns:
#     float: The maximum stress in the beam (Pa).
#     """
#     return moment / section_modulus
#
# # Function to determine the stress in the rectangular beam using second moment of area
# def stress_from_second_moment_rectangular(moment, second_moment, beam_thickness):
#     """
#     Determines the maximum stress in a rectangular beam using the second moment of area.
#
#     Parameters:
#     moment (float): The bending moment applied to the beam (Nm).
#     second_moment (float): The second moment of area of the beam (m^4).
#     beam_thickness (float): The thickness of the beam (m).
#
#     Returns:
#     float: The maximum stress in the beam (Pa).
#     """
#     return moment * (beam_thickness / 2) / second_moment
#
# # Function to determine the stress in the circular rod using the second moment of area
# def stress_from_second_moment_circular(moment, second_moment, diameter):
#     """
#     Determines the maximum stress in a circular rod using the second moment of area.
#
#     Parameters:
#     moment (float): The bending moment applied to the rod (Nm).
#     second_moment (float): The second moment of area of the rod (m^4).
#     diameter (float): The diameter of the rod (m).
#
#     Returns:
#     float: The maximum stress in the rod (Pa).
#     """
#     return moment * (diameter / 2) / second_moment
#
# # Function to determine the stress in the hollow circular rod using the second moment of area
# def stress_from_second_moment_hollow_circular(moment, second_moment, external_diameter):
#     """
#     Determines the maximum stress in a hollow circular rod using the second moment of area.
#
#     Parameters:
#     moment (float): The bending moment applied to the rod (Nm).
#     second_moment (float): The second moment of area of the hollow rod (m^4).
#     external_diameter (float): The external diameter of the hollow rod (m).
#
#     Returns:
#     float: The maximum stress in the hollow rod (Pa).
#     """
#     return moment * (external_diameter / 2) / second_moment
#
#
# # Initial conditions
# radius = 0.08  # Radius of the disc (m)
# thickness = 0.60  # Thickness of the disc (m)
# beam_length = 0.2  # Length of the beam (m)
# external_diameter = 0.05  # External diameter of the hollow circular rod (m)
# internal_diameter = 0.03  # Internal diameter of the hollow circular rod (m)
# beam_angle_deg = 45  # Angle of inclination of the beam (degrees)
# youngs_modulus = 210e9  # Young's modulus of the material (Pa)
# stagnation_pressure = 40e6  # Stagnation pressure at the disc (Pa)
#
# # Calculations
# area = disc_area(radius)
# force = stagnation_force(stagnation_pressure, area)
# angle_rad = np.radians(beam_angle_deg)
# moment = bending_moment(force, beam_length, angle_rad)
# second_moment = second_moment_area_hollow_circular(external_diameter, internal_diameter)
# stress = stress_from_second_moment_hollow_circular(moment, second_moment, external_diameter)
#
# # Output results
# print(f"Stagnation Force: {force:.2f} N")
# print(f"Bending Moment at the base of the beam: {moment:.2f} Nm")
# print(f"Maximum Stress in the hollow circular rod: {stress:.2f} Pa")


import numpy as np

# Constants
PI = np.pi

# Function to determine the area of the disc
def disc_area(radius):
    """
    Determines the area of the disc.

    Parameters:
    radius (float): The radius of the disc (m).

    Returns:
    float: The area of the disc (m^2).
    """
    return PI * radius**2

# Function to determine the force at the stagnation point
def stagnation_force(stagnation_pressure, area):
    """
    Determines the force at the stagnation point on the disc.

    Parameters:
    stagnation_pressure (float): The stagnation pressure acting on the disc (Pa).
    area (float): The area of the disc (m^2).

    Returns:
    float: The force at the stagnation point (N).
    """
    return stagnation_pressure * area

# Function to determine the bending moment
def bending_moment(force, beam_length, angle_rad):
    """
    Determines the bending moment at the base of the beam.

    Parameters:
    force (float): The force applied to the disc (N).
    beam_length (float): The length of the beam (m).
    angle_rad (float): The angle of the beam with respect to the horizontal, in radians.

    Returns:
    float: The bending moment at the base of the beam (Nm).
    """
    moment_arm = beam_length * np.cos(angle_rad)
    return force * moment_arm

# Function to determine the section modulus for a rectangular beam
def section_modulus_rectangular(beam_width, beam_thickness):
    """
    Determines the section modulus of a rectangular beam.

    Parameters:
    beam_width (float): The width of the beam (m).
    beam_thickness (float): The thickness of the beam (m).

    Returns:
    float: The section modulus of the rectangular beam (m^3).
    """
    return (beam_width * beam_thickness**2) / 6

# Function to determine the second moment of area for a rectangular beam
def second_moment_area_rectangular(beam_width, beam_thickness):
    """
    Determines the second moment of area (area moment of inertia) for a solid rectangular beam.

    Parameters:
    beam_width (float): The width of the beam (m).
    beam_thickness (float): The thickness of the beam (m).

    Returns:
    float: The second moment of area of the rectangular beam (m^4).
    """
    return (beam_width * beam_thickness**3) / 12

# Function to determine the second moment of area for a hollow rectangular beam
def second_moment_area_hollow_rectangular(external_width, external_height, thickness):
    """
    Determines the second moment of area (area moment of inertia) for a hollow rectangular beam.

    Parameters:
    external_width (float): The external width of the beam (m).
    external_height (float): The external height of the beam (m).
    thickness (float): The thickness of the walls of the beam (m).

    Returns:
    float: The second moment of area of the hollow rectangular beam (m^4).
    """
    internal_width = external_width - 2 * thickness
    internal_height = external_height - 2 * thickness
    return (external_width * external_height**3 - internal_width * internal_height**3) / 12

# Function to determine the second moment of area for a circular rod
def second_moment_area_circular(diameter):
    """
    Determines the second moment of area (area moment of inertia) for a solid circular rod.

    Parameters:
    diameter (float): The diameter of the circular rod (m).

    Returns:
    float: The second moment of area of the circular rod (m^4).
    """
    return (PI * diameter**4) / 64

# Function to determine the second moment of area for a hollow circular rod
def second_moment_area_hollow_circular(external_diameter, internal_diameter):
    """
    Determines the second moment of area (area moment of inertia) for a hollow circular rod.

    Parameters:
    external_diameter (float): The external diameter of the hollow circular rod (m).
    internal_diameter (float): The internal diameter of the hollow circular rod (m).

    Returns:
    float: The second moment of area of the hollow circular rod (m^4).
    """
    return (PI / 64) * (external_diameter**4 - internal_diameter**4)

# Function to determine the stress in the beam using section modulus
def stress_from_section_modulus(moment, section_modulus):
    """
    Determines the maximum stress in the beam using the section modulus.

    Parameters:
    moment (float): The bending moment applied to the beam (Nm).
    section_modulus (float): The section modulus of the beam (m^3).

    Returns:
    float: The maximum stress in the beam (Pa).
    """
    return moment / section_modulus

# Function to determine the stress in the rectangular beam using second moment of area
def stress_from_second_moment_rectangular(moment, second_moment, beam_thickness):
    """
    Determines the maximum stress in a rectangular beam using the second moment of area.

    Parameters:
    moment (float): The bending moment applied to the beam (Nm).
    second_moment (float): The second moment of area of the beam (m^4).
    beam_thickness (float): The thickness of the beam (m).

    Returns:
    float: The maximum stress in the beam (Pa).
    """
    return moment * (beam_thickness / 2) / second_moment

# Function to determine the stress in the circular rod using the second moment of area
def stress_from_second_moment_circular(moment, second_moment, diameter):
    """
    Determines the maximum stress in a circular rod using the second moment of area.

    Parameters:
    moment (float): The bending moment applied to the rod (Nm).
    second_moment (float): The second moment of area of the rod (m^4).
    diameter (float): The diameter of the rod (m).

    Returns:
    float: The maximum stress in the rod (Pa).
    """
    return moment * (diameter / 2) / second_moment

# Function to determine the stress in the hollow circular rod using the second moment of area
def stress_from_second_moment_hollow_circular(moment, second_moment, external_diameter):
    """
    Determines the maximum stress in a hollow circular rod using the second moment of area.

    Parameters:
    moment (float): The bending moment applied to the rod (Nm).
    second_moment (float): The second moment of area of the hollow rod (m^4).
    external_diameter (float): The external diameter of the hollow rod (m).

    Returns:
    float: The maximum stress in the hollow rod (Pa).
    """
    return moment * (external_diameter / 2) / second_moment

# Function to determine the deflection of the beam
def deflection(beam_length, force, youngs_modulus, second_moment):
    """
    Determines the deflection of a cantilever beam with a point load at the end.

    Parameters:
    beam_length (float): The length of the beam (m).
    force (float): The force applied to the beam (N).
    youngs_modulus (float): The Young's modulus of the material (Pa).
    second_moment (float): The second moment of area of the beam (m^4).

    Returns:
    float: The deflection of the beam (m).
    """
    return (force * beam_length**3) / (3 * youngs_modulus * second_moment)

# Function to determine the factor of safety
def factor_of_safety(yield_strength, stress):
    """
    Determines the factor of safety based on the yield strength and maximum stress.

    Parameters:
    yield_strength (float): The yield strength of the material (Pa).
    stress (float): The maximum stress in the beam (Pa).

    Returns:
    float: The factor of safety.
    """
    return yield_strength / stress

# Initial conditions
radius = 0.08  # Radius of the disc (m)
thickness = 0.60  # Thickness of the disc (m)
beam_length = 0.2  # Length of the beam (m)
external_diameter = 0.06  # External diameter of the hollow circular rod (m)
internal_diameter = 0.04  # Internal diameter of the hollow circular rod (m)
beam_angle_deg = 45  # Angle of inclination of the beam (degrees)
youngs_modulus = 210e9  # Young's modulus of the material (Pa)
yield_strength = 250e6  # Yield strength of the material (Pa)
stagnation_pressure = 4e6  # Stagnation pressure at the disc (Pa)

# Calculations
area = disc_area(radius)
force = stagnation_force(stagnation_pressure, area)
angle_rad = np.radians(beam_angle_deg)
moment = bending_moment(force, beam_length, angle_rad)
second_moment = second_moment_area_hollow_circular(external_diameter, internal_diameter)
stress = stress_from_second_moment_hollow_circular(moment, second_moment, external_diameter)
beam_deflection = deflection(beam_length, force, youngs_modulus, second_moment)
fos = factor_of_safety(yield_strength, stress)
yielding = "Yes" if stress > yield_strength else "No"
deflection_percentage = (beam_deflection / beam_length) * 100
allowable_deflection = beam_length / 360
exceeds_allowable_deflection = "Yes" if beam_deflection > allowable_deflection else "No"

# Calculate the new beam angle
deflection_angle_rad = beam_deflection / beam_length

# Calculate the new beam angle (in radians and degrees)
new_beam_angle_rad = angle_rad + deflection_angle_rad
new_beam_angle_deg = np.degrees(new_beam_angle_rad)



# Output results
print(f"Stagnation Force: {force:.2f} N")
print(f"Bending Moment at the base of the beam: {moment:.2f} Nm")
print(f"Maximum Stress in the hollow circular rod: {stress:.2f} Pa")
print(f"Beam Deflection: {beam_deflection:.6f} m")
print(f"Will the beam yield? {yielding}")
print(f"Factor of Safety: {fos:.2f}")
print(f"Deflection as a percentage of beam length: {deflection_percentage:.2f}%")
print(f"Exceeds allowable deflection limit ({allowable_deflection:.6f} m)? {exceeds_allowable_deflection}")
print(f"Deflection-induced angle: {np.degrees(deflection_angle_rad):.4f} degrees")
print(f"New beam angle after deflection: {new_beam_angle_deg:.2f} degrees")


# 20 x pitot pressure?