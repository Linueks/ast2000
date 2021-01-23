from __future__ import division
from ast2000solarsystem_27_v5 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import random as random
from matplotlib import style

np.random.seed(1)

style.use('ggplot')
star_system_seed = 11466
star_system = AST2000SolarSystem(star_system_seed)
gravitational_constant = 6.67408e-11                                                                            # [m**3 * kg**-1 * s**-2]
hydrogen_gass_molar_mass = 2.01588                                                                              # [gr/mol]
boltzmann_constant = 1.380648e-23                                                                               # [m**2 * kg * s**-2 * K**-1]
avogadro_constant = 6.02214179e23                                                                               # [mol**-1]
particle_mass = hydrogen_gass_molar_mass / (avogadro_constant * 1e3)                                            # [kg]
number_of_particles = 100000
temperature = 10000                                                                                             # [K]
initial_fuel_mass = 100000                                                                                      # [kg]
satellite_mass = 1100                                                                                           # [kg]
total_mass = initial_fuel_mass + satellite_mass                                                                 # [kg]
box_dim = 1e-6                                                                                                  # [m]

planets_radii = star_system.radius                                                                              # [km]
planets_mass = star_system.mass                                                                                 # [Solar Masses]
planets_orbital_period = star_system.period * 24 * 60 * 60                                                      # [s]

home_planet_mass= 1.98855e30 * planets_mass[0]                                                                  # [kg]
home_planet_radius = planets_radii[0] * 1000                                                                    # [m]
home_planet = [home_planet_radius, home_planet_mass]
home_planet_gravity = gravitational_constant * home_planet[1] / home_planet[0]**2
home_planet_initial_position = np.array([star_system.x0[0], star_system.y0[0]])
home_planet_initial_velocity = np.array([star_system.vx0[0], star_system.vy0[0]])
home_planet_rotational_speed = (2 * np.pi * planets_orbital_period[0]/ (24 * 60 * 60))


def simulate_engine(number_of_particles, timesteps, dt):
    standard_deviation = np.sqrt(boltzmann_constant * temperature / particle_mass)
    mean = 0.0
    dimensions = 3
    hole_size = box_dim / 2
    hole = (box_dim / 2 - hole_size / 2, box_dim / 2 + hole_size / 2)

    particle_positions = np.zeros((number_of_particles, dimensions, timesteps))
    particle_velocities = np.zeros((number_of_particles, dimensions))

    for particle in xrange(number_of_particles):
        for dimension in xrange(dimensions):
            particle_positions[particle, dimension, 0] = np.random.uniform(0, box_dim)
            particle_velocities[particle, dimension] = np.random.normal(mean, standard_deviation)

    number_of_escapes = 0.0
    partition_momentum = 0.0


    for i in xrange(timesteps-1):
        is_outside_box_1 = particle_positions[:, :, i] > box_dim
        is_outside_box_2 = particle_positions[:, :, i] < 0

        particle_velocities[is_outside_box_1] *= -1             #swaps velocity direction if particle is outside box
        particle_velocities[is_outside_box_2] *= -1

        hole_mask_x = particle_positions[:, 0, i] < 0
        hole_mask_y = np.logical_and(particle_positions[:, 1, i] > hole[0], particle_positions[:, 1, i] < hole[1])
        hole_mask_z = np.logical_and(particle_positions[:, 2, i] > hole[0], particle_positions[:, 2, i] < hole[1])

        hole_mask_yz = np.logical_and(hole_mask_y, hole_mask_z)
        hole_mask_xyz = np.logical_and(hole_mask_x, hole_mask_yz)

        hole_hit = np.sum(hole_mask_xyz)
        number_of_escapes += hole_hit

        particle_positions[:, :, i+1] = particle_positions[:, :, i] + particle_velocities[:, :] * dt
        partition_momentum += 2 * particle_mass * sum(particle_velocities[hole_mask_xyz, 0])                    # the 2 right here is to check if fuel use is better part1a


    partition_force = partition_momentum / (timesteps * dt)
    particles_per_second = number_of_escapes / (timesteps * dt)


    return partition_force, particles_per_second


def simulate_launch(engine_partitions, partition_force, partition_mass_difference):
    launch_duration = 20. * 60
    time_steps = 100000
    dt = launch_duration / time_steps
    time = 0

    engine_force = partition_force * engine_partitions
    engine_fuel_consumption = partition_mass_difference * engine_partitions

    rocket_velocity = 0
    radial_position = home_planet_radius
    current_mass = total_mass

    home_planet_gravity_force = lambda r: gravitational_constant * home_planet_mass * current_mass / (r)**2
    home_planet_escape_velocity = lambda r: np.sqrt(2 * gravitational_constant * home_planet_mass / (r))
    force_array = np.zeros(time_steps)
    mass_array = np.zeros(time_steps)



    for i in xrange(time_steps):
        force_on_rocket = engine_force - home_planet_gravity_force(radial_position)
        force_array[i] = force_on_rocket
        rocket_velocity += force_on_rocket / current_mass * dt
        radial_position += rocket_velocity * dt
        current_mass -= engine_fuel_consumption * dt
        mass_array[i] = current_mass
        time += dt

        if current_mass <= satellite_mass:
            print 'Out of fuel'
        elif rocket_velocity < 0:
            print 'Gravity sucks'
        elif rocket_velocity >= home_planet_escape_velocity(radial_position):
            print 'Rocket escaped, will not get pulled back ever'
            fuel_spent = initial_fuel_mass - (current_mass - satellite_mass)
            return fuel_spent, time, radial_position, rocket_velocity, force_array, mass_array, i

    return 'ran through loop without entering ifs'


def plot(number_of_particles, timesteps, dt):
    def update_lines(num, data_lines, lines):
        for line, data in zip(lines, data_lines):
            line.set_data(data[0 : 2, num-1 : num])
            line.set_3d_properties(data[2, num-1 : num])
        return lines

    fig = plt.figure()
    ax = p3.Axes3D(fig)
    frames = 100


    total_momentum, data_x = simulate_engine(number_of_particles, timesteps, dt)
    lines = []
    data = []
    for i in range(number_of_particles):
        data.append([data_x[i]])
        lines.append([ax.plot(data[i][0][0,0:1], data[i][0][1,0:1], data[i][0][2,0:1], 'o')[0]])

    ax.set_title('Particle Animation')
    ax.set_xlim3d([0.0, box_dim])
    ax.set_xlabel('x')

    ax.set_ylim3d([0.0, box_dim])
    ax.set_ylabel('y')

    ax.set_zlim3d([0.0, box_dim])
    ax.set_zlabel('z')

    ani = list(range(number_of_particles))
    for i in ani:
        ani[i] = animation.FuncAnimation(fig, update_lines, frames, fargs=(data[i], lines[i]),
                                interval=50, blit=False)

    engine_hole = patches.Circle((box_dim/2, box_dim/2), box_dim/4)
    ax.add_patch(engine_hole)
    art3d.pathpatch_2d_to_3d(engine_hole, z=0, zdir='z')
    engine_hole.set_facecolor('black')

    plt.show()


def calculate_sat_position_relative_to_sun(launch_duration, launch_end_velocity, distance_from_planet_center):
    total_time = launch_duration * 3.17098e-8                                                                   # [s]
    time_steps = 10000
    dt = total_time / time_steps

    planet_position = home_planet_initial_position * 149597870700                                               # [m]
    planet_velocity = home_planet_initial_velocity * 4743.7173611                                               # [m/s]
    launch_angle = np.pi/2 - np.arctan2(home_planet_rotational_speed, launch_end_velocity)

    initial_satellite_position_relative_to_planet = [home_planet_radius * np.cos(launch_angle), home_planet_radius * np.sin(launch_angle)]
    satellite_position_relative_to_planet = [distance_from_planet_center * np.cos(launch_angle), distance_from_planet_center * np.sin(launch_angle)]
    satellite_velocity_relative_to_planet = [0, np.sqrt(launch_end_velocity**2 + home_planet_rotational_speed**2)]

    distance_from_rotation = home_planet_rotational_speed * launch_duration
    satellite_position_relative_to_planet += np.array([distance_from_rotation * np.cos(np.pi/2 - launch_angle, distance_from_rotation * np.sin(np.pi/2 - launch_angle))])


    for t in xrange(time_steps):
        planet_position += planet_velocity * dt


    satellite_position = (planet_position + satellite_position_relative_to_planet) / 149597870700

    return satellite_position, initial_satellite_position_relative_to_planet


analytical_acceleration_per_box = number_of_particles * boltzmann_constant * temperature / (2 * box_dim * hydrogen_gass_molar_mass)

force_to_beat_gravity = home_planet_gravity * total_mass + 1e-6

#partition_force, particles_per_sec_per_partition = print simulate_engine(number_of_particles, 1000, 1e-12)

# test parameters
partition_force = 3.3985997328656428e-09
partition_mass_difference = 63135999999999.99 * particle_mass

number_of_partitions = force_to_beat_gravity / partition_force

fuel_spent, time, radial_position, rocket_velocity, force_array, mass_array, i = simulate_launch(number_of_partitions, partition_force, partition_mass_difference)


def make_accel_graph(force, mass, i):
    accel = force / mass

    plt.xlabel('mass [kg]')
    plt.ylabel('acceleration [m/s**2]')
    times = np.linspace(0, i, 100000)

    plt.subplot(111)
    plt.plot(times, accel)
    plt.subplot(211)
    plt.plot(mass_array, accel)
    plt.show()


make_accel_graph(force_array, mass_array, i)
