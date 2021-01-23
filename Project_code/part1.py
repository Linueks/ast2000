import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.patches as patches
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
from ast2000solarsystem_27_v6 import AST2000SolarSystem
import matplotlib.pyplot as plt

np.random.seed(1)
star_system_seed = 11466

H2_molar_mass = 2.01588            #[g/mol]
avogadros_const = 6.02214179e23
boltzmann_const = 1.380648e-23
particle_mass = H2_molar_mass/(avogadros_const * 1e3) #[kg]
gravitational_constant = 6.67408e-11  # [m**3 * kg**-1 * s**-2]

star_system = AST2000SolarSystem(star_system_seed)


planets_radii = star_system.radius  # [km]
planets_mass = star_system.mass  # [Solar Masses]
home_planet_mass = 1.98855e30 * planets_mass[0]  # [kg]
home_planet_radius = planets_radii[0] * 1000  # [m]
home_planet_initial_position = np.array([star_system.x0[0], star_system.y0[0]])
home_planet_initial_velocity = np.array([star_system.vx0[0], star_system.vy0[0]])
home_planet_rotational_speed = (2 * np.pi * star_system.radius[0]*1000 / (star_system.period[0] * 24 * 60 * 60))
home_planet_rotational_angular_speed = (2 * np.pi * star_system.period[0] /(24 * 60 * 60))

fuel_mass = 100000
satellite_mass = 1100   # [kg]
total_mass = satellite_mass + fuel_mass


ms_to_AUyr = 0.000210945021
m_to_AU = 6.68458712226844e-12


def simulate_engine(number_of_particles, number_of_time_steps, dt, box_dimension, temperature):
    standard_deviation = np.sqrt(boltzmann_const * temperature / particle_mass)
    mean = 0
    number_of_escaping_particles = 0
    partition_momentum = 0.0
    hole_size = box_dimension / 2
    hole = (box_dimension/2 - hole_size/2, box_dimension/2 + hole_size/2)

    print "nkT = ", number_of_particles*boltzmann_const*temperature

    ##Initial Values
    positions = np.zeros((3, number_of_particles, number_of_time_steps))
    velocities = np.zeros((3, number_of_particles))

    for particle in range(number_of_particles):
        for dimension in range(3):
            positions[dimension, particle, 0] = np.random.uniform(0, box_dimension)
            velocities[dimension, particle] = np.random.normal(mean, standard_deviation)

    for t in range(number_of_time_steps-1):
        positions[:, :, t+1] = positions[:, :, t] + velocities[:, :] * dt

        ##Masking for collision with walls
        outside_box_mask1 = positions[:, :, t + 1] < 0
        outside_box_mask2 = positions[:, :, t + 1] > box_dimension
        velocities[outside_box_mask1] *= -1
        velocities[outside_box_mask2] *= -1

        ##Masking for particles leaving hole
        x_mask = positions[0, :, t + 1] < 0
        y_mask = np.logical_and(positions[1, :, t + 1] > hole[0], positions[1, :, t + 1] < hole[1])
        z_mask = np.logical_and(positions[2, :, t + 1] > hole[0], positions[2, :, t+ 1] < hole[1])
        xy_mask = np.logical_and(x_mask, y_mask)
        xyz_mask = np.logical_and(xy_mask, z_mask)

        particles_leaving_box = np.sum(xyz_mask)
        number_of_escaping_particles += particles_leaving_box

        partition_momentum += 2 * particle_mass * sum(velocities[0, xyz_mask]) # *particles_leaving_box...

    return partition_momentum/(number_of_time_steps*dt), number_of_escaping_particles/(number_of_time_steps*dt)
    #return partition_momentum, positions

def plot(number_of_particles, time_steps, dt, box_dim, T):
    def update_lines(num, data_lines, lines):
        for line, data in zip(lines, data_lines):
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2,num-1:num])
        return lines

    fig = plt.figure()
    ax = p3.Axes3D(fig)
    frames = 1000

    partition_momentum, position_data = simulate_engine(number_of_particles, time_steps, dt, box_dim, T)
    lines = []
    data = []

    for i in xrange(number_of_particles):
        data.append([position_data[:, i]])
        lines.append([ax.plot(data[i][0][0,0:1], data[i][0][1,0:1], data[i][0][2,0:1], 'o')[0]])

    ax.set_title('Particle Animation')
    ax.set_xlim3d([0.0, box_dim])
    ax.set_xlabel('Z')

    ax.set_ylim3d([0.0, box_dim])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, box_dim])
    ax.set_zlabel('X')

    hole = patches.Circle((box_dim/2, box_dim/2), box_dim/4)
    hole.set_facecolor('black')
    ax.add_patch(hole)
    art3d.pathpatch_2d_to_3d(hole, z=0, zdir="z")

    ani = [i for i in xrange(number_of_particles)]
    for i in ani:
        ani[i] = animation.FuncAnimation(fig, update_lines, frames, fargs=(data[i], lines[i]),
                                interval=50, blit=False)
    plt.savefig("engine_sim.png")
    plt.show()


def simulate_launch(number_of_partitions, partition_force, number_of_particles_per_partition_per_second):
    total_time = 20.0*60
    number_of_time_steps = 10000
    dt = total_time/number_of_time_steps

    engine_force = partition_force*number_of_partitions
    fuel_loss_per_second = number_of_particles_per_partition_per_second * particle_mass * number_of_partitions
    print engine_force, fuel_loss_per_second
    velocity = 0.
    time = 0.
    radial_position = home_planet_radius

    current_mass = total_mass

    plot_pos = np.zeros(number_of_time_steps)

    for t in range(number_of_time_steps):
        plot_pos[t] = radial_position
        velocity += (engine_force/current_mass - gravitational_constant*home_planet_mass/radial_position**2)*dt
        radial_position += velocity*dt
        escape_velocity = np.sqrt(2 * gravitational_constant * home_planet_mass / home_planet_radius) #radial_position)
        current_mass -= fuel_loss_per_second*dt
        time += dt

        if current_mass <= satellite_mass:
            print ("Out of fuel")
            return
        elif velocity < 0:
            print ("Gravity stronger than engine")
            return
        elif velocity >= escape_velocity:
            fuel_used = fuel_mass - (current_mass-satellite_mass)
            print "Reached escape velocity"
            x = np.linspace(0, time, t)
            plt.plot(x, plot_pos[:t])#[:t]
            plt.legend(["Satellite height"])
            plt.xlabel("Time[s]")
            plt.ylabel("Height above planet center[m]")
            plt.savefig("launch_sim.png")
            plt.show()
            return fuel_used, time, radial_position, velocity
    return 0., 0.

def calculate_position_relative_to_sun(launch_time, distance_from_planet_center, escape_velocity):
    initial_satellite_position_relative_to_planet = np.array([0, home_planet_radius])        #[m]
    initial_satellite_position_relative_to_star = home_planet_initial_position*149597870700 + initial_satellite_position_relative_to_planet                 #[m]

    planet_velocity = home_planet_initial_velocity * 4743.7173611
    radial_unit_vector = initial_satellite_position_relative_to_planet / np.linalg.norm(initial_satellite_position_relative_to_planet)
    satellite_position_relative_to_planet = distance_from_planet_center * radial_unit_vector    #[m]
    """
    angular_unit_vector = np.array([-radial_unit_vector[1], radial_unit_vector[0]])
    satellite_velocity_due_to_rotation = home_planet_rotational_speed*angular_unit_vector
    satellite_position_relative_to_planet += satellite_velocity_due_to_rotation * launch_time
    #satellite_position_relative_to_planet += planet_velocity*launch_time        #[m]
    satellite_position_relative_to_star = home_planet_initial_position*149597870700 + satellite_position_relative_to_planet
    satellite_position_relative_to_star += planet_velocity*launch_time
    """

    sat_pos_copy = satellite_position_relative_to_planet.copy()
    omega = 2.0 * np.pi * star_system.period[0] / 86400.0
    theta = omega * launch_time
    satellite_position_relative_to_planet[0] = sat_pos_copy[0] * np.cos(theta) - sat_pos_copy[1] * np.sin(theta)
    satellite_position_relative_to_planet[1] = sat_pos_copy[0] * np.sin(theta) + sat_pos_copy[1] * np.cos(theta)

    initial_rotation = np.array((-satellite_position_relative_to_planet[1], satellite_position_relative_to_planet[0]))
    rotation_unit = initial_rotation / np.linalg.norm(initial_rotation)
    rotation_vel = rotation_unit * omega * distance_from_planet_center
    end_velocity = escape_velocity * satellite_position_relative_to_planet/np.linalg.norm(satellite_position_relative_to_planet) + rotation_vel + planet_velocity


    satellite_position_relative_to_planet += planet_velocity * launch_time
    satellite_position_relative_to_star = satellite_position_relative_to_planet*m_to_AU + home_planet_initial_position


    initial_satellite_position_relative_to_star *= m_to_AU
    end_velocity *= ms_to_AUyr


    return initial_satellite_position_relative_to_star, satellite_position_relative_to_star, end_velocity, satellite_position_relative_to_planet*m_to_AU

def fuel_needed_for_deltav(deltav):
    engine_force = 1073357.35327                #[N]
    fuel_loss_per_second = 66.7474270958        #[kg/s]

    exhaust_velocity = engine_force/fuel_loss_per_second

    fuel_used = np.exp(deltav - exhaust_velocity)

    return fuel_used

partition_force = 3.4021029823078977e-09
particles_per_sec_per_partition = 63200999999999.992

gravity_force_at_surface = gravitational_constant*total_mass*home_planet_mass/home_planet_radius**2
force_wanted = gravity_force_at_surface + 100


engine_partitions_needed = force_wanted/partition_force
simulate_launch(engine_partitions_needed, partition_force, particles_per_sec_per_partition)
