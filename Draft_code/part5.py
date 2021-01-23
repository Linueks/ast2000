from ast2000solarsystem_27_v4 import AST2000SolarSystem

star_system = AST2000SolarSystem(11466)

star_system.engine_settings(3.4021029823078977e-09, 3.15498195925e14, 63200999999999.992, 100000, 1120,  [4.18992917, -1.46374675], 11.3)
star_system.mass_needed_launch([ 4.18998484, -1.46357376], test=True)

star_system.send_satellite("instructions.txt")


# 53
# x: 0.69970958 y: 5.21180747
# Input x-position: 4.42711
# Input y-position: -0.435386
