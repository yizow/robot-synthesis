from robot import *

class fourBar(object):

	"""
	Give the orientations of the two cranks. Coupler is calculated. 

	angles should be between 0 and 180 degrees. 
	angles given are measured CCW as in polar angles. 
	"""
	def __init__(self, length1, angle1, length2, angle2, distance):
		
		# Make sure lengths are positive
		if length1 <= 0 or length2 <= 0:
			print "Lengths must be positive"
			raise AttributeError

		# Ensure that angles are between 0 and 180 to prvent intersection of beams
		angle1 %= pi
		angle2 %= pi
		if (angle1 >= pi or angle1 <= 0) or (angle2 >= pi or angle2 <= 0):
			print "Angles must be between 0 and pi"
			raise AttributeError

		l1 = Link(0,length1,0,0)
		l2 = Link(0,length2,0,0)
		base = Link (0,distance,0,0)

		# TODO: 
		# For now, I'll be ignoring the intersection of coupler with beams
		# I will probably add a function that evaluates robot state for intersections

		# Calculate length of coupler
		# Position of crank endpoint
		pos1 = [length1 * trig for trig in (cos(angle1), sin(angle1))]
		# Position of other beam endpoint
		pos2 = [length2 * trig for trig in (cos(angle2), sin(angle2))]
		pos2[0] += distance
		dx = pos2[0]-pos1[0]
		dy = pos2[1]-pos1[1]

		# Pythagorean formula
		couplerLength = sqrt(dx**2 + dy**2)
		coupler = Link(0,couplerLength,0,0)

		# Calculate angles between cranks and coupler. 
		self.turn1 = -(angle1 + arctan(dy/dx))
		self.turn2 = pi + angle1 - angle2 - self.turn1

		# Construct the robot
		self.fourBar = Robot([l1, coupler, l2, base], name = 'fourBar', manuf = 'Yizhe Liu', comment = 'test')
		self.l1 = l1
		self.l2 = l2
		self.coupler = coupler
		self.base = base
		self.angle1 = angle1
		self.angle2 = angle2

		return None

	def __repr__(self):

		return "angle1=%f, length1=%f, turn1=%f, coupler=%f, turn2=%f, length2=%f, angle2=%f, base=%f" % (self.angle1, self.l1.A, self.turn1, self.coupler.A, self.turn2, self.l2.A, self.angle2, self.base.A)

def test():
	f = fourBar(1, pi/4, 1, pi/4, 4)
	print f
	print fkine(f.fourBar, [f.angle1, f.turn1, f.turn2, -f.angle2])
