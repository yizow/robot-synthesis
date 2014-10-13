# Contains functions for the input and output of data to the program
# Reading/writing test results to a file, so we dont have to recompute every time program is run. 
# Filter traces to only keep ones above a certain threshold of points. 

def printResults(results):
	with open('results.txt', 'w') as f:
		for items in results:
			s = ""
			for coordinate in items:
				s += "%s:" % coordinate
			s = s[:-1] + "\n" # drop the last colon and add new line
			f.write(s)

def loadResults():
	with open('results.txt', 'r') as f:
		results = []
		
		while True:
			line = f.readline()
			
			# File End. Return
			if line == '':	
				return results

			trace = []
			for coordinate in line.split(':'):
				# Remove newline character, brackets, and commas; then split each coordinate into values
				coordinate = coordinate[:-1].translate(None, '[],').split()
				# Convert to floats and transform back into start-end pairs; then add to trace
				pairs = []
				for i in range(len(coordinate)/6):
					index = i * 6
					try:
						for add in range(6):
							coordinate[index + add] = float(coordinate[index+add])
					except ValueError:
						print coordinate[i]
					startEnd = [[coordinate[index], coordinate[index+1], coordinate[index+2]], [coordinate[index+3],coordinate[index+4],coordinate[index+5]]]
					pairs += [startEnd]

				trace += [pairs]

			results += [trace]

	return results

def filterResults(results):
	"""We only want traces that have more than a minimum number of points and are almost a closed loop
	"""
	with open('results_filtered.txt', 'w') as f:
		for trace in results:
			if len(trace) > 80 and closeEnough(trace):
				s = ""
				for coordinate in trace:
					s += "%s:" % coordinate
				s = s[:-1] + "\n" # drop the last colon and new line
				f.write(s)

def closeEnough(trace, distance = 1.0):
	items = getMids(trace)
	length = len(items) - 1
	dx = items[0][0] - items[length][0]
	dy = items[0][1] - items[length][1]
	dz = items[0][2] - items[length][2]
	d = sqrt(dx**2 + dy**2 + dz**2)
	if d > distance:
		return False
	return True				

def printPoI(PoI):
	with open('PoI.txt', 'w') as f:
		for items in results:
			s = ""
			for coordinate in items:
				s += "%s:" % coordinate
			s = s[:-1] + "\n" # drop the last colon and add new line
			f.write(s)