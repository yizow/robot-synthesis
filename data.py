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

