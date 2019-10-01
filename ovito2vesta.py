# coding=utf-8
import sys

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Error: enter the file name")
		sys.exit()
	fileName = sys.argv[1]
	
	with open(fileName, 'r') as f:
		coords = []
		natom = 0
		for line in f.readlines():
			words = line.strip().split()
			if len(words) >= 3:
				point = [float(words[0]), float(words[1]), float(words[2])]
				coords.append(point)
				natom += 1

	ofileName = ""
	if len(sys.argv) == 2:
		print("Using the default output file name: MoS2.xyz")
		ofileName = "MoS2.xyz"
	else:
		ofileName = sys.argv[2]

	with open(ofileName, 'w') as f:
		f.write("%d\n"%natom)
		for point in coords:
			if point[2] != 0:
				f.write("S\t")
			else:
				f.write("Mo\t")
			for coord in point:
				f.write("%12.6f\t"%coord)
			f.write("\n")


