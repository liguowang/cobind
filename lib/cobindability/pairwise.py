import sys,os

def find_bed_files (root_path,verbose=False):
	'''return list of bed files'''
	
	bed_files = []
	for root, dirs, files in os.walk(root_path,followlinks=True):
		for f in files:
			if f.endswith('.BED'):
				bed_files.append(f)
	return bed_files
all_bed_files = find_bed_files(sys.argv[1])

j = 0
for i in range(0, len(all_bed_files)):
	j = i +1
	for n in range(j,len(all_bed_files)):
		print "python2.7 overlap_coefficient.py " + all_bed_files[i] + ' ' + all_bed_files[n] + " >>output" 	