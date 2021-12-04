import sys
from numpy import mean, median, std


infile = sys.argv[1]
size = []
total_size = 0
for l in open(infile,'r'):
	l = l.strip()
	f = l.split()
	chrom = f[0]
	start = int(f[1])
	end = int(f[2])
	total_size += (end -start)
	#size.append(end - start)
	
#print '\t'.join([infile,str(len(size)), str(mean(size)), str(median(size)), str(std(size))])	
print total_size	
