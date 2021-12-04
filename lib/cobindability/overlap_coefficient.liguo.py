import sys,os
import math
from warnings import warn
from bx.bitset import *
from bx.bitset_builders import *
from bx.cookbook import doc_optparse


def find_bed_files (root_path,verbose=False):
	'''return list of bed files'''
        
	bed_files = []
	for root, dirs, files in os.walk(root_path,followlinks=True):
		for f in files:
			if f.endswith('.BED'):
				bed_files.append(f)
	return bed_files

def intersectBed2(file1,file2):
	"""
	calculate total number of bases overlapped between 2 bed files
	"""
	overlap_size = 0.0
	bits1 = binned_bitsets_from_file( open( file1) )
	bits2 = binned_bitsets_from_file( open( file2) )

	bitsets = dict()
	#if file1 == file2:
	#	return 1.0
	for key in bits1:
		if key in bits2:
			bits1[key].iand( bits2[key] )
			bitsets[key] = bits1[key]

	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0	
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			overlap_size += end - start
			#print "%s\t%d\t%d" % ( chrom, start, end )
	return overlap_size


def total_bases(bedfile1, bedfile2):
	bed_file1_size = 0.0
	bed_file2_size = 0.0
	
	#print >>sys.stderr, "reading %s ..." %  bedfile1
	for l in open(bedfile1,'r'):
		l = l.strip()
		f = l.split()
		bed_file1_size += (int(f[2]) - int(f[1]))
	
	#print >>sys.stderr, "reading %s ..." %  bedfile2
	for l in open(bedfile2,'r'):
		l = l.strip()
		f = l.split()
		bed_file2_size += (int(f[2]) - int(f[1]))
	
	
	return 	(bed_file1_size, bed_file2_size)

"""
#############################################
#
# usage: python2.7 overlap_coefficient.py .
#
##############################################



# get all bed files from current folder
all_bed_files = find_bed_files(sys.argv[1])


print "Name\t" + '\t'.join(all_bed_files)

for f in all_bed_files:
	print f + '\t',
	tmp = []
	for ff in all_bed_files:
		file1_size, file2_size = total_bases(f,ff)
		overlap = intersectBed2(f,ff)
		tmp.append(overlap/min(file1_size, file2_size))
	print '\t'.join([str(i) for i in tmp])
		

##test
"""


file1_size, file2_size = total_bases(sys.argv[1],sys.argv[2])
overlap = intersectBed2(sys.argv[1],sys.argv[2])
print(sys.argv[1] + '\t' + sys.argv[2] + '\t' + str(overlap/math.sqrt(file1_size*file2_size)))
