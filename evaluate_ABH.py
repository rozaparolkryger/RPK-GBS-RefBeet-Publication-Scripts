##########################################
### Roza Parol-Kryger                  ###
### parolkrr@cebitec.uni-bielefeld.de  ###
### v0.1                               ###
##########################################


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# evaluate_ABH_physical_distance

__version__ = "v1.0"
__version_data__ = "20140428"
__usage__ = """

USAGE:
	evaluate_ABH <options> <ABHfilename>

	Scans fastafile for ids and coordinates and prints extracted 
        sequence to STDOUT.

  OPTIONS:

     -h     (help:) display this message.

     FILTERING:

     -ws     (windowsize:) size of the sliding window in bp (default:10000 = 10kb).

     -st     (step:) overlap of the sliding windows in bp (default:5000 = 5kb,
            must be shorter than the windowsize).

     -A     (ABH header:) the ABH table contains a header line, which 
            should be ignored.

     OUTPUT OPTIONS:

     -B     (basename:) basename for the outfile (default:basename of 
            infile).

     OTHER OPTIONS:

     -v     (verbose:) display processing information.

Script version: %s (%s)
""" % (__version__, __version_data__)

#  python ~/Skripts/evaluate_ABH.py -ws 100000 -st 50000 ../Covered_ABH_Chr1.txt 

#################
# CONFIGURATION #
#################

# default values for commandline options
options = {
	'windowsize': 10001,
	'step': 5000,
	'with_header': False,
	'abh_header': False,
	'verbose': False,
	'basename': '',
	}


##########
# IMPORT #
##########

import sys, re, os, pdb, cPickle as pickle, numpy


###############################
# CONSTANTS, GLOBAL VARIABLES #
###############################

# formatted string for error message
__error_msg__ = "ERROR %02d: %s\n"

# defined errors
__errorcodes__ = {
	0: "unrecognized error",
	1: "'%s' is not a file",
	2: "wrong or missing option for '%s'",
	13: "need filename to process",
	41: "wrong id definition: '%s'",
	100: "not yet supported"
	}


#########################
# CLASSES AND FUNCTIONS #
#########################

# Uses an ABH-table (tab-delimited) in the format
# ID POLYMORPHIC_POSITION GENOTYPE1_ABH GENOTYPE2_ABH ...
# The iterator returns all lines with the same ID as
# a list.


class ABHIterator:
	def __init__(self, filename, header=True):
		self.filename = filename
		self.fi = open(self.filename)
		if header:
			d = self.fi.readline()
		self.c_id = ""
		self.cache = []
	def next(self):
		if not self.fi:
			return None
		while 1:
			l = self.fi.readline()
			if not l:
				self.close()
				return self.cache
			l = re.sub("#.*", "", l.strip())
			if not l:
				continue
			#line split and strip
			pp = l.strip().split()
			if pp[0] != self.c_id:
				self.c_id = pp[0]
				tmpr = self.cache[:]
				self.cache = [pp]
				if tmpr:
					return tmpr
			else:
				self.cache.append(pp)
	def close(self):
		self.fi.close()
		self.fi = None

# error reporting 


def exit_on_error(number, *args):

	sys.stderr.write("evaluate_ABH_version=%s\n" % __version__)
	if args:
		error_txt = __errorcodes__[number] % tuple(args)
	else:
		error_txt = __errorcodes__[number]
	sys.stderr.write(__error_msg__ % (number, error_txt))

	sys.exit()

# clip off the first two fields (id and position) of each line
# No! it`s important!

####################################################################################################
def convert_abh_scaffold(sc):
	abh_lines, snp_pos = [], []
	for pp in sc:
		abh_lines.append("".join(pp[2:]))
		snp_pos.append(pp[1])
	return abh_lines, snp_pos
#?
getdiffs = lambda s1, s2: filter(lambda (m, c): ('-' not in (m, c) and m != c), zip(s1, s2))
listsum = lambda l: reduce(lambda a, b: a + b, l, 0)
distlist = lambda m, l: map(lambda v: len(getdiffs(m, v)), l)

####################################################################################################

		
##############################################################
# def of data writing and fast reading
##############################################################

def store(data, fn):
	pickle.dump(data, open(fn, "wb"))

def load(fn):
	return pickle.load(open(fn, "rb"))
###############################################################
# parse reference lengths and read
##############################################################

def ref(sc_id):
	global reflen
	reflist = load('reflist.p')
	for i in range(len(reflist)):
		if sc_id == reflist[i][0]:
			reflen= reflist[i][1]
			break
	return reflen	

#################################################
# parse cov_pileup and save variable as variable for each scaffold
#################################################

def read_cov(sc_id):
	cf = open("coverage_Bvchr1_all60_f.txt", 'r')
	tmp=[]
	tmp_f = []

	while 1:
		l = cf.readline().strip().split()
		if not l:
			break
		else:
			if l[0]==sc_id:		
				tmp_f.append(l)

	store(tmp_f, 'tmp_file.p')
	del tmp_f
	cf.close()
############################################
#this part goes to chec_cov fuction
############################################
def find_cov(sc_id, s, e):
	read_cov(sc_id)
	#print "ich lese coverage..."
	cov = []
	tmp_f=load('tmp_file.p')	
	for i in range(len(tmp_f)):
		if tmp_f[i][0] == sc_id:
			if int(tmp_f[i][1]) == s:
				if int(tmp_f[i][2])==e:
					#print "eintrag gefunden fuer", s, e
					cov = tmp_f[i][3:]
	
	print sc_id, s, e, cov
	return cov

def check_cov(sc_id, s, e):
	#pattern_aa = 60*['a'] 
	patterns = []
	cov = find_cov(sc_id, s, e)
	
	for i in range(len(cov)):
		if float(cov[i])> 0:
			patterns.append('a')
		else:
			patterns.append('-')
	
	pattern = "".join(patterns)		
	return pattern	

####################################################################################################

# evaluate the "average pattern" in a list of ABH-lines

#pattern = check_cov(sc_id,s,e, cov)

def pfreq2(l, p, ws, s, e, sc_id):
## l is abh lines
## p is snp pos list of those l lines
## ws windowsize
## s step	
	#print "wir haben ein fenster mit", len(p), " SNPs erwischt"
	tmp = []
	pattern = []
	# snpfreq wird berechnet, zuerst ist es oben definiert
	snpfreq = float(len(p))/ws
	#print "snpfreq is ", snpfreq "%.2f" % a
	# process each column
	# 1 SNPs in 10 kb = 0.0001 adjust value for Bvchr!!!
	# 0.01 SNPs in 1 kb = 0.00001 adjust value for Scaffolds!!!
	# 1Mb 6 SNPs freq 0.000006
	# 500kb 6 SNP freq 0.000012
	# 300kb (6 SNPs) freq < 0.00002
	# 100kb (10 SNPs) freq < 0.0001
	# 10kb (10 SNPs) freq < 0.001
	if snpfreq < 0.000012:
		tmp = 60*'a'
	else:
		for i in range(len(l[0])):
			# get the values of a column in l line
			cc = map(lambda p: p[i], l)
			# count the occurence of the values
			abhv = [cc.count('A'), cc.count('B'), cc.count('H'), cc.count('-')]
			# numbers of positions with information
			wsi = len(p) - abhv[3]
			# wsi = 72 - 16 = 56
			# assignment to categories
			if abhv[3] > len(p)*0.6:
				# not enough information
				tmp.append('-')
			elif abhv[0] > wsi*0.9:
				# at least 80% are A's
				tmp.append('A')
			elif abhv[1] > wsi*0.8:
				# at least 80% are B's
				tmp.append('B')
			#elif abhv[2] > wsi*0.2:
				# at least half are H's
				#tmp.append('H')
			else:
				# unclear distribution
				tmp.append('H')
	return [s, e, len(p), "%.5f" % snpfreq, "".join(tmp)]

############################################################################
# process each scaffold with a sliding window
def walk2(abh_lines, snp_pos, reflen, fo1, fo3, options, sc_id):	
	
	#
	# adaption of windowsize, if scaffold is shorter than the default
	ws = options['windowsize']
	st = options['step']
	
	# set some values
	#output
	start_pattern = None
	end_pattern = None
	all_lines = []
	
	#containers
	tmp_lines=[]
	tmp_pos=[]
	tmp_lines_all=[]
	tmp_pos_all=[]
	
	
	m = ''
	i = 1
	snp_pos_i = 0
	snplistlen = len(snp_pos)
	
	if ws > int(reflen):
		print "The length of ", sc_id, reflen, "is smaller than windowsize ", ws
		print "The windowsize will be set to reference length", reflen
		ws = int(reflen)
		st = int(reflen)-2
	
	while i+st < int(reflen):
		# processing each window 
		s, e, next_s = i, min(reflen, i+ws-1), i+st
		window_snp = []
		window_abh = []
		#print "WINDOW %s-%s" % (s, e)
		#print "Inspecting SNPs from %s to %s" % (snp_pos_i, snplistlen)
		for j in range(snp_pos_i, snplistlen, 1):
			# processing window_snp with relevant SNPs
			if s > int(snp_pos[j]):
				# sollte NICHT auftreten!
				snplist_i = j
				continue
			elif int(snp_pos[j]) <= e:
				window_snp.append(snp_pos[j])
				window_abh.append(abh_lines[j])
				if int(snp_pos[j]) < next_s:
					snp_pos_i = j+1
			else:
				break
		#print "SNPs in window: %s" % window_snp
		#print "abh in windows: %s" % window_abh
		if len(window_snp) > 0:
			tmp_lines= pfreq2(window_abh, window_snp, ws, s, e, sc_id)
			
		else:
			snpfreq = 0
			#pattern = check_cov(sc_id, s, e)
			#tmp_lines = [s, e, len(window_snp), snpfreq, pattern]
			tmp_lines = [s, e, len(window_snp), "%.2f" % snpfreq, 60*'a' ]
		
		tmp_lines_all.append(tmp_lines)	
		i += st
		
		#write in file
		if fo1:
			fo1.write("%s\n" % ("\t".join(map(str,tmp_lines))))
		elif verbose:
			#print "\t".join(map(str,tmp_lines"%s(%s)" % ))
			print '.'
		if fo3:
			fo3.write("%s\t%s\n" % (("\t".join(map(str,tmp_lines[0:4]))),("\t".join(map(str,list(tmp_lines[4]))))))
			#list('abcdb')
			#'/t'.join(list('abcgd')'
		elif verbose:
			#print "\t".join(map(str,tmp_lines"%s(%s)" % ))
			print '.'
			
		#start_pattern = tmp_lines_all[0]
		#end_pattern = tmp_lines_all[-1]
			
		w = 1
		y = 3
		#v = 1
		if len(tmp_lines_all) <= 2:
			start_pattern = tmp_lines_all[0]
			end_pattern = tmp_lines_all[-1]
			#mid_pattern = start_pattern
		else:
			start_pattern = tmp_lines_all[0]
			while 'a' in list(start_pattern[4]) and w<len(tmp_lines_all):
				start_pattern = tmp_lines_all[w]
				w+=1
			
			end_pattern = tmp_lines_all[-2]
			while 'a' in list(end_pattern[4]) and y<len(tmp_lines_all):
				end_pattern = tmp_lines_all[-y]
				y+=1
				
			#mid_pattern = 	tmp_lines_all[int(len(tmp_lines_all)/2)]
			#while 'a' in list(mid_pattern[4]) and y<len(tmp_lines_all):
				#mid_pattern = tmp_lines_all[int(((len(tmp_lines_all))/2)+v)]
				#v+=1
				
				
		#except IndexError: 
			#print "Error: scaffold to short:", sc_id
			#end_pattern = tmp_lines_all[-1]
	print "The length of",sc_id, "is",reflen ,"with",len(snp_pos),"SNPs. # evaluated windows is	:",len(tmp_lines_all)	
	
		
	if not end_pattern:
		end_pattern = m
	
	return tmp_lines_all, start_pattern, end_pattern

# tmp_lines_all:
#[[1, 3000, 0.0013333333333333333, 'hBAABHhhBAhHhhhBAABAAAhBBhhhBBhAAABAhhAAAABBAAAhAhAhhBAhHhAB'],
#[1001, 4000, 0.002, 'hBAABHhhBAhHhhhBAABAhAhBBhHhBBhAAABAAAAAAABBAhAhAhAhhBAAhhAB']]

#  - <basename>.sc_report.ws<windowsize>.dv<dv_cutoff>.csv, containing part of the
#    data in csv format, each value of ABH-data in a separate column.
#     # Bvchr1.sca001 (16118 SNPs)
#     0 	501	A	A	A	A	A	A	A	A	A	H	H	A	A	A	H
#     100	501	A	A	A	A	H	A	A	A	A	H	H	A	A	A	A
#     200	501	A	A	A	A	A	A	A	A	A	H	H	A	A	A	A
#     ...
#
#  - <basename>.sc_overview.ws<windowsize>.dv<dv_cutoff>.txt, showing the first and 
#    last pattern of each scaffold and some additional processing data.
#     Bvchr1.sca001	16118	AAAAAAAAAHHAAAH	0	AAAAAAAAAHHAAAA	501	501	2020	2831
#     Bvchr1.sca002	7352	AAAAAAAAAHHAAAH	0	HAAAAAAAAHHAAAA	501	501	1896	2310
#     Bvchr1.sca003	4424	HAAAAAAAAHHAAAA	0	AAAAAAAAAHHAAAA	501	501	1846	2189
#
##########################################################################################################3
def evaluate_patterns_per_scaffold(abh_file, options):
	########################################
	# open file handles, write header lines
	########################################
	# create File 1 (txt with pattern)
	fno = "%(basename)s.report.ws%(windowsize)s.st%(step)s.txt" % options
	if options['verbose']:
		print "OPENING FILEHANDLE '%s'" % fno
	fo1 = open(fno, 'w')
	if fo1:
		fo1.write("WIN_start\tWIN_end\tSnpsNr in win\tPATTERN\n")
	
	# File 2 with start und stop pattern
	if options['verbose']:
		print "OPENING FILEHANDLE '%s'" % "%(basename)s.overview.ws%(windowsize)s.st%(step)s.txt" % options
	fo2 = open("%(basename)s.overview.ws%(windowsize)s.st%(step)s.txt" % options, 'w')
	fo2.write("SCAFF_ID\tWin_start\tWin_end\tSTART_P\tEND_P\tSNP_Freq(?)\n")
	
	#file 3 = file 1 with abh lines tab sep
	if options['verbose']:
		print "OPENING FILEHANDLE '%s'" % "%(basename)s.report.ws%(windowsize)s.st%(step)s.csv" % options
	fo3 = open("%(basename)s.report.ws%(windowsize)s.st%(step)s.csv" % options, 'w')

	#####################################################
	# open ABH iterator, read and parse scaffoldwise
	#####################################################
	abh = ABHIterator(abh_file, options['abh_header'])

	while 1:

		# scaffold-wise iteration
		rr = abh.next()
		if not rr:
			break
		sc_id = rr[0][0]

		#reflen parsen
		reffile = open("RefBeet-1.2_lengths.txt", 'r')
		reflist = []

		while 1:
			l = reffile.readline().strip().split()
			if not l:
				break
			else:
				reflist.append([l[0],l[1]])

		store(reflist, 'reflist.p')
		del reflist
		reffile.close()
		
		reflen = ref(sc_id) 
		
		# prepare rr for following processing and write scaffold name to outfiles
		abh_lines, snp_pos = convert_abh_scaffold(rr)
		fo1.write("# %s (%s SNPs) (length: %s) \n" % (sc_id, len(abh_lines), reflen))
		fo3.write("# %s (%s SNPs) (length: %s) \n" % (sc_id, len(abh_lines), reflen))

		
		
		# process rr with a sliding window	
		sc_r = walk2(abh_lines, snp_pos,reflen, fo1, fo3, options, sc_id)
		
		# write processing data to overview file
		#fo2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % sc_id, len(abh_lines)] + list(sc_r)))
		fo2.write("%s\t%s\t%s\t%s\t%s\n" % (sc_id, reflen,len(abh_lines), sc_r[1], sc_r[2]))
		
		# report processing progress
		if options['verbose']:
			print ".",

	# close file handles
	fo1.close()
	fo2.close()
	fo3.close()



##############
# PROCESSING #
##############

def main():

	#try:

		# help needed?

		if len(sys.argv) == 1 or '-h' in sys.argv:
			sys.stdout.write(__usage__)
			sys.exit()

		# get / check infiles

		if len(sys.argv) < 2 or sys.argv[-1].startswith("-"):
			exit_on_error(13)

		if os.path.isfile(sys.argv[-1]):
			abh_file = sys.argv[-1]
		else:
			exit_on_error(1, sys.argv[-1])

		# parameters for processing

		if '-A' in sys.argv:
			options['abh_header'] = True

		if '-ws' in sys.argv:
			try:
				options['windowsize'] = int(sys.argv[sys.argv.index('-ws')+1])
			except (KeyError, ValueError):
				exit_on_error(2, '-ws')

		if '-st' in sys.argv:
			try:
				options['step'] = int(sys.argv[sys.argv.index('-st')+1])
				if options['step'] >= options['windowsize']:
					exit_on_error(2, '-st')
			except (KeyError, ValueError):
				exit_on_error(2, '-st')

		if '-B' in sys.argv:
			try:
				options['basename'] = sys.argv[sys.argv.index('-B')+1]
			except (KeyError, ValueError):
				options['basename'] = os.path.basename(abh_file)
		else:
			options['basename'] = os.path.basename(abh_file)

		# other options

		if '-v' in sys.argv:
			options['verbose'] = True

		# do processing

		if options['verbose']:
			print "Evaluating patterns by sliding window..."
		
		evaluate_patterns_per_scaffold(abh_file, options)
		if options['verbose']:
			print


### main ###

if __name__ == "__main__":
	main()

