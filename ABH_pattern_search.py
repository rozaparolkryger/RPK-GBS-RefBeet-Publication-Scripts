##########################################
### Roza Parol-Kryger                  ###
### parolkrr@cebitec.uni-bielefeld.de  ###
### v0.1                               ###
##########################################


# import
import sys, re, os


def read_file(abh_file, header=False):
	fi = open(abh_file)
	if header:
		d = fi.readline()
		#print "This is header:\t", d
	
# def of global variables, directed to main ()
	
	global corepatterns	
	global querypatterns 
	
	query_name = []
	query_pattern = []
	core_name = []
	core_pattern = []
	corepattern = []
	querypattern = []
	corepatterns = []
	querypatterns = []
	
	while 1:
		l = fi.readline().strip().split()
		if not l:
			break
		elif l == "":
			break		
		else:
			un = re.compile('\ABvchr[1-9].sca')
			if not un.findall(l[0]):
				query_name.append(l[0])
				query_pattern.append(''.join(map(str, l[1:])))
			else:
				core_name.append(l[0])
				core_pattern.append(''.join(map(str, l[1:])))
				
	#print len(query_name)			
			
	for jj in range(len(query_name)):
		querypattern.append([query_name[jj],query_pattern[jj]])
		querypatterns = sorted(querypattern, key=lambda tup: tup[0])
		
	
	#print len(core_name)
			
	for ii in range(len(core_name)):
		corepattern.append([core_name[ii],core_pattern[ii]])
		corepatterns = sorted(corepattern, key=lambda tup: tup[0])
		
		
	#print querypatterns
		
def close(abh_file):
		fi.close()
		fi = None
						
###########################
#Bvchr1.sca001_01	H	B	H	H	B	H	H	H	B	A	H	H	H	H	H	B	A	A	B	A	H	A	H	B	B	H	H	H	B	B	H	A	A	A	B	H	H	H	A	A	A	A	B	B	A	H	A	H	A	H	A	H	H	B	A	H	H	H	A	B
#Bvchr1.sca001_02	H	B	H	H	B	H	H	H	B	A	H	H	H	H	H	B	A	A	B	A	H	A	H	B	B	H	H	H	H	B	H	A	A	A	B	H	H	H	A	A	A	A	B	B	A	H	A	H	A	H	A	H	H	B	A	H	H	H	A	B
######################
# CLASS AND FUNCTION #
######################

#scoring: penalties with edit distance(treshold max60-20=40)
scoring_matrix = {
	('A', 'A'): +1,
	('A', 'B'): 0,
	('A', 'H'): 0,
	('A', '-'): 0,
	('A', '?'): 0,
	('B', 'A'): 0,
	('B', 'B'): +1, 
	('B', 'H'): 0,
	('B', '-'): 0,
	('B', '?'): 0,
	('H', 'A'): 0,
	('H', 'B'): 0,
	('H', 'H'): +1,
	('H', '-'): 0,
	('H', '?'): 0,
	('-', 'A'): 0,
	('-', 'B'): 0,
	('-', 'H'): 0,
	('-', '-'): 1,
	('-', '?'): 0,
	('?', 'A'): 0,
	('?', 'B'): 0,
	('?', 'H'): 0,
	('?', '-'): 0,
	}

# variante liste
def compare(querypatterns, corepatterns, scoring_matrix):
	scoring = []
	for i in range(len(corepatterns)):
		score_list = map(lambda p, s: scoring_matrix[(p, s)], corepatterns[i][1], querypatterns)
		score_sum = reduce(lambda a, b: a+b, score_list, 0)
		scoring.append([score_sum, corepatterns[i][0], corepatterns[i][1]])
	scoring.sort()
	scoring.reverse()
	return scoring


###########
# EXAMPLE #
###########

# variante liste
#patternlist = [
	#('Chr1_sc001', 'ABAH-BBA'),
	#('Chr2_sc003', 'ABBBAHH-'),
	#('Chr3_sc006', 'AB--BBAA'),
	#('Chr4_sc002', 'ABHHBABA'),
	#('Chr5_sc003', 'ABAA-HHB'),
	#('Chr6_sc001', 'ABHBBAAB'),
	#]

#searchpatterns = [
	#('Chr1_un_sc001', 'ABAH-BBA'),
	#('Chr2_un_sc003', 'ABAAAHH-'),
	#]

def main():

		if os.path.isfile(sys.argv[-1]):
			abh_file = sys.argv[-1]
		else:
			exit_on_error(1, sys.argv[-1])

		### do processing
		read_file(abh_file)
		
		
		lq = len(querypatterns)
		lc = len(corepatterns)
		all = lq+lc
		print "all pattern:", all
		print "query patterns:", lq
		print "core patterns", lc
		
	
		
		### Output only names of search sequences 
		#print "###### query sequences #####"
		#for s in range(lq):
			#print querypatterns[s][0],"\t".join(list(querypatterns[s][1]))
			#print querypatterns[s][0:1]
		
		### Output of only pattern names
		#print "###### core patterns  #####"
		#for s in range(lc):
			#print corepatterns[s][0],"\t".join(list(corepatterns[s][1]))
			#print corepatterns[s][0:2]
			
		
		
		
		### counters of diff. matches
		counter_perfect = 0
		counter_optimal = 0
		counter_suboptimal = 0
		counter_muell = 0
		perfect_list = {}
		optimal_list = {}
		suboptimal_list = {}
		
		#cutoff for 1.matrix=55; 2.matrix=120; 3.matrix=40)
		
		cutoff = 57 #60
		sc = 48
		
		#cutoff = 47 #50 
		#sc = 40
		
		#cutoff = 38 #40
		#sc = 32
		
		#cutoff = 28 #30
		#sc = 24
		
		#cutoff = 19 #20
		#sc =16
		
		#cutoff = 9 #10
		#sc = 8
		
		for [sp_name, sp] in querypatterns:
			scoring = compare(sp, corepatterns, scoring_matrix)  # for each pattern one scoring will be created: together for both 20 _un.scaffold_start & _un.scaffold_end
			found_perfect_match = False
			found_optimal_match = False
			for [score, p_name, p_pattern] in scoring:  # for each entry in scoring: [114, 'Bvchr1.sca006_start', 'AAAAhBhAAAhhhhhhhAhAhBBhBAAhAB']
				if score > 0:
					if sp == p_pattern:
						#print sp_name
						#print 'PERFECT_MATCH', sp_name, sp, "is on", p_name, p_pattern, score
						#print 'PERFECT_MATCH', sp_name,  "is on", p_name, p_pattern, score
						counter_perfect = counter_perfect + 1
						found_perfect_match = True
						#erstelle ein eintrag in perfect_list:
						if sp_name in perfect_list:
							perfect_list[sp_name].append(p_name)
						else:
							perfect_list[sp_name]=[p_name]
						
					else: 
						if not found_perfect_match:
							if score >= cutoff:
								#print sp_name
								#print 'OPTIMAL_MATCH', sp_name, "is on", p_name, p_pattern,score
								#print 'OPTIMAL_MATCH', sp_name, sp, "is on", p_name, p_pattern,score
								counter_optimal = counter_optimal + 1
								found_optimal_match = True
								if sp_name in optimal_list:
									optimal_list[sp_name].append(p_name)
								else:
									optimal_list[sp_name]=[p_name]
								
				
							elif (score >= sc and score < cutoff):
								if not found_optimal_match:
									#print sp_name
									#print 'SUBOPTIMAL_MATCH', sp_name, "is on", p_name, score
									#print sp_name,"\t",score
									#print 'SUBOPTIMAL_MATCH', sp_name, sp, "is on", p_name, p_pattern, score
									counter_suboptimal = counter_suboptimal +1
									if sp_name in suboptimal_list:
										suboptimal_list[sp_name].append(p_name)
									else:
										suboptimal_list[sp_name]=[p_name]
				else: 
					break


#######################
# Perfect matches 
#######################		
		#for i in perfect_list:
			#print i, perfect_list[i]
			
		#Output for uniqe und multimatches:
		counter_perf_unique = 0
		counter_perf_multi = 0
		for i in perfect_list:
			if len(perfect_list[i])<2:
				counter_perf_unique = counter_perf_unique+1
			elif len(perfect_list[i])==2:
				#print i, perfect_list[i]
				if perfect_list[i][0][:-3]!=perfect_list[i][1][:-3]:
					#print i, perfect_list[i]
					counter_perf_unique = counter_perf_unique+1			
			else:
				counter_perf_multi = counter_perf_multi+1
		
		print "Unique perfect matches",counter_perf_unique
		print "Multimatches of perfect matches",counter_perf_multi

#######################
# Optimal matches 
#######################		
		#for i in optimal_list:
			#print i, optimal_list[i]
			
		#Output for uniqe und multimatches:
		counter_opt_unique = 0
		counter_opt_multi = 0
		for i in optimal_list:
			if len(optimal_list[i])<2:
				counter_opt_unique = counter_opt_unique+1
			elif len(optimal_list[i])==2:
				#print i, optimal_list[i]
				if optimal_list[i][0][:-3]!=optimal_list[i][1][:-3]:
					#print i, optimal_list[i]
					counter_opt_unique = counter_opt_unique+1			
			else:
				counter_opt_multi = counter_opt_multi+1
		
		print "Unique optimal matches",counter_opt_unique
		print "Multimatches of optimal matches",counter_opt_multi
		
#######################
# Suboptimal matches 
#######################		
		#for i in suboptimal_list:
			#print i, suboptimal_list[i]
			
		#Output uniqe und multimatches:
		counter_subopt_unique = 0
		counter_subopt_multi = 0
		for i in suboptimal_list:
			if len(suboptimal_list[i])<2:
				counter_subopt_unique = counter_subopt_unique+1
			elif len(suboptimal_list[i])==2:
				#print i, suboptimal_list[i]
				if suboptimal_list[i][0][:-3]!=suboptimal_list[i][1][:-3]:
					#print i, suboptimal_list[i]
					counter_subopt_unique = counter_subopt_unique+1			
			else:
				counter_subopt_multi = counter_subopt_multi+1
		
		print "Unique suboptimal matches",counter_subopt_unique
		print "Multimatches of suboptimal matches",counter_subopt_multi
				
	
		#print ""
		#print "Number of search sequences:", ll
		#print "Number of used patterns (containing start and end patterns):", lp
		#print '\n'
		#print "Count of perfect matches:", counter_perfect
		#print "Count of optimal matches:", counter_optimal
		#print "Count of suboptimal matches:", counter_suboptimal
		#print "Count of junky sequences (not possible to anchor) :", counter_muell
		#print ""
		
			
if __name__ == "__main__":
	main()
	

