################################import statements###############################
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from optparse import OptionParser

################################################################################
##############################Function Definitions##############################
################################################################################

##################################Program Admin#################################

################################################################################
#Option Parsing
################################################################################
def do_args(parser):
	'''A function that parses the command-line arguments supplied by the user
	and returns the options object and the leftover args list in a tuple.'''


	'''Can't be unit-tested.'''
		
	#define all the options
	##file options
	parser.add_option('--in','--in-file',type='string',help='A required '+\
						'option. Takes the name of the fasta file containing'+\
						' the sequence to annotate. Only the first sequence '+\
						'in the file will be used.',dest='in_file')
	parser.add_option('--out','--out-file',type='string',help='A required '+\
						'option. Takes the name of the file to which the '+\
						'output annotations should be written. Silently over-'+\
						'writes the provided file if it exists, so be careful!'\
						,dest='out')
	parser.add_option('--frmt','--format',type='string',default='csv',help=\
						'Specifies the desired output file type. Currently '+\
						'only csv is supported, so this option is not used, '+\
						'but we hope to eventually implement genbank.',\
						dest='format')
	parser.add_option('--ci','--cass-inner',action='store_true',default=False,\
						help='Use this flag if you want the sequences in the '+\
						'cassette fastas to exclude the attC sites and only '+\
						'include the non-attC parts of the cassettes. Off by '+\
						'default.')

	##attI options
	parser.add_option("--attISeq",help='Allows the user to provide their own '+\
						'attI sequence if they have one that is not one of '+\
						'the ones provided. The attI sequence must end in GTT'+\
						'. Defaults to None.')
	parser.add_option('--nattI','--no-attI',action='store_true',default=False,\
						help='Use this flag to instruct the program not to '+\
						'look for an attI site, but just find attCs. This is '+\
						'best combined with start position (--sp) and '+\
						'direction (--dir) arguments (if you know the values '+\
						'for them) so that the computer doesn\'t waste its '+\
						'time.',dest='nattI')

	##attC options
	parser.add_option('--co1','--score1-cutoff',type='int',help='The minimum '+\
						'allowable score for the structure-based scoring of '+\
						'attC. Must be an integer between 0 and 59. Defaults'+\
						' to 44.',default = 44, dest='co1')
	parser.add_option('--co2','--score2-cutoff',type='int',help='The minimum '+\
						'allowable score for the sequence-matching scoring of'+\
						'attC. Must be an integer between 0 and 50. Defaults '+\
						'to 17.',dest = 'co2',default=17)
	parser.add_option('--attC','--attC-only',action='store_true',default=False,\
						dest='attCo',help='Use this flag when your input '+\
						'fasta file only contains attC sequences (one per '+\
						'entry in the fasta file) and you want each one to be'+\
						'annotated individually. Defaults to False.')
	
	##length options
	parser.add_option('--aml','--attCminLen',type='int',default=50,help=\
						'The minimum allowable length of attC. Must be at '+\
						'least 20. Defaults to 50.')
	parser.add_option('--aMl','--attCmaxLen',type='int',default=150,help=\
						'The maximum allowable length of attC. Defaults to '+\
						'150.')
	parser.add_option('--cml','--cassminLen',type='int',default=0,help=\
						'The minimum allowable length of a cassette, '+\
						'excluding the attCs. Defaults'+\
						' to 0.')
	parser.add_option('--cMl','--cassmaxLen',type='int',default=4000,help=\
						'The maximum allowable length of a cassette, '+\
						'excluding the attCs. 0 means no maximum. Defaults'+\
						' to 4000.')

	##searching options
	parser.add_option('--sp','--start-position',type='int',default=1,help=\
						'The position at which you want the program to start '+\
						'looking for attCs. Defaults to 1 (the first base in '+\
						'the sequence). Interacts with dir such that if dir '+\
						'is set to -1, the value 1 will refer to the first '+\
						'base in the reverse complement of the provided '+\
						'sequence (the last base of the provided sequence).'+\
						' It is best to leave this at 0 if --dir is going to '+\
						'be 0. If a value less than 1 is given, 1 will be '+\
						'used.')
	parser.add_option('--ep','--end-position',type='int',default=None,help=\
						'The position at which you want to program to stop '+\
						'looking for attCs. If not provided, or if the value '+\
						'provided is greater than the length of the sequence '+\
						'(minus --sp) then the program will seach to the end '+\
						'of the sequence. Interacts with --dir like --sp does.')
	parser.add_option('--dir','--direction',type='int',default=0,help=\
						'The direction in which the program will search for '+\
						'attIs (or attCs if --nattI is used). By default the '+\
						'program searches for attI in both directions and '+\
						'the for attC only in the direction where an attI '+\
						'search was successful. If there is no attI found (or'+\
						' if -nattI was used) then use of this option is '+\
						'highly recommended if you know which direction you '+\
						'want to search in, to save computer time. If you '+\
						'have no idea which direction you want to search in '+\
						'then leave this to the default value of 0 and both '+\
						'directions will be searched. If you leave this at 0 '+\
						'it is best not to give values for --sp and --ep. '+\
						'--dir Must be given the value 1 (forward), -1 '+\
						'(reverse), or 0 (both).',dest='direction')

	##behaviour options
#	parser.add_option('--q','--quiet',action=store_true,default=False,help=\
#						'Use of this option is required for batch use of '+\
#						'this software such as running from a Make rule. It'+\
#						' suppresses all prompts to the user.')

	(options,args) = parser.parse_args()

	return (options,args)

################################################################################
#File handling
################################################################################
def open_fasta(fasta,aml):
	'''A function that takes the name of a fasta file, fasta, and returns a
	sequence record object of the first sequence in that fasta file. Raises an
	exception if the sequence is shorter than the minimum attC length.'''

	''' Final Writing and unit tested. '''

	try:
		rec = SeqIO.read(fasta,'fasta')

	except ValueError,e: 
		#if there is more than one record in the file use the first and warn
		if e.args[0] == 'More than one record found in handle':
			recs = SeqIO.parse(fasta,'fasta')
			rec = recs.next()
			print 'The file %s ' % fasta +'contains more than one record. '+\
					'Using only the first.'

			#other exceptions get raised directly
		else:
			raise e
	if len(rec) < aml:
		ex = 'The sequence must be at least as long as the minimum attC length.'
		raise ValueError(ex)

	return rec

def get_seqs(drn,rec):
	'''A function that takes a direction variable, drn, and a sequence record
	object, rec, and returns a list of sequences taken from rec in the
	directions specified by drn (1 for forward, -1 for backward, 0 for
	both.)'''


	''' Final Writing and unit tested. '''

	#Create the output list:
	seqs = []

	if drn == 0:
		seqs.append(str(rec.seq))
		seqs.append(str(rec.seq.reverse_complement()))
	elif drn == 1:
		seqs.append(str(rec.seq))
	elif drn == -1:
		seqs.append(str(rec.seq.reverse_complement()))

	return seqs

def list_attIs(attIfnm,attIseq=''):
	'''A function that takes the file name where the attIs are stored, attIfnm,
	and the attI sequence provided by the user, attISeq, and returns a list of
	properly-formatted and checked attIs that can be used for searching.'''

	''' Final Writing and unit tested. '''

	#fetch default attIs
	'''Default attIs will be in a text file in the directory with the
	rest of the program. To fetch them we open the file and put all the
	lines into a list. The user can add more attIs to this directory
	directly if they want.'''

	attIf = open(attIfnm)
	attIs = attIf.readlines()
	attIf.close()

	#attI sequence provided?
	if attIseq:
		#if yes, add to list
		attIs.append(attIseq)

	#make sure everything is uppercase and stripped
	for i in range(len(attIs)):
		attIs[i] = attIs[i].strip()
		attIs[i] = attIs[i].upper()

	#check the sequences for GTT
	msg = ''
	for i in attIs:
		if i[-3:] != 'GTT':
			msg += 'The attI sequence %s doesn\'t end in GTT. Check the ' % i\
			+'default attI list in %s.' % attIfnm
	
	if msg:
		raise ValueError(msg)

	return attIs
	
################################Program Function################################

################################################################################
#Initial Search
################################################################################
def find_all_attIs(seqstr,attIfnm,aml,attIseq=''):
	'''A function that takes a sequence string (seqstr), a name of a file of
	attIs (attIfnm), and an optional string of an attI sequence (attIseq) and
	returns a list of starting dictionaries for integrons within that
	sequence.'''


	'''Written and unit tested.'''

	attIs = list_attIs(attIfnm,attIseq) #get list of attIs
	this_dir = [] #to store all integrons in this direction

	'''Now that I have a list of all the attIs I'm using, iterate
	through all the attIs in the list and store the position of each one
	in a list of start points.'''

	for attI in attIs: #check each attI sequence
		################################################################
		#iterate through attIs and find each instance of each one
		################################################################
		combo = [] #create lst of output data for this seq/attI combo

		#prep for the while loop below
		start = 0
		attI_pos = 0
		while (start != -1) and (start < len(seqstr)-len(attI)-aml):
			
			#search for the next instance of the attI
			this_int = find_attI(start,attI,seqstr)
			if this_int:
				combo.append(this_int) #store this instance of this attI
				start = this_int['start']
			else:
				start = -1

		this_dir += combo #store all instances of a given attI

	return this_dir

def find_first_attC(seq,mco1,co2,aml,aMl,rbd,bds,final=0):
	'''A function that takes a sequence string (seq), a score1 cutoff modified
	because we're not looking for a full match (mco1), a score2 cutoff (co2),
	the minimum and maximum attC lengths (aml and aMl) and base frequency
	dictionaries (rbd,bds) and searches for a single, intact attC. Returns a
	list containing a single starting dictionary with the first attC.'''


	'''Final Writing and unit tested.'''

	#Look for the first attC
	first_attC = attC_search(seq,0,len(seq),mco1,co2,aml,aMl,rbd,bds,\
								final=final)
	
	#if there was no attC, return 0
	if not first_attC:
		return 0

	#Create the "first_attC" dictionary
	start_d = {}
	for k in first_attC:
		key = 'First'+k
		start_d[key] = first_attC[k]
	
	#add the data needed to find the rest of the attCs if not final
	if not final:
		start_d['start'] = start_d['Firstend']
		start_d['R1'] = start_d['FirstattC'][-4:]

	return start_d


def find_attI(start,attI,seqstr):
	'''A function that takes a start position, start, an attI sequence, attI,
	and a sequence string, seqstr. Searches the seqstr for the first occurence
	of the given attI and returns a dictionary with the keys attIStart,
	attIEnd, attISeq, start, and R1. The value of attISeq will include the
	entire R' it is fused with and R1 will contain the last four bases of that
	R' sequence.'''


	''' Final Writing and unit tested'''

	#look for the attI
	attI_pos = seqstr.find(attI,start)

	#if found
	if attI_pos != -1:
		start = attI_pos + len(attI) + 4

		#create the dictionary for this integron
		this_int = {}
		this_int['attISeq'] = seqstr[attI_pos:start] #the whole attI incl R'
		this_int['attIStart'] = attI_pos
		this_int['attIEnd'] = start
		this_int['R1'] = this_int['attISeq'][-4:]
		this_int['start'] = start #the place to start looking for the attC

	#if not found
	else:
		start = -1
		this_int = None

	return this_int

	
################################################################################
#attC Searching
################################################################################
def find_R(seq,se,beg=0,end=0):
	'''A function that takes a sequence string, seq, an integer representing 
	the index to start searching at, beg, a string, se, that must be either 
	"start" or "end", and an optional argument, end, that gives the index at 
	which to stop looking, and returns a tuple whose first element is the 
	starting index of the first R site it encounters in the sequence and whose 
	second element is that site's sequence. If se is "start" it searches for 
	'AAC' and gives an index four bases before that. If se is "end" it searches 
	for 'GTT' and gives the index of the G.  Raises an exception if se isn't 
	"start", "s", "end", or "e". Returns -1 if no R site was found.'''

	''' Final Writing and unit tested'''

	#Check inputs
	msg = ''
	if (se.lower() != "s") and (se.lower() != "start") and (se.lower() != "e") \
		and (se.lower() != "end"):
		msg += '\nThe second argument of find_R() must be one of \'start\', '+\
				'\'s\', \'end\', or \'e\'.\n'

	if seq == '':
		msg += '\nThe sequence may not be null.\n'

	if beg >= len(seq):
		msg += '\nThe beginning search position must be less than the length '+\
					'of the sequence.\n'
	
	if msg:
		raise ValueError(msg)

	if end:
		if end <= beg+7:
			return (-1,None)
		if end > len(seq):
			end = len(seq)

	#Start or end?
	if (se.lower() == 's') or (se.lower() == 'start'):
		r = 'AAC'
		#the value to adjust the positioning of the R site
		ap = 4
	
	elif (se.lower() == 'e') or (se.lower() == 'end'):
		r = 'GTT'
		ap = 0

	#check if we need to include the stop value
	if not end:
		pos = seq.find(r,beg)
	else:
		pos = seq.find(r,beg,end)

	#is there room for a whole R site in the searched area?
	if pos - ap >= 0:
		pos = pos - ap
		if end:
			if pos + 7 > end:
				#there isn't room for the whole R site
				return (-1,None)
		else:
			if pos + 7 > len(seq):
				#there isn't room for the whole R site
				return (-1,None)
	else:
		#there isn't room for the whole R site
		return (-1,None)

	return (pos,seq[pos:pos+7])

def check_scores(score_dict,co1,co2,bds,attC):
	'''A function that takes a score dictionary, score_dict, an attC sequence,
	attC, the dictionary of base frequency dictionaries for score2 (bds) and
	the two score cutoffs, co1 and co2. Returns the score dictionary with
	score2 added to it if both scores meet the cutorr. Returns a dictionary
	with no score2 value otherwise.'''


	'''Final Writing and not unit tested (score2)'''

	#check first score
	score1 = score_dict['score']
	if score1 >= co1:
		#check second score
		s2 = score2(score_dict,attC,bds)
		if s2 >= co2:
			#add it to the dictionary and return the dictionary
			score_dict['score2'] = s2
			return score_dict
		else:
			return score_dict
	else:
		return score_dict
	
def gtt_search(seqstr,aac,min_end,max_end,co1,co2,rbd,bds,R4='',final=0):
	'''A function that takes a sequence string, seqstr, the starting position
	of the current potential attC (aac), a gtt position for minimum length
	(min_end), a gtt position for maximum length (max_end), a score 1 cutoff
	(co1), a score 2 cutoff (co2), base frequency dictionaries (rbd,bds), an R'
	sequence to match to that has the empty string as a default value if there
	is no R' to match to, and a switch, final, that takes the value 0 if we are
	not searching for a final cassette (default) and 1 if we are. Searches for
	an acceptable attC in the sequence.  Returns the score data for the first
	acceptable attC it finds, unless final is 1, in which case it returns all
	potentially-acceptable attCs. If it doesn't find an acceptable attC,
	returns 0.'''


	''' Final Writing but not unit tested (score2)'''
	
	#Make list of good enough attCs
	ge = []

	#Set gtt to start
	gtt = 0
	while (gtt != -1) and (min_end < max_end):
		'''This loop searches through gtts until either the score cutoff is met
		or the maximum attC length is reached.'''

		if not final:
			#find the next gtt
			gtt,R2 = find_R(seqstr,'e',min_end,max_end)
		elif final:
			#find the next g
			gtt = seqstr.find('G',min_end,max_end)
						
		if gtt != -1:
			#Score the attC (fn not done)
			if not final:
				pot_attC = seqstr[aac:gtt+7]
			else:
				pot_attC = seqstr[aac:gtt+1]
			score_dict = score_whole_attC(pot_attC,rbd,final=final,R4=R4)
			score_dict['dir'] = 1

			#Check the scores against the cutoffs
			score_dict = check_scores(score_dict,co1,co2,bds,pot_attC)
			if 'score2' in score_dict: #this indicates it passed
				#add the start and end positions
				score_dict['start'] = aac
				score_dict['attC'] = pot_attC
				if not final:
					score_dict['end'] = gtt+7
				else:
					score_dict['end'] = gtt+1
				ge.append(score_dict)

			else: #it failed
				#check this attC's reverse complement
				rev_attC = str(Seq(pot_attC).reverse_complement())
				score_dict = score_whole_attC(rev_attC,rbd,final=final,R4=R4,\
												rev=1)
				score_dict['dir'] = -1

				#Check the scores against the cutoffs
				score_dict = check_scores(score_dict,co1,co2,bds,rev_attC)
				if 'score2' in score_dict: #indicates it passed
					#add start and end positions
					score_dict['start'] = aac
					score_dict['attC'] = rev_attC
					if not final:
						score_dict['end'] = gtt+7
					else:
						score_dict['end'] = gtt+1
					ge.append(score_dict)

			#look for the next gtt
			if not final:
				min_end = gtt + 3
			else:
				min_end = gtt+1
	
	#Once out of the while loop, check if any attCs were found
	if len(ge) == 0: #no attCs found
		return 0

	else: #attCs found
		#get the best one
		if len(ge) == 1:
			best_attC = ge[0]
		else:
			best_attC = get_best(ge)

		return best_attC

def get_best(ge):
	'''A function that takes a list of score dicts, ge, and returns the dict
	with the highest score. If there is more than one then the first is
	taken.'''

	
	'''under construction and not unit tested'''

	
	#Get the best scores by score 1
	best = 0
	best_ds = []
	for d in ge:
		if d['score'] == best:
			best_ds.append(d)
		elif d['score'] > best:
			best_ds = [d]
			best = d['score']

	#If there is more than one, sort by score 2
	if len(best_ds) > 1:
		best = 0
		bestest = []
		for d in best_ds:
			if d['score2'] == best:
				bestest.append(d)
			elif d['score'] > best:
				bestest = [d]
				best = d['score2']
		#If there is still more than one, just take the first
		best_attC = bestest[0]

	else:
		best_attC = best_ds[0]

	return best_attC

	
				

def old_gtt_search(seqstr,aac,min_end,max_end,co1,co2,rbd,lbd1,lbd2,R4='',final=0):
	'''A function that takes a sequence string, seqstr, the starting position
	of the current potential attC (aac), a gtt position for minimum length
	(min_end), a gtt position for maximum length (max_end), a score 1 cutoff
	(co1), a score 2 cutoff (co2), base frequency dictionaries for R, L'' and
	L' (rbd, lbd1, and lbd2 respectively), an R' sequence to match to that has
	the empty string as a default value if there is no R' to match to, and a
	switch, final, that takes the value 0 if we are not searching for a final
	cassette (default) and 1 if we are. Searches for an acceptable attC in the
	sequence.  Returns the score data for the first acceptable attC it finds,
	unless final is 1, in which case it returns all potentially-acceptable
	attCs. If it doesn't find an acceptable attC, returns 0.'''


	''' Final Writing and unit tested'''
		
	gtt = 0
	while gtt != -1:
		'''This loop searches through gtts until either the score cutoff is met
		or the maximum attC length is reached.'''

		if min_end > max_end:
			#we have exceeded the maximum allowable attC len
			#without finding a good attC
			return 0

		if not final:
			#find the next gtt
			gtt,R2 = find_R(seqstr,'e',min_end,max_end)
		elif final:
			#find the next g
			gtt = seqstr.find('G',min_end,max_end)
						
		if gtt != -1:
			#Score the attC (fn not done)
			if not final:
				pot_attC = seqstr[aac:gtt+7]
			else:
				pot_attC = seqstr[aac:gtt+1]
			score_dict = score_whole_attC(pot_attC,rbd,final=final,R4=R4)
			score_dict['dir'] = 1

			#Check the scores against the cutoffs
			score_dict = check_scores(score_dict,co1,co2,bds,pot_attC)
			if 'score2' in score_dict: #this indicates it passed
				#add the start and end positions
				score_dict['start'] = aac
				score_dict['attC'] = pot_attC
				if not final:
					score_dict['end'] = gtt+7
				else:
					score_dict['end'] = gtt+1
				return score_dict

			else: #it failed
				#check this attC's reverse complement
				rev_attC = str(Seq(pot_attC).reverse_complement())
				score_dict = score_whole_attC(rev_attC,rbd,final=final,R4=R4,\
												rev=1)
				score_dict['dir'] = -1

				#Check the scores against the cutoffs
				score_dict = check_scores(score_dict,co1,co2,bds,rev_attC)
				if 'score2' in score_dict: #indicates it passed
					#add start and end positions
					score_dict['start'] = aac
					score_dict['attC'] = rev_attC
					if not final:
						score_dict['end'] = gtt+7
					else:
						score_dict['end'] = gtt+1
					return score_dict

				else: #the revcomp also failed
					#look for the next gtt
					if not final:
						min_end = gtt + 3
					else:
						min_end = gtt+1
	
	'''If the function gets out of the while loop, it has reached the end of the
	possible attC lengths without finding an acceptable attC.'''

	return 0

def attC_search(seqstr,min_start,max_start,co1,co2,aml,aMl,rbd,bds,R4='',final=0):
	'''A function that takes a sequence string, seqstr, an aac position for if
	the cassette has minimum length, min_start, an aac position for if the
	cassette has maximum length, max_start, cutoff values for scores 1 and 2
	(co1 and co2), minimum and maximum attC lengths (aml and aMl,
	respectively), base frequency dictionaries for R, L'', and L' (rbd,bds), a
	four-letter string that is the R' to match to (R1), and a switch, final.
	final will be 1 if we are searching for a final cassette and 0 otherwise.
	R1 will be an empty string if we are searching for the first cassette.
	Searches the entire range of possible cassette lengths for an acceptable
	attC. Returns the first good attC it finds. If no attC can be found within
	the allowable lengths, returns 0.'''


	''' Final Writing not unit tested (score2)'''

	aac = 0

	while aac != -1:
		'''This loop searches through aacs until either the score cutoff is met
		or the maximum cassette length is reached.'''

		if max_start > len(seqstr):
			max_start = len(seqstr) #dont' go beyond the end of the seq

		if min_start > max_start:
			#we have exceeded the maximum allowable cass len
			#without finding a good attC
			return 0
		
		#look for the next aac within the allowable cass lens
		aac,R1 = find_R(seqstr,'s',min_start,max_start)
					
		if aac != -1: #if aac was found
			#the place a gtt will be if the attC is min len
			min_end = aac + aml - 7
			#the place a gtt will be if the attC is max len
			max_end = aac + aMl - 7
			if max_end > len(seqstr):
				max_end = len(seqstr) #don't go beyond the end of the seq

			this_aac = \
			gtt_search(seqstr,aac,min_end,max_end,co1,co2,rbd,bds,R4=R4,final=final)

			if this_aac:
				'''An attC has been found for this cassette. Return it.'''

				return this_aac
			else:
				min_start = aac+7

		else:
			'''A suitable attC has not yet been found. Increment min_start and
			find the next AAC that could mark the end of this cassette. (i.e.
			continue with this loop).'''
			
			min_start = aac+7

	'''If the function gets out of the while loop, it has reached the end of the
	possible cassette lengths without finding an attC.'''

	return 0
							
def add_attC(attCs,this_attC):
	'''A function that takes a dictionary of found attCs (attCs) and a current
	attC dict (this_attC) and adds the current dict to the dict of found attCs.
	Returns the modified dict of found attCs, but also modifies it in place.'''

	
	'''Written and unit tested'''
	
	k = attCs.keys()
	if k == []:
		attCs[0] = this_attC
	else:
		num = max(k) + 1
		attCs[num] = this_attC

	return attCs

def all_intact_attCs(intd,cml,cMl,seqstr,co1,co2,aml,aMl,rbd,bds):
	'''A function that takes a dictionary corresponding to a single integron,
	intd, the minimum and maximum cassette lengths (cml and cMl), the sequence
	string, seqstr, the score 1 and 2 cutoffs (co1 and co2), the minimum and
	maximum attC lengths (aml and aMl), and the base frequency dictionaries
	(rbd,bds) and finds all the intact genomic attCs associated with that
	integron and adds them to the dictionary of attCs.  Returns a tuple whose
	first element is the modified intd and whose second and third elements are
	updated values for min_start and max_start that can be used to find the
	final attC. Also changes it in place.'''


	'''Written and not unit tested (score2)'''
	
	attCs = intd['attCs']

	#Get the first R' sequence
	R4 = intd['R1']
	start = intd['start']
	
	#Check bit, will be set to 0 when the last intact attC is found:
	this_attC = True 
	while this_attC:
		
		#the place aac will be if the cass is min len
		min_start = start + cml + 4
		#the place aac will be if the cass is max len
		max_start = start + cMl + 4
					
		#Search for an attC between min_start and max_start
		this_attC = attC_search(seqstr,min_start,max_start,co1,co2,aml,aMl,\
										rbd,bds,R4=R4)
		'''If no attC is found, the above will return 0 and we will
		exit the while loop.'''
		
		if this_attC: #if an attC is found
			#increment start and update R4
			start = this_attC['end']
			R4 = this_attC['attC'][-4:]
			
			#add the current attC to the attCs
			attCs = add_attC(attCs,this_attC)

	'''The above loop will end when the last attC has been found. Once
	it is found we must search for the final attC.'''

	return (intd,min_start,max_start)

def find_starts(nattI,seqstr,attIfnm,attIseq,co1,co2,aml,aMl,rbd,bds):
	'''A function that takes a binary switch, nattI, that is 0 if we are using
	attI and 1 if we are not, a sequence string (seqstr), the name of the file
	where the default attIs are stored (attIfnm), a user-provided sequence for
	a non-default attI sequence (attIseq), score 1 and 2 cutoffs (co1,co2), min
	and max attC lengths (aml,aMl), frequency dictionaries for R, (rbd,bds) and
	returns a tuple whose first element is a switch, done, that is False if
	there are potentially more attCs to find in this direction and True if
	there are not, and whose second element is a list containing all the
	integron dictionaries for this direction.'''


	'''Written and unit tested.'''

	#create a check bit for identifying when we're done searching for attCs
	done = False

	#Check if we're using attI
	if not nattI:
		####################################################################
		#yes using attI
		####################################################################
		print 'Searching for attI(s)...'

		'''Find all the attIs associated with this direction & store them
		in integron dictionaries in a list.'''

		this_dir = find_all_attIs(seqstr,attIfnm,aml,attIseq=attIseq)
		if this_dir == []: #if there were no attIs
			done = True

	else:
		####################################################################
		#not using attI
		####################################################################
		'''Since we're not using attI we find the first attC after the start
		point and use it as our only start in this sequence.'''

		print 'Searching for first attC...'

		mco1 = co1 -4 +2
		#create list to hold the int dicts
		this_dir = []
		
		#Get an integron dictionary of the first attC
		first_attC = find_first_attC(seqstr,mco1,co2,aml,aMl,rbd,bds)
	
		if first_attC:
			#put in in the list for this direction
			this_dir.append(first_attC)
			
		else:
			#no first attC was found. Search for an incomplete
			#co1 mod'ed for 1st cass & incomplete: co1-(4-mean)-mean
			mco1 = co1 - 4
			first_attC = find_first_attC(seqstr,mco1,co2,aml,aMl,rbd,bds,\
											final=1)
				
			if first_attC:
				#Only 1 attC in this direction.
				#Put it in the list
				this_dir.append(first_attC)
				
				#Move on to file writing
				done = True
			else:
				#No attCs in this direction
				#Move on to file writing
				done = True
	return (done,this_dir)

def find_all_attCs(this_dir,cml,cMl,seqstr,co1,co2,aml,aMl,rbd,bds):
	'''A function that takes the list of integron dicts in this direction
	(this_dir), min and max cassette lengths (cml,cMl), the sequence string
	(seqstr), score 1 and 2 cutoffs (co1,co2), min and max attC lengths
	(aml,aMl), base frequency dictionaries (rbd,bds), and iterates through
	this_dir to find all the attCs associated with each integron in a given
	direction. Returns the this_dir list, which has been modified so that all
	the integron dictionaries in it now contain all the attCs associated with
	that integron.  Also modifies this_dir in place.'''


	'''Written and not unit tested (score2)'''

	for i in range(len(this_dir)):
		print 'Searching for attCs...'
		intd = this_dir[i]
		
		#Create a dictionary to store the attCs in
		intd['attCs'] = {}
		attCs = intd['attCs']
		
		#Find all the intact attCs associated with this integron and get start 
		#values for finding the final attC.
		intd,min_start,max_start = all_intact_attCs(intd,cml,cMl,seqstr,\
														co1,co2,aml,aMl\
														,rbd,bds)
	return this_dir

def all_but_write(seq,nattI,attIfnm,attIseq,co1,co2,aml,aMl,rbd,bds,\
					cml,cMl,sp,ep):
	'''A function that takes A BUNCH OF STUFF and returns a list of integron
	annotation dictionaries to be written to a file.'''

	if ep == 0:
			ep = len(seq)
	seqstr = seq.upper()
	seqstr = seqstr[sp:ep]
	
	#Find all the starting integrons in this dir
	done,this_dir = find_starts(nattI,seqstr,attIfnm,attIseq,co1,co2,\
									aml,aMl,rbd,bds)
	
	'''Staying within a given sequence direction, we now iterate through all
	the potential integrons identified above and find all the attCs.'''
		
	if not done: #if it looks like there are more attCs
		this_dir = find_all_attCs(this_dir,cml,cMl,seqstr,co1,co2,aml\
										,aMl,rbd,bds)

	return this_dir,seqstr


def search_all_directions(seqs,nattI,attIfnm,attIseq,co1,co2,aml,aMl,rbd,bds,\
							cml,cMl,sp,ep,outfnm,inner=0):
	'''A function that takes a list of sequences (seqs), the nattI switch
	(nattI), the name of the default attI file (attIfnm), a user-provided attI
	sequences (attIseq), score 1 and 2 cutoffs (co1,co2), max and min attC
	lengths (aml,aMl), base freq dicts (rbd,bds), min and max cassette lengths
	(cml,cMl), start and end points for the sequence to be searched (sp,ep),
	and the output file name (outfnm), A BUNCH OF STUFF. Iterates through all
	the sequences in seqs and finds all the integrons and attCs associated with
	them, then writes a csv file for each one. Returns None.'''


	'''Written and unit tested.'''

	############################################################################
	#Iterate through the sequences
	############################################################################
	for seqnum in range(len(seqs)): #iterate through the sequence directions
		print 'Analyzing direction %i...' % (seqnum+1,)

		#get a plain, upper case string of the sequence
		seq = seqs[seqnum]

		#Find all attachment sites
		this_dir,seqstr = all_but_write(seq,nattI,attIfnm,attIseq,co1,co2,aml,\
									aMl,rbd,bds,cml,cMl,sp,ep)

		#Once searching for attC is over, write the files
		#Write the files
		print 'Writing attC csv file for direction %i...' % (seqnum+1,)
		write_files(this_dir,outfnm,seqnum,sp)

		########################################################################
		#Find the cassettes and write the files
		########################################################################
		print 'Writing cassette fasta file for direction %i...' % (seqnum+1,)
		get_and_write_all_cass(this_dir,seqstr,sp,outfnm,seqnum,inner)

	return
				
################################################################################
#attC Scoring
################################################################################
def seq_match(s1,s2):
	'''A function that takes two dna sequences (one of which has been revcomped 
	if necessary), s1 and s1, and assigns points based on the number of matches 
	between the two. Returns the integer number of points.'''


	''' Written and unit tested.'''
 
	#check the inputs
	if len(s1) != len(s2):
		msg = 'The two sequences are not the same length.'
		raise ValueError(msg)
 
	#match
	points = 0
	for i in range(len(s1)):
		if s1[i] == s2[i]:
			points += 1
	return points
 
def base_pat_score(s,bd):
	'''A function that takes a dna sequence, s, and a base
	frequency dictionary, bd, and calculates a score out of 4 for the s. The
	score is calculated by giving a score out of 1 for each base in the
	following way: the score of a given base at a given position is the
	frequency of that base at that position divided by the frequency of the most
	common base at that position. Returns the score as a float.'''


	''' Final Writing and unit tested.'''
	
	#check inputs
	msg = ''
	for i in s:
		if (i != 'A') and (i != 'T') and (i != 'G') and (i != 'C'):
			msg += '\nThe sequence can only have A, T, G, and C.\n'
			break

	k = bd.keys()
	k.sort()
	if 'tot' not in k:
		m = '\nThe keys of the base frequency dictionary must be the '+\
				'numbers from 0 to n-1 where n is the length of sequence to '+\
				'be tested, plus \'tot\'.\n'
		raise ValueError(m)

	else:
		k.remove('tot')
		if k != range(len(s)):
			msg += '\nThe keys of the base frequency dictionary must be the '+\
					'numbers from 0 to n-1 where n is the length of sequence to '+\
					'be tested, plus \'tot\'.\n'

	for i in k:
		j = bd[i].keys()
		j.sort()
		if j != ['A','C','G','T']:
			msg += '\nThe keys of the base frequency subdictionary at '+\
					'position %i' %i +' are not A, C, G, T.\n'

		freqs = bd[i].values()
		if sum(freqs) != bd['tot']:
			msg += '\nThe base frequencies at position %i' %i +' don\'t sum '+\
					'to 1.\n'
	if msg:
		raise ValueError(msg)

	#score
	points = 0

	for i in range(len(s)):
		vals = bd[i].values() #get all frequencies at this position
		max_val = max(vals) #get the maximum values
		
		points += bd[i][s[i]]/float(max_val)

	return points

def trip_not_same(trip):
	'''A function that takes a triplet string, trip. Returns 0 if the three
	letters in that triplet are all the same, 1 if not.'''

	
	''' Unit tested.'''

	if len(trip) != 3:
		msg = 'This triplet is not 3 letters long.'
		raise ValueError(msg)

	if (trip[0]==trip[1]==trip[2]):
		return 0
	else:
		return 1

def score_triplet_pair(t1,t2):
	'''A function that takes two potentially-matching triplets from L sites and
	returns their triplet score.'''


	''' Unit tested. '''

	
	#each triplet pair starts with a score of 0
	score = 0

	#add one point for each match between the pairs
	score += seq_match(t1,t2) #written and tested

	#add one point for each of the two triplets that isn't the same base 3 times
	score += trip_not_same(t1) #written and tested
	score += trip_not_same(t2)

	return score

def score_L(L1,L2):
	'''A function that takes two potential L sites (L1 and L2) and finds their
	match score'''


	''' Unit tested. '''


	#L1 should have an extrahelical base and therefore be longer than L2 by 1
	if len(L1) != len(L2)+1:
		msg = 'L\'\' should be one bp longer than L\'.'
		raise ValueError(msg)
	
	#score the triplets
	trip1 = score_triplet_pair(L1[:3],L2[:3]) #written and tested
	trip2 = score_triplet_pair(L1[-3:],L2[-3:])

	score = trip1*trip2

	return score


def score_Ls(att):
	'''A function that takes a dna string from an attC that starts immediately
	after R'' and ends immediately before R' (seq) and checks all 4 possible L
	scores.  Returns a score dictionary for the best-scoring L. The dictionary
	will have keys L1,L2,L1-pos,L2-pos,L-score,ehb.'''


	''' written and unit tested. '''
	
	#A dictionary to store all the scores in
	all_scores = {}

	#Take the complement of the sequence
	revatt = str(Seq(att).reverse_complement())

	for e in range(5,7): #iterates through 5 and 6 for L2
		for s in range(5,7): #iterates through 5 and 6 for L1
			pos = '%i-%i' % (s,e) #the key for all_scores for this combo
			all_scores[pos] = {}
			
			#add the L sequences to the dictionary
			L1 = att[s:s+7]
			all_scores[pos]['L1'] = L1
			L2 = revatt[e:e+6]
			all_scores[pos]['L2'] = str(Seq(L2).reverse_complement())
			
			all_scores[pos]['L-score'] = score_L(L1,L2) #written and tested

	#Choose the best score (if there is more than one max, takes at random)
	scores = []
	for k in all_scores:
		scores.append(all_scores[k]['L-score'])
	max_score = max(scores)
	for pos in all_scores:
		if all_scores[pos]['L-score'] == max_score:
			best_score = all_scores[pos]
			best_score['L1-pos'] = int(pos[0])
			best_score['L2-pos'] = int(pos[2])
			break

	#Double the L structure score to increase its weighting
	best_score['L-score'] = best_score['L-score'] * 2
	L1 = best_score['L1']
	L2 = str(Seq(best_score['L2']).reverse_complement())

	#check the extrahelical base
	if L1[3] != L2[3]:
		best_score['L-score'] += 1
		best_score['ehb'] = 1
	else:
		best_score['ehb'] = 0

	return best_score


def score_whole_attC(attC,rbd,final=0,R4='',rev=0):
	'''A function that takes an attC sequence (attC), a switch (final) that
	will be 0 if this is not a final attC and 1 if it is, a switch (rev) that
	takes 0 (default) if the attC being tested is not reversed from the overall
	sequence, and 1 if it is, a four-character string, R4 that is the R' of the
	previous attI/attC to match with, and a base frequency dictionary for R
	(rbd).  If this is the first attC with no previous R' then R4 will be an
	empty string.  Scores the attC.  Returns a dictionary that has all the
	features of the attC so they can be properly annotated.'''


	''' Final Writing and unit tested '''

	#Check the R sites
	##Check R''
	if rev: #this will be the last four bases, and is already revcomped
		r1r = attC[-4:]
	else:
		r1 = attC[:4]
		r1r = str(Seq(r1).reverse_complement())
	if R4:
		#match R'' to the old R'
		r1_points = seq_match(r1r,R4) #written and tested
	else:
		#use R base pattern matching
		r1_points = base_pat_score(r1r,rbd) #written, not tested

	##Check R'
	if not final: #there will be an R'
		if rev: #this will be the first four bases and needs to be revcomped
			r2r = attC[:4]
			r2 = str(Seq(r2r).reverse_complement())
		else:
			r2 = attC[-4:]
		r2_points = base_pat_score(r2,rbd)

	else: #there is no R'
		r2_points = 0
	

	#Check the L sites
	if not final:
		Ls = attC[7:-7]
	else:
		if not rev:
			Ls = attC[7:-1]
		else:
			Ls = attC[1:-7]

	score_d = score_Ls(Ls) #written, not tested
	score_d['score'] = score_d['L-score']+r1_points+r2_points
	score_d['r1_score'] = r1_points
	if final:
		score_d['r2_score'] = 'final'
	else:
		score_d['r2_score'] = r2_points

	return score_d


def old_score2(score_dict,lbd1,lbd2):
	'''A function that takes an attC score dictionary, score_dict, an L1 base
	frequency dictionary that has the base frequencies for L'' (lbd1) and an L2
	base frequency dictionary that has the base frequencies for L' (lbd2) and
	scores the base frequencies in the L.  Returns a score.'''


	'''Final Writing and unit tested.'''
	
	score = 0

	#score L1
	L1 = score_dict['L1']
	score += base_pat_score(L1,lbd1)

	#score L2
	L2 = score_dict['L2']
	score += base_pat_score(L2,lbd2)

	return score

def score2(score_dict,attC,bds,final=0):
	'''A function that takes an attC score dictionary, score_dict, the attC
	sequence, attC, and a dictionary of base dicts, bds. Returns a score based
	on nucleotide frequencies at a number of positions.'''

	
	'''Written and unit tested'''

	score = 0
	
	#score Ls
	L1 = score_dict['L1']
	score += base_pat_score(L1,bds['L1'])
	L2 = score_dict['L2']
	score += base_pat_score(L2,bds['L2'])

	#score first gap and first internal
	if score_dict['L1-pos'] == 5:
		G1 = attC[7:12]
		I1 = attC[19:25]
		score += base_pat_score(G1,bds['G1-5'])
		score += base_pat_score(I1,bds['I1'])
	elif score_dict['L1-pos'] == 6:
		G1 = attC[7:13]
		I1 = attC[20:26]
		#normalize this score so it doesn't have an advantage
		g1s = base_pat_score(G1,bds['G1-6'])
		g1s = (g1s/6)*5
		score += g1s
		score += base_pat_score(I1,bds['I1'])
	else:
		msg = '\nL1 must be either 5 or 6 bases from R1'
		raise ValueError(msg)
	
	#score second gap and internal
	if not final:
		if score_dict['L2-pos'] == 5:
			G2 = attC[-12:-7]
			I2 = attC[-24:-18]
			score += base_pat_score(G2,bds['G2-5'])
			score += base_pat_score(I2,bds['I2'])
		elif score_dict['L2-pos'] == 6:
			G2 = attC[-13:-7]
			I2 = attC[-25:-19]
			#normalize this score so it doesn't have an advantage
			g2s = base_pat_score(G2,bds['G2-6'])
			g2s = (g2s/6)*5
			score += g2s
			score += base_pat_score(I2,bds['I2'])
		else:
			msg = '\nL2 must be either 5 or 6 bases from R2'
			raise ValueError(msg)
	else:
		if score_dict['L2-pos'] == 5:
			G2 = attC[-6:-1]
			I2 = attC[-18:-12]
			score += base_pat_score(G2,bds['G2-5'])
			score += base_pat_score(I2,bds['I2'])
		elif score_dict['L2-pos'] == 6:
			G2 = attC[-7:-1]
			I2 = attC[-19:-13]
			#normalize this score so it doesn't have an advantage
			g2s = base_pat_score(G2,bds['G2-6'])
			g2s = (g2s/6)*5
			score += g2s
			score += base_pat_score(I2,bds['I2'])
		else:
			msg = '\nL2 must be either 5 or 6 bases from R2'
			raise ValueError(msg)
		

	#score Rs
	R1 = attC[:4]
	score += base_pat_score(R1,bds['R1'])
	if not final:
		R2 = attC[-4:]
		score += base_pat_score(R2,bds['R2'])

	return score
	
################################################################################
#attC Mode
################################################################################
def check_seq_pats(seq):
	'''A function that takes a sequence string (seq) and checks that it follows
	the attC pattern (i.e. that it starts with NNNNAAC, ends with GTTNNNN, and
	is at least 31 bases long. Returns 1 if the conditions are met, 0 if not.'''
	
	'''Written and unit tested.'''

	if seq[-7:-4].upper() != 'GTT':
		return 0
	if seq[4:7].upper() != 'AAC':
		return 0
	if len(seq) < 31:
		return 0

	return 1

def score_all_attCs(in_file,rbd,bds):
	'''A function that takes the name of a fasta file of attCs (in_file) and
	the base dictionaries (rbd,bds) and returns a list of tuples whose first
	element is an annotated attC and whose second element is either an attC
	score dictionary or, if the sequence didn't appear to be an attC, 0.'''


	'''Written and unit tested.'''

	#Open the fasta file
	recs = SeqIO.parse(in_file,'fasta')
	
	#Create a list to store the results
	score_lst = []

	#Iterate through the records
	for rec in recs:
		seq = str(rec.seq)
		#check the sequence
		ok_seq = check_seq_pats(seq)

		#if the sequence looks like an attC, score it
		if ok_seq:
			score_dict = score_whole_attC(seq,rbd)
			
			#Add score 2
			score_dict['score2'] = score2(score_dict,seq,bds)
			
			#Add the sequence to the dictionary
			score_dict['attC'] = seq
			
			#Add this dictionary to the list
			score_lst.append((rec.id,score_dict))

		else: 
			score_lst.append((rec.id,0))

	return score_lst
	
################################################################################
#Writing Output
################################################################################
def write_attC_files(score_lst,fnm):
	'''A function that takes a list of attC score dictionary tuples and a file
	name and writes a csv annotating each attC.'''


	'''Written and unit tested'''

	#open the file for writing
	fnm = fnm+'_attC_scores.csv'
	outf = open(fnm,'w')

	#Create the header line
	L = 'ID,R\'\',R\'\'-score,5/6,L\'\',ehb,L\''+\
				',L-score,5/6,R\',R\'-score,final_score,score2\n'
	outf.write(L)
	#Write all the scores
	for stpl in score_lst:	
		#id
		L = str(stpl[0])+','
		#attC
		attCd = stpl[1]
		if attCd: #if this passed the attC pattern test
			L += str(attCd['attC'][:4])+','
			L += str(attCd['r1_score'])+','
			#for 5/6:
			pos = attCd['L1-pos']
			L += str(attCd['attC'][7:7+pos])+','
			L += str(attCd['L1'])+','
			if attCd['ehb']:
				L += str(attCd['L1'][3])+','
			else:
				L += 'X,'
			L += str(attCd['L2'])+','
			L += str(attCd['L-score'])+','
			#for next 5/6
			pos = attCd['L2-pos']
			if attCd['r2_score'] == 'final':
				L += str(attCd['attC'][-1-pos:-1])+','
			else:
				L += str(attCd['attC'][-7-pos:-7])+','
			if attCd['r2_score'] == 'final':
				L+= 'final,'
			else:
				L += str(attCd['attC'][-4:])+','
			L += str(attCd['r2_score'])+','
			L += str(attCd['score'])+','
			L += str(attCd['score2']) +'\n'
		else: #if this didn't pass the pattern test
			L += 'This doesn\'t appear to be an attC.\n'
		outf.write(L)
	
	outf.close()

	return


def write_files(int_lst,fnm,seqnum,sp):
	'''A function that takes a list of integron dictionaries (int_lst), a file
	name (fnm), a number indicating whether this is the first (0) or second (1)
	direction (seqnum), and the start position of the sequence (sp) and writes
	the attC annotations to a file.'''

	
	'''Written and unit tested.'''

	#Create a name for the file based on the sequence number
	fnm = fnm+str(seqnum+1)+'.csv'

	#Open the file for writing
	outf = open(fnm,'w')

	#Iterate through the integrons
	for intnum in range(len(int_lst)):
		int_dict = int_lst[intnum]
		L = 'Integron %i:\n' % (intnum+1,)
		outf.write(L)

		#Check if there's an attI
		if 'attISeq' in int_dict:
			#Write the attI lines
			L = 'attI:\n'
			L += 'attIStart,attIEnd,R\',attISeq\n'
			L += str(int_dict['attIStart']+sp+1)+','
			L += str(int_dict['attIEnd']+sp)+','
			L += str(int_dict['R1'])+','
			L += str(int_dict['attISeq'])+'\n'
			outf.write(L)
			#Create the header for the attCs
			L = 'attCs:\n'
			L += 'StartPos,EndPos,dir,R\'\',R\'\'-score,5/6,L\'\',ehb,L\''+\
					',L-score,5/6,R\',R\'-score,final_score,attCSeq\n'
			outf.write(L)

		elif 'FirstattC' in int_dict:
			#Write the attC headers
			L = 'No attI. attCs:\n'
			L += 'StartPos,EndPos,dir,R\'\',R\'\'-score,5/6,L\'\',ehb,L\''+\
					',L-score,5/6,R\',R\'-score,final_score,attCSeq\n'
			outf.write(L)
			#Write the line for the first attC
			if int_dict['Firstdir'] == -1:
				int_dict['FirstattC'] = \
						str(Seq(int_dict['FirstattC']).reverse_complement())
			L = ''
			L += str(int_dict['Firststart']+sp+1)+','
			L += str(int_dict['Firstend']+sp)+','
			L += str(int_dict['Firstdir'])+','
			L += str(int_dict['FirstattC'][:4])+','
			L += str(int_dict['Firstr1_score'])+','
			#for 5/6:
			pos = int_dict['FirstL1-pos']
			L += str(int_dict['FirstattC'][7:7+pos])+','
			L += str(int_dict['FirstL1'])+','
			if int_dict['Firstehb']:
				L += str(int_dict['FirstL1'][3])+','
			else:
				L += 'X,'
			L += str(int_dict['FirstL2'])+','
			L += str(int_dict['FirstL-score'])+','
			#for next 5/6
			pos = int_dict['FirstL2-pos']
			if int_dict['Firstr2_score'] == 'final':
				L += str(int_dict['FirstattC'][-1-pos:-1])+','
			else:
				L += str(int_dict['FirstattC'][-7-pos:-7])+','
			if int_dict['Firstr2_score'] == 'final':
				L += 'final,'
			else:
				L += str(int_dict['FirstattC'][-4:])+','
			L += str(int_dict['Firstr2_score'])+','
			L += str(int_dict['Firstscore'])+','
			L += str(int_dict['FirstattC'])+'\n'
			outf.write(L)

		#Write the rest of the attCs
		attCs = int_dict['attCs']
		if len(attCs) > 0:
			keys = attCs.keys()
			keys.sort()
			for k in keys:
				d = attCs[k]
				if d['dir'] == -1:
					d['attC'] = str(Seq(d['attC']).reverse_complement())
				L = ''
				L += str(d['start']+sp+1)+','
				L += str(d['end']+sp)+','
				L += str(d['dir'])+','
				L += str(d['attC'][:4])+','
				L += str(d['r1_score'])+','
				#for 5/6:
				pos = d['L1-pos']
				L += str(d['attC'][7:7+pos])+','
				L += str(d['L1'])+','
				if d['ehb']:
					L += str(d['L1'][3])+','
				else:
					L += 'X,'
				L += str(d['L2'])+','
				L += str(d['L-score'])+','
				#for next 5/6
				pos = d['L2-pos']
				if d['r2_score'] == 'final':
					L += str(d['attC'][-1-pos:-1])+','
				else:
					L += str(d['attC'][-7-pos:-7])+','
				if d['r2_score'] == 'final':
					L+= 'final,'
				else:
					L += str(d['attC'][-4:])+','
				L += str(d['r2_score'])+','
				L += str(d['score'])+','
				L += str(d['attC'])+'\n'
				outf.write(L)
		else:
			L = 'No (more) attCs.\n'
			outf.write(L)

	if int_lst == []:
		L = 'No integrons found in this direction.\n'
		outf.write(L)

	outf.close()
	return

def write_cass_files(outfnm,seqnum,intnum,cass_lst,sp):
	'''A function that takes the output file name base (outfnm), the sequence
	number (seqnum), the integron number (intnum), a list of cassettes
	(cass_lst), and the start position of the sequence, sp, and writes a fasta
	file for all the cassettes in that integron in order.'''


	'''under construction and not unit tested'''

	#Open the file for writing
	outfnm = outfnm+'_'+str(seqnum+1)+'_integron_'+str(intnum+1)+'_cassettes.fna'
	outf = open(outfnm,'w')

	#Iterate through the cassettes
	for i in range(len(cass_lst)):
		tpl = cass_lst[i]
		
		#Create the sequence record object
		seq = Seq(tpl[2])
		rec = SeqRecord(seq)
		rec.id = str(i)
		rec.description = str(tpl[0]+sp+1)+' - '+str(tpl[1]+sp)

		#Write to the file
		L = rec.format('fasta')
		outf.write(L)

	outf.close()
	return

###############################Cassette Functions###############################

def get_whole_cass(attC1,attC2,seqstr,final=0):
	'''A function that takes the end points of two attCs (or one attI and one
	attC) in order (attC1, attC2), a sequence string (seqstr), and a switch that
	is 1 if this is the final cassette and 0 if it is not and returns a tuple
	whose first two elements are the start and end positions of the cassette and
	whose third element is the cassette sequence.'''

	#Get the cassette boundaries
	cass_start = attC1 - 6
	if not final:
		cass_end = attC2 - 6
	else:
		cass_end = attC2

	#Get the sequence
	cass_seq = seqstr[cass_start:cass_end]

	return (cass_start,cass_end,cass_seq)

def get_inner_cass(attC1,attC2,seqstr):
	'''A function that takes the end point of an attI or attC (attC1) and the
	start point of the next attC (attC2), and a sequence string (seqstr).
	Returns a tuple whose first two elements are the start and end points of the
	inner section of the cassette bounded by these attCs (inner section means
	the cassette excluding the attCs) and whose final element is the sequence of
	that inner section.'''
	
	cass = seqstr[attC1:attC2]
	
	return (attC1,attC2,cass)

def get_cass_from_csv(csv):
	'''A function that takes a csv of an annotated integron and the associated
	sequence (seqstr) and extracts all but the last cassette.'''
	
	#Open the file
	pass

def get_all_cass(intd,seqstr,inner=0):
	'''A function that takes an integron dictionary that has all the annotated
	attI/attCs of a given integron in it (intd), the sequence this was taken from
	(seqstr), and a switch that is 1 if we are looking for inner cassettes and 0
	otherwise (inner) and returns a list of tuples in order where each tuple
	contains the start position, end position, and sequence of a cassette.'''


	'''Written and unit tested'''
	
	attCs = intd['attCs']

	#Create the list
	cass_lst = []
	
	#Check if there are any cassettes
	if len(attCs) == 0: 
		#there are no cassettes 
		return []
	
	#Get the first cassette
	#Check if there is more than one cassette (i.e. if the first is also last):
	if len(attCs) == 1:
		#there is only one cassette
		#inner or outer?
		if inner:
			#get cassette boundaries
			s = intd['start']
			e = attCs[0]['start']
			cass = get_inner_cass(s,e,seqstr)
		else:
			s = intd['start']
			e = attCs[0]['end']
			cass = get_whole_cass(s,e,seqstr)
		cass_lst.append(cass)
		return cass_lst
	else:
		#there is more than one cassette
		#Get the first cassette
		#inner or outer?
		if inner:
			#get cassette boundaries
			s = intd['start']
			e = attCs[0]['start']
			cass = get_inner_cass(s,e,seqstr)
		else:
			s = intd['start']
			e = attCs[0]['end']
			cass = get_whole_cass(s,e,seqstr)
		cass_lst.append(cass)
		#iterate through all the attCs
		for i in range(len(attCs)-1):
			#inner or outer
			if inner:
				s = attCs[i]['end']
				e = attCs[i+1]['start']
				cass = get_inner_cass(s,e,seqstr)
			else:
				s = attCs[i]['end']
				e = attCs[i+1]['end']
				cass = get_whole_cass(s,e,seqstr)

			cass_lst.append(cass)

	return cass_lst

def get_and_write_all_cass(this_dir,seqstr,sp,outfnm,seqnum,inner=0):
	'''A function that takes a list of integrons in this direction (this_dir),
	a string of the search sequence (seqstr), the start position (sp), the base
	name of the output files, outfnm, the sequence number, seqnum, and a
	switch, inner, that is 1 if the cassettes should exclude the attC and 0 if
	not. Writes annotated cassettes to fasta files.'''
	

	'''Written but not unit tested.'''

	#Iterate through the list of integrons in this direction
	for i in range(len(this_dir)):
		int_dict = this_dir[i]
		
		#Get a list of all the cassettes
		cass_lst = get_all_cass(int_dict,seqstr,inner)

		#Write the cassettes to a file
		write_cass_files(outfnm,seqnum,i,cass_lst,sp)

	return
