#!/usr/bin/env python

from optparse import OptionParser
import ISAAC_funcs as apf
import pickle

#function definitions might go here once development is done

################################################################################
#Define Variables
################################################################################
attIfnm = 'attIs.txt'
rbd = pickle.load(open('rbd.pickle'))
#lbd1 = pickle.load(open('lbd1.pickle'))
#lbd2 = pickle.load(open('lbd2.pickle'))
bds = pickle.load(open('all_bfds.pickle'))

################################################################################
#Start
################################################################################
if __name__ == '__main__':

  ############################################################################
	#Parse all the options with optparse
	############################################################################
	print 'Parsing inputs...'
	parser = OptionParser()

	options,args = apf.do_args(parser)


	############################################################################
	#Check all the inputs and set variables
	############################################################################
	print 'Checking input values...'
	msg = ''
	#an input fasta file must be provided
	if not options.in_file: 
		msg += '\nA fasta file must be provided using the -in or --in-file'+\
					' flag.\n'
	else:
		fasta = options.in_file

	if not options.out: #an output file name must be provided
		msg += '\nA name for the output file must be provided using the'+\
					'-out or --out-file flag.\n'

	else:
		outfnm = options.out

	inner = options.ci #this doesn't get checked

	if options.attISeq: #attI sequences must end in GTT
		if options.attISeq[-3:].upper != 'GTT':
			msg += '\nThe attI sequence that you provide using the -attISeq'+\
					' flag must end in GTT. You may choose not to use '+\
					'this flag, in which case the default will be used.\n'
		else:
			attIseq = options.attISeq
	else:
		attIseq = False

	nattI = options.nattI #this doesn't get checked

	if (options.co1 < 0) or (options.co1) > 59:
		msg += '\nThe score 1 cutoff provided using the flag -co1 or --score1'+\
				'-cutoff must be an integer between 0 and 59. You may choose'+\
				' not to use this flag, in which case the default value of 75'+\
				' will be used.'
	else:
		co1 = options.co1

	if (options.co2 < 0) or (options.co2) > 50:
		msg += '\nThe score 2 cutoff provided using the flag -co2 or --score2'+\
				'-cutoff must be an integer between 0 and 50. You may choose'+\
				' not to use this flag, in which case the default value of 75'+\
				' will be used.'
	else:
		co2 = options.co2

	attCo = options.attCo #doesnt get checked

	if options.aml < 20: #attC minimum length can't be less than 20
		msg += '\nYour minimum attC length provided using the flag --aml or '+\
				'--attCminLen must be at least 20. You may choose not to '+\
				'use this flag, in which case the default will be used.\n'
	else:
		aml = options.aml

	if options.aMl < aml: #attC max length can't be less than min
		msg += '\nThe maximum attC length may not be shorter than the '+\
				'minimum.\n'
	else:
		aMl = options.aMl

	if (options.cml > options.cMl) and (options.cMl != 0): 
		msg += '\nYour minimum attC length must be less than the max if the '+\
				'max is not set to 0 (no limit)\n'
	else:
		cml = options.cml 
		cMl = options.cMl 

	#skipping sp and ep because they require information about the sequence

	if (options.direction != 1) and (options.direction != 0) and \
			(options.direction != -1):
		msg += '\nThe provided sequence direction (using the -dir or '+\
				'--direction flag) must be 1 (for forward), -1 (for '+\
				'reverse), or 0 (for both). You can choose not to use '+\
				'this flag, in which case the default 0 will be used.\n'
	else:
		drn = options.direction

	if msg: #If anything failed, raise this exception
		raise ValueError(msg)
	
	############################################################################
	#Check if running in attC mode or regular mode
	############################################################################
	if not attCo: #not running in attC mode

		########################################################################
		#Open the fasta file
		########################################################################
		print 'Opening fasta file...'
		rec = apf.open_fasta(fasta,aml)
		
		########################################################################
		#Check the input values that couldn't be checked without the file 
		########################################################################
		print 'Checking fasta-related inputs...'
		sp = options.sp
		ep = options.ep
		
		msg = ''
		if not ep:
			ep = len(rec)
			end_name = 'the sequence length'
		elif ep > len(rec):
			ep = len(rec)
			end_name = 'the sequence length'
			print 'The specified end point is greater than the sequence lengt'+\
					'h. Using the end of the sequence.'
			
		else:
			end_name = 'your end position'
		
		if sp > ep - aml: #must start at least 1 attC length before end
			msg += 'The start position provided using the -sp or '+\
			'--start-position flag must be at least one attC length less '+\
			'than %s.\n' % end_name
		if msg:
			raise ValueError(msg)
		
		#adjust the indexing of sp to 0-indexed now that it's been checked
		sp = sp-1
		if sp < 0:
			sp = 0
	
		########################################################################
		#generate a list of the sequences to use based on which direction we're
		#checking
		########################################################################
		print 'Getting all requested directions...'
		seqs = apf.get_seqs(drn,rec)
		
		########################################################################
		#Find the features and write the files
		########################################################################
		print 'Executing...'
		apf.search_all_directions(seqs,nattI,attIfnm,attIseq,co1,co2,aml,aMl,\
									rbd,bds,cml,cMl,sp,ep,outfnm,inner)
		

	############################################################################
	#If running in attC mode
	############################################################################
	else:
		#Score all the attCs in the fasta file
		score_lst = apf.score_all_attCs(fasta,rbd,bds)

		#Write the output file
		apf.write_attC_files(score_lst,outfnm)

################################################################################
#End
################################################################################
