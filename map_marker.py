### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###


__usage__ = """  python3 map_marker.py
							--in <GENETIC_MARKER_TABLE>
							--ref <REFERENCE_SEQUENCE_FILE>
							--target <TARGET_SEQUENCE_FILE>
							--out <OUTPUT_FOLDER>
					"""

import os, sys, subprocess

# --- end of imports --- #


def load_marker_positions( input_table ):
	"""! @brief load marker positions """
	
	markers = []
	with open( input_table, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			#markers.append( { 'id': parts[0], 'chr': parts[1], 'pos': int( ( int( parts[9] ) + int( parts[10] ) ) / 2.0 ) } )
			if len( parts ) > 2:
				markers.append( { 'id': parts[0], 'chr': parts[1], 'pos': int( parts[2] ) } )
			line = f.readline()
	return markers


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if "_" in header:
			header = header.split("_")[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if "_" in header:
						header = header.split("_")[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					best_hits[ parts[0] ] = { 	'score': float( parts[-1] ),
																'chr': parts[1],
																'pos': int( ( int( parts[8] ) + int( parts[9] ) ) / 2.0 ),
																'id': parts[0]
															}
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ),
																	'chr': parts[1],
																	'pos': int( ( int( parts[8] ) + int( parts[9] ) ) / 2.0 ),
																	'id': parts[0]
																	 } } )
			line = f.readline()
	return best_hits


def main( arguments ):
	
	input_table = arguments[ arguments.index('--in')+1 ]
	ref_seq_file = arguments[ arguments.index('--ref')+1 ]
	target_seq_file = arguments[ arguments.index('--target')+1 ]

	output_folder = arguments[ arguments.index('--out')+1 ]
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	marker_positions = load_marker_positions( input_table )
	assembly = load_sequences( ref_seq_file )
	marker_blast_file = output_folder + "marker.fasta"
	with open( marker_blast_file, "w" ) as out:
		for marker in marker_positions:
			try:
				out.write( '>' + marker['id'] + "\n" + assembly[ marker['chr'] ][ marker['pos']-100: marker['pos']+100 ] + "\n" )
			except KeyError:
				pass


	blastdb = output_folder + "blastn_db"
	blast_result_file = output_folder + "blast_hits.txt"

	if not os.path.isfile( blast_result_file ):
		p = subprocess.Popen( args= "makeblastdb -in " + target_seq_file + " -out " + blastdb + " -dbtype nucl", shell=True )
		p.communicate()

		p = subprocess.Popen( args= "blastn -query " + marker_blast_file + " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads 10", shell=True )
		p.communicate()


	best_hits = load_best_blast_hit( blast_result_file )

	final_output_file = output_folder + "final_result_file.txt"
	with open( final_output_file, "w" ) as out:
		for hit in list( best_hits.values() ):
			out.write( "\t".join( list( map( str, [ hit['id'], hit['chr'], hit['pos'] ] ) ) ) + "\n" )


if '--in' in sys.argv and '--ref' in sys.argv and '--target' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
