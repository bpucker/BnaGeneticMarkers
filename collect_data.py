### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """ python3 collect_data.py
						--gff <GFF_FILE>
						--athrbhs <ATH_RBH_FILE>
						--bnarbhs <BnaDarmorRBH_FILE>
						--anno <ANNO_FILE>
						--out <OUTPUT_FILE>
					"""


import sys, os

# --- end of imports --- #

def load_gene_positions( gff_file ):
	"""! @brief load gene positions """
	
	gene_positions = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 2:
				if parts[2] == "mRNA":
					gene_positions.append( { 'chr': parts[0], 'start': parts[3], 'end': parts[4], 'id': parts[-1].split('ID=')[1].split('_BnaEXP')[0] } )
			line = f.readline()
	return gene_positions


def load_ath_rbhs( rbhs_vs_ath ):
	"""! @brief load RBHs """
	
	ath_rbhs = {}
	with open( rbhs_vs_ath, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			ath_rbhs.update( { parts[0].replace( "_BnaEXP", "" ): parts[1].split('.')[0] } )
			line = f.readline()
	return ath_rbhs


def load_darmor_rbhs( rbhs_vs_darmorv51 ):
	"""! @brief load RBHs """
	
	raps_rbhs = {}
	with open( rbhs_vs_darmorv51, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			raps_rbhs.update( { parts[0].replace( "_BnaEXP", "" ): parts[1] } )
			line = f.readline()
	return raps_rbhs
	

def load_annotation( ath_anno_file ):
	"""! @brief load annotation file """
	
	anno = {}
	with open( ath_anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno.update( { parts[0]: ";".join( parts[1:] ) } )
			line = f.readline()
	return anno


def main( arguments ):
	"""! @brief run everything """
	
	gff_file = arguments[ arguments.index('--gff')+1 ]
	rbhs_vs_ath = arguments[ arguments.index('--athrbhs')+1 ]
	rbhs_vs_darmorv51 = arguments[ arguments.index('--bnarbhs')+1 ]
	ath_anno_file = arguments[ arguments.index('--anno')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]


	gene_positions = load_gene_positions( gff_file )
	raps_to_ath_mapping = load_ath_rbhs( rbhs_vs_ath )
	raps_to_darmor_mapping = load_darmor_rbhs( rbhs_vs_darmorv51 )
	annotation_mapping_table = load_annotation( ath_anno_file )

	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "GeneID", "Chromosome", "Start", "End", "AGI", "DarmorBzhV5.1", "AGI_Annotation" ] ) + "\n" )
		for gene in gene_positions:
			try:
				try:
					new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], raps_to_ath_mapping[ gene['id'] ], raps_to_darmor_mapping[ gene['id'] ], annotation_mapping_table[ raps_to_ath_mapping[ gene['id'] ] ] ]
				except KeyError:
					new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], raps_to_ath_mapping[ gene['id'] ], raps_to_darmor_mapping[ gene['id'] ], "n/a" ]
			except KeyError:
				try:
					new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], "n/a", raps_to_darmor_mapping[ gene['id'] ], "n/a" ]
				except KeyError:
					try:
						try:
							new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], raps_to_ath_mapping[ gene['id'] ], "n/a", annotation_mapping_table[ raps_to_ath_mapping[ gene['id'] ] ] ]
						except KeyError:
							new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], raps_to_ath_mapping[ gene['id'] ], "n/a", "n/a" ]
					except KeyError:
						new_line = [ gene['id'], gene['chr'], gene['start'], gene['end'], "n/a", "n/a", "n/a" ]
			
			out.write( "\t".join( new_line ) + "\n" )


if '--gff' in sys.argv and '--athrbhs' in sys.argv and '--bnarbhs' in sys.argv and '--anno' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
