/*
 * Check input samplesheet or folder and get read channels
 */

params.options = [:]

include { parse_samplesheet } from '../../modules/local/parse_samplesheet' addParams( options: params.options )


if (params.barcodes) {
   	//PacBio data
	cutadapt_options_args       = " --rc -g file:${params.barcodes} -o data/{name}.fastq"	
}

def cutadapt_options 			= modules['cutadapt']
cutadapt_options.args          += cutadapt_options_args
cutadapt_options.args          += params.retain_untrimmed ? '' : " --discard-untrimmed"

include { CUTADAPT_WORKFLOW             } from '../subworkflows/local/cutadapt_workflow'        addParams( standard_options: cutadapt_options, readthrough_options: cutadapt_readthrough_options,doubleprimer_options: cutadapt_doubleprimer_options,summary_options: modules['cutadapt_summary'],summary_merge_options: modules['cutadapt_summary_merge'] )


workflow DEMULTIPLEX {
    take:
    input // file.tsv or folder
	single_end
    multiple_sequencing_runs
    extension

    main:

	if ( input.toString().toLowerCase().endsWith(".fasta") || input.toString().toLowerCase().endsWith(".fna") || input.toString().toLowerCase().endsWith(".fa") ) {
		// Fasta input not applicable for demultiplexing
		exit 1, "Fasta input can not be used for demultiplexing"
	} else if (!single_end) {
	  exit 1, "Demultiplexing only implemented for single end data,  so far"
	} else {
		ch_fasta = Channel.empty()


		if ( input.toString().toLowerCase().endsWith("tsv") ) {
			// Sample sheet input

			tsvFile = file(input).getName()
			// extracts read files from TSV and distribute into channels
			Channel
				.fromPath(input)
				.ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
				.splitCsv(header:true, sep:'\t')
				.map { parse_samplesheet(it, single_end) }
				.set { ch_reads }
		} else {
			// Folder input
			exit 1, "Demultiplexing not implemented yet for folder input"

			//Check folders in folder when multiple_sequencing_runs
		}








    /*
     * MODULE: Cutadapt
     */
    CUTADAPT_WORKFLOW ( 
		RENAME_RAW_DATA_FILES.out,
		params.illumina_pe_its,
		params.double_primer
	).reads.set { ch_trimmed_reads }
    ch_software_versions = ch_software_versions.mix(CUTADAPT_WORKFLOW.out.version.first().ifEmpty(null))






		//Check whether all sampleID = meta.id are unique
		ch_reads
			.map { meta, reads -> [ meta.id ] }
			.toList()
			.subscribe { 
				if( it.size() != it.unique().size() ) {
					ids = it.take(10);
					exit 1, "Please review data input, sample IDs are not unique! First IDs are $ids" 
				}
			}

		//Check that no dots "." are in sampleID
		ch_reads
			.map { meta, reads -> [ meta.id ] }
			.subscribe { if ( "$it".contains(".") ) exit 1, "Please review data input, sampleIDs may not contain dots, but \"$it\" does." }
	}

    emit:
    reads   = ch_reads
	fasta   = ch_fasta
}