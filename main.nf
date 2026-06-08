#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* ---------------------------------------------------
 * SRA & FASTQ PROCESSING PROCESSES
 * --------------------------------------------------- */

process SRA_DOWNLOAD {
    tag "$sra_id"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("*.fastq")

    script:
    """
    # Download SRA and split to fastq
    prefetch ${sra_id} --max-size 100G
    fastq-dump -I --split-files ${sra_id}
    """
}

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp_reports", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_trim.fq.gz"), emit: trimmed_reads
    path "*.html"
    path "*.json"

    script:
    // Groovy dynamically checks if input is a single file (SE) or a list (PE)
    def is_paired = reads instanceof List || reads.getClass().isArray() ? true : false
    
    if (is_paired) {
        """
        fastp \
            -w ${task.cpus} \
            -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \
            --detect_adapter_for_pe --n_base_limit 5 \
            --cut_front --cut_front_window_size 1 --cut_front_mean_quality 30 \
            --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30 \
            -l 20 \
            -i ${reads[0]} -I ${reads[1]} \
            -o ${sample_id}_1_trim.fq.gz -O ${sample_id}_2_trim.fq.gz
        """
    } else {
        """
        fastp \
            -w ${task.cpus} \
            -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \
            --n_base_limit 5 \
            --cut_front --cut_front_window_size 1 --cut_front_mean_quality 30 \
            --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30 \
            -l 20 \
            -i ${reads} \
            -o ${sample_id}_1_trim.fq.gz
        """
    }
}

process STAR_INDEX {
    storeDir "${file(params.genome).toAbsolutePath().getParent()}"

    input:
    path genome

    output:
    path "star_index_${genome.baseName}"

    script:
    """
    mkdir star_index_${genome.baseName}
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeSAindexNbases 12 \
         --genomeDir star_index_${genome.baseName} \
         --genomeFastaFiles ${genome}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    // Copy FASTQ files instead of symlinking to avoid Docker volume
    // sync lag on ext4 external HDDs serving stale data through symlinks
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(reads)
    path star_index

    output:
    tuple val(sample_id), path("*.bam")

    script:
    """
    # -------------------------------------------------------
    # Guard: validate that every input FASTQ is a readable gzip
    # stream before launching STAR.  On external HDDs with Docker
    # volume mounts, even physically copied files can briefly be
    # incomplete if the kernel write-back hasn't finished yet.
    # -------------------------------------------------------
    uncompressed_reads=""
    for f in ${reads}; do
        retries=0
        while ! gzip -t "\$f" 2>/dev/null && [ \$retries -lt 30 ]; do
            echo "Waiting for \$f gzip stream to be valid (attempt \$((retries+1)))..."
            sleep 2
            retries=\$((retries+1))
        done
        if ! gzip -t "\$f" 2>/dev/null; then
            echo "ERROR: \$f is not a valid gzip file after 60 s — aborting." >&2
            exit 1
        fi
        # Quick sanity check: first line must start with @ (FASTQ)
        first_char=\$(zcat "\$f" | head -c1)
        if [ "\$first_char" != "@" ]; then
            echo "ERROR: \$f does not look like a FASTQ file (first char='\$first_char')." >&2
            exit 1
        fi
        
        # Decompress to disk to avoid STAR's pipe short-read bugs
        uncompressed="\${f%.gz}"
        echo "Decompressing \$f to \$uncompressed ..."
        zcat "\$f" > "\$uncompressed"
        uncompressed_reads="\$uncompressed_reads \$uncompressed"
    done

    STAR --runMode alignReads \
         --genomeDir ${star_index} \
         --twopassMode Basic \
         --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
         --outSAMstrandField intronMotif \
         --runThreadN ${task.cpus} \
         --readFilesIn \$uncompressed_reads \
         --outSAMunmapped Within \
         --outSAMtype BAM Unsorted \
         --alignSJoverhangMin 20 \
         --alignSJDBoverhangMin 1 \
         --outFilterMultimapNmax 20 \
         --alignIntronMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --outFileNamePrefix ${sample_id}.
         
    # Cleanup uncompressed files to save disk space
    rm -f \$uncompressed_reads
    """
}

process PROCESS_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam")

    script:
    // Piped samtools commands to drastically save IO operations and disk space
    """
    samtools sort -n -@ ${task.cpus} ${bam} | \
    samtools fixmate -m -@ ${task.cpus} - - | \
    samtools sort -@ ${task.cpus} - | \
    samtools markdup -r -@ ${task.cpus} - - | \
    samtools view -@ ${task.cpus} -b -q 20 -o ${sample_id}.filtered.bam -
    """
}

process ASSEMBLE_STRINGTIE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "*_stringtie*.gtf"

    script:
    """
    stringtie -p ${task.cpus} -o ${sample_id}_stringtie_default.gtf ${bam}
    stringtie -f 0.99 -m 120 -a 15 -j 3 -c 3 -s 4.75 -g 50 -p ${task.cpus} -t -o ${sample_id}_stringtie_morus.gtf ${bam}
    """
}

process ASSEMBLE_PSICLASS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "*_psiclass.gtf_vote.gtf"

    script:
    """
    psiclass -p ${task.cpus} -b ${bam} -o ${sample_id}_psiclass.gtf
    """
}

/* ---------------------------------------------------
 * ANNOTATION & FILTERING PROCESSES (Original)
 * --------------------------------------------------- */

process DEDUP_PROTEINS {
    tag "$fasta"
    publishDir "${params.outdir}/dedup_proteins", mode: 'copy'

    input:
    path fasta

    output:
    path "dedup/${fasta.name}"

    script:
    """
    mkdir -p dedup
    seqkit rmdup -s ${fasta} 2> ${fasta.baseName}.rmdup_seq.log \
        | seqkit rmdup -n -o dedup/${fasta.name} 2> ${fasta.baseName}.rmdup_name.log
    """
}

process MERGE_PROTEINS_FOR_BRAKER {
    tag "merge_proteins"
    publishDir "${params.outdir}/dedup_proteins", mode: 'copy', pattern: "braker_proteins.fasta"

    input:
    path proteins

    output:
    path "braker_proteins.fasta"

    script:
    // Concatenate in source_priority order so cross-database ID duplicates
    // keep the highest-priority entry (seqkit rmdup -n retains the first occurrence).
    def prot_files = proteins instanceof List ? proteins : [proteins]
    def priority = params.source_priority.split(',').collect { it.trim() }
    def ordered = []
    priority.each { p ->
        def match = prot_files.find { f -> f.baseName == p || f.name.startsWith(p) }
        if (match) ordered << match
    }
    def remaining = prot_files.findAll { f -> !ordered.contains(f) }
    def cat_files = (ordered + remaining).join(' ')
    """
    cat ${cat_files} > all_proteins.fasta
    seqkit rmdup -s all_proteins.fasta 2> rmdup_seq.log \
        | seqkit rmdup -n -o braker_proteins.fasta 2> rmdup_name.log
    """
}

process MAKE_DIAMOND_DB {
    tag "$fasta"
    storeDir "${params.fasta_databases_dir}"

    input:
    path fasta

    output:
    path "${fasta.baseName}.dmnd"

    script:
    """
    diamond makedb -p ${task.cpus} --in ${fasta} -d ${fasta.baseName}
    """
}

process RUN_BRAKER {
    tag "braker"
    publishDir "${params.outdir}/braker", mode: 'copy'

    input:
    path genome
    path bams
    path proteins

    output:
    path "braker/augustus.hints.gtf", emit: augustus_gtf
    path "braker/**/genemark.gtf", emit: genemark_gtf

    script:
    def bam_list = (bams instanceof List ? bams : [bams]).findAll { it.name.endsWith('.bam') }
    def prot_list = proteins instanceof List ? proteins : [proteins]
    def bam_input = bam_list ? "--bam=" + bam_list.join(',') : ""
    def prot_input = prot_list ? "--prot_seq=" + prot_list.join(',') : ""
    
    """
    braker.pl  \
        --genome=${genome}  \
        ${bam_input}  \
        ${prot_input}  \
        --threads=${task.cpus}  \
        --workingdir=braker  \
        --softmasking_off  \
        --gff3
    """
}

process AEGIS_MERGE_ANNOTATIONS {
    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path genome
    path ab_initio_files
    path transcriptome_files
    val version

    output:
    path "**/merge_${version}.gff3", emit: gff
    path "*_unique_proteins.fasta", emit: proteins

    script:
    """
    cat << 'EOF' > merge.py
    import argparse
    import os
    from aegis.annotation import Annotation
    from aegis.genome import Genome

    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', required=True)
    parser.add_argument('--ab_initio', nargs='+', required=True)
    parser.add_argument('--transcriptome', nargs='+', required=True)
    parser.add_argument('--version', required=True)
    args = parser.parse_args()

    genome_obj = Genome('RefGenome', args.genome)

    def get_name(filepath):
        return os.path.splitext(os.path.basename(filepath))[0]

    print('Processing Ab Initio evidences...')
    ab_initio_annotations = [
        Annotation(name=get_name(f), annot_file_path=f, genome=genome_obj, rename_source=get_name(f)) 
        for f in args.ab_initio
    ]

    print('Processing Transcriptome evidences...')
    transcriptome_annotations = [
        Annotation(name=get_name(f), annot_file_path=f, genome=genome_obj, rename_source=get_name(f), rework_all_CDSs=True) 
        for f in args.transcriptome
    ]

    print('Merging Transcriptome evidences...')
    transcript_evidence = transcriptome_annotations[0].copy()
    for annotation in transcriptome_annotations[1:]:
        transcript_evidence.merge(annotation)

    print('Merging all evidences...')
    AEGIS_MERGE_ANNOTATIONS = ab_initio_annotations[0].copy()
    for annotation in ab_initio_annotations[1:]:
        AEGIS_MERGE_ANNOTATIONS.merge(annotation)
    
    AEGIS_MERGE_ANNOTATIONS.merge(transcript_evidence)
    AEGIS_MERGE_ANNOTATIONS.make_alternative_transcripts_into_genes()

    AEGIS_MERGE_ANNOTATIONS.id = f'merge_{args.version}_on_genome'
    AEGIS_MERGE_ANNOTATIONS.name = f'merge_{args.version}'

    print('Exporting results...')
    AEGIS_MERGE_ANNOTATIONS.export.gff('.', f'merge_{args.version}.gff3')
    AEGIS_MERGE_ANNOTATIONS.export.unique_proteins(custom_path='.')
    EOF

    python3 merge.py \
        --genome ${genome} \
        --ab_initio ${ab_initio_files} \
        --transcriptome ${transcriptome_files} \
        --version ${version}
    """
}

process RUN_DIAMOND {
    publishDir "${params.outdir}/blast_results", mode: 'copy'

    input:
    path proteins
    path db

    output:
    path "${db.baseName}_merged.diamond"

    script:
    """
    diamond blastp \
        --threads ${params.threads} \
        --db ${db} \
        --ultra-sensitive \
        --out ${db.baseName}_merged.diamond \
        --outfmt 6 \
        --query ${proteins} \
        --max-target-seqs 1 \
        --evalue 1e-3
    """
}

process AEGIS_FILTER_REDUNDANCY {
    publishDir "${params.outdir}/final_annotation", mode: 'copy'

    input:
    path genome
    path masked_genome
    path merged_gff
    path blast_results
    val priority_list
    val version

    output:
    path "**/final_annotation_${version}.gff3"
    path "*.csv"
    path "*.fasta"
    path "*.pkl"
    path "**/merge_overlaps_${version}.gff3"

    script:
    """
    cat << 'EOF' > finalize.py
    import argparse
    import os
    from aegis.utils.misc import pickle_save
    from aegis.annotation import Annotation
    from aegis.genome import Genome

    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', required=True)
    parser.add_argument('--masked_genome', required=True)
    parser.add_argument('--merged_gff', required=True)
    parser.add_argument('--blast_results', nargs='+', required=True)
    parser.add_argument('--priority', required=True)
    parser.add_argument('--version', required=True)
    args = parser.parse_args()

    genome_obj = Genome('RefGenome', args.genome)
    masked_genome_obj = Genome('MaskedGenome', args.masked_genome)
    
    merge = Annotation(
        name='merge', 
        annot_file_path=args.merged_gff, 
        genome=genome_obj, 
        hard_masked_genome=masked_genome_obj
    )

    merge.stats.calculate_transcript_masking()
    merge.update()
    merge.overlaps.detect()

    print('Adding BLAST results...')
    for blast_file in args.blast_results:
        source_name = blast_file.replace('_merged.diamond', '')
        print(f"Loading BLAST hits for {source_name}...")
        merge.add_blast_hits(source_name, blast_file)

    print('Saving intermediate overlaps...')
    pickle_save(f'merge_overlaps_{args.version}.pkl', merge)
    merge.update()
    merge.export.gff('.', f'merge_overlaps_{args.version}.gff3')

    print('Reducing redundancy...')
    priority_sources = args.priority.split(',')
    merge.redundancy.filter(source_priority=priority_sources)

    merge.id = f'final_annotation_{args.version}_on_genome'
    merge.name = f'final_annotation_{args.version}'

    print('Exporting final outputs...')
    merge.export.gff('.', f'final_annotation_{args.version}.gff3')
    merge.overlaps.export(overlap_threshold=0, export_csv=True, export_self=True, NAs=False, custom_path='.')
    merge.export.proteins(only_main=True, verbose=False, custom_path='.')
    merge.export.proteins(only_main=False, verbose=False, custom_path='.')
    EOF

    python3 finalize.py \
        --genome ${genome} \
        --masked_genome ${masked_genome} \
        --merged_gff ${merged_gff} \
        --blast_results ${blast_results} \
        --priority "${priority_list}" \
        --version ${version}
    """
}

/* ---------------------------------------------------
 * MAIN WORKFLOW LOGIC
 * --------------------------------------------------- */

workflow {
    // Log SRA download environment configuration
    if (params.sra_list) {
        if (params.has_prefetch) {
            log.info "SRA Toolkit detected on host. Running SRA_DOWNLOAD natively (without Docker)."
        } else {
            log.info "SRA Toolkit not found on host. Running SRA_DOWNLOAD inside Docker container."
        }
    }

    // Core Reference Channels
    genome_ch        = Channel.fromPath(params.genome, checkIfExists: true)
    masked_genome_ch = Channel.fromPath(params.masked_genome, checkIfExists: true)
    fasta_dbs_ch      = Channel.fromPath("${params.fasta_databases_dir}/*.{fa,fasta}", checkIfExists: true)
    dedup_proteins_ch = DEDUP_PROTEINS(fasta_dbs_ch)
    diamond_dbs_ch    = MAKE_DIAMOND_DB(dedup_proteins_ch)
    local_ab_initio_ch = params.ab_initio ? Channel.fromPath(params.ab_initio, checkIfExists: true) : Channel.empty()
    
    // 1. Gather Initial Local Annotations (if any)
    local_transcriptomes_ch = params.transcriptome ? Channel.fromPath(params.transcriptome) : Channel.empty()

    // 2. Gather RAW reads (SRA and/or Local FASTQ)
    raw_reads_ch = Channel.empty()

    if (params.sra_list) {
        // Read lines from file, ignoring header/blanks. Assuming CSV with SRR on the first column
        Channel.fromPath(params.sra_list)
               .splitCsv(sep: '\t', header: true) // adjust if comma separated
               .map { row -> row.Run } // Assuming column is named 'Run', adjust if needed
               | SRA_DOWNLOAD
               | set { sra_reads_ch }
        raw_reads_ch = raw_reads_ch.mix(sra_reads_ch)
    }

    if (params.fastq_dir) {
        // Handle both a single string and a list of directories
        def fastq_patterns = params.fastq_dir instanceof List 
            ? params.fastq_dir.collect { "${it}/*{1,2}.{fastq,fq}*" }
            : "${params.fastq_dir}/*{1,2}.{fastq,fq}*"

        // Using fromFilePairs handles both PE (e.g. _1.fq, _2.fq) and SE (-1 size) 
        Channel.fromFilePairs(fastq_patterns, size: -1)
               | set { local_fastq_ch }
        raw_reads_ch = raw_reads_ch.mix(local_fastq_ch)
    }

    // 3. Process the Transcriptomic raw data if it exists
    assembled_transcriptomes_ch = Channel.empty()
    
    // Create STAR index automatically (cached if already exists)
    star_idx_ch = STAR_INDEX(genome_ch)

    if (params.sra_list || params.fastq_dir) {
        // Run preprocessing and alignment
        trimmed_ch  = FASTP(raw_reads_ch).trimmed_reads
        aligned_ch  = STAR_ALIGN(trimmed_ch, star_idx_ch.first())
        filtered_ch = PROCESS_BAM(aligned_ch)
        
        // Assemble Transcripts
        gtfs_stringtie = ASSEMBLE_STRINGTIE(filtered_ch)
        gtfs_psiclass  = ASSEMBLE_PSICLASS(filtered_ch)
        
        // Combine outputs of assemblies
        assembled_transcriptomes_ch = gtfs_stringtie.mix(gtfs_psiclass)
        
        // Collect BAMs for Braker
        bams_to_braker_ch = filtered_ch.map { it[1] }.collect()
    } else {
        bams_to_braker_ch = Channel.fromPath(params.genome) // Dummy to avoid empty path error
    }

    // 4. Merge All Transcriptomes (Local annotations + Newly Assembled Data)
    all_transcriptomes_ch = local_transcriptomes_ch.mix(assembled_transcriptomes_ch).collect()

    // 5. Run Braker and collect all ab initio annotations
    braker_proteins_ch = MERGE_PROTEINS_FOR_BRAKER(dedup_proteins_ch.collect())
    RUN_BRAKER(
        genome_ch,
        bams_to_braker_ch,
        braker_proteins_ch
    )

    all_ab_initio_ch = local_ab_initio_ch.mix(RUN_BRAKER.out.augustus_gtf, RUN_BRAKER.out.genemark_gtf).collect()

    // 6. Run the Annotation Pipeline
    AEGIS_MERGE_ANNOTATIONS(
        genome_ch,
        all_ab_initio_ch,
        all_transcriptomes_ch, // Feeds dynamically from 0 to N inputs
        params.version
    )

    RUN_DIAMOND(
        AEGIS_MERGE_ANNOTATIONS.out.proteins,
        diamond_dbs_ch
    )

    AEGIS_FILTER_REDUNDANCY(
        genome_ch,
        masked_genome_ch,
        AEGIS_MERGE_ANNOTATIONS.out.gff,
        RUN_DIAMOND.out.collect(),
        params.source_priority,
        params.version
    )
}
