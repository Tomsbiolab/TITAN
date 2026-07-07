#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* ---------------------------------------------------
 * SRA & FASTQ PROCESSING PROCESSES
 * --------------------------------------------------- */

process SRA_DOWNLOAD {
    tag "$group_id:$sra_id"

    input:
    tuple val(group_id), val(sra_id)

    output:
    tuple val(group_id), val(sra_id), path("*.fastq")

    script:
    """
    # Download SRA and split to fastq
    prefetch ${sra_id} --max-size 100G
    fastq-dump -I --split-files ${sra_id}

    # Remove the prefetched .sra (and its folder): it is dead weight once the
    # FASTQ files have been extracted and is never used by downstream steps.
    rm -rf ${sra_id}
    """
}

process FASTP {
    tag "$group_id:$sample_id"
    publishDir "${params.outdir}/fastp_reports", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(group_id), val(sample_id), path(reads)

    output:
    tuple val(group_id), val(sample_id), path("*_trim.fq.gz"), emit: trimmed_reads
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
    tag "$group_id:$sample_id"
    // Copy FASTQ files instead of symlinking to avoid Docker volume
    // sync lag on ext4 external HDDs serving stale data through symlinks
    stageInMode 'copy'

    input:
    tuple val(group_id), val(sample_id), path(reads)
    path star_index

    output:
    tuple val(group_id), val(sample_id), path("*.bam")

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
    tag "$group_id:$sample_id"

    input:
    tuple val(group_id), val(sample_id), path(bam)

    output:
    tuple val(group_id), val(sample_id), path("${sample_id}.filtered.bam")

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

process MERGE_SORT_INDEX_BAM {
    // Merge all per-sample filtered BAMs that share the same group_id (i.e. the same
    // input source: one fastq_dir folder or one SRA CSV file) into a single
    // coordinate-sorted, indexed BAM. samtools merge already does a merge-sort of
    // coordinate-sorted inputs, so no extra sort pass is needed.
    tag "$group_id"

    input:
    tuple val(group_id), path(bams)

    output:
    tuple val(group_id), path("${group_id}.merged.bam"), path("${group_id}.merged.bam.bai")

    script:
    def bam_list = (bams instanceof List ? bams : [bams]).join(' ')
    """
    samtools merge -f -@ ${task.cpus} ${group_id}.merged.bam ${bam_list}
    samtools index -@ ${task.cpus} ${group_id}.merged.bam
    """
}

process ASSEMBLE_STRINGTIE {
    tag "$group_id"
    publishDir "${params.outdir}/transcriptomes/stringtie", mode: 'copy'

    input:
    tuple val(group_id), path(bam)

    output:
    path "*_stringtie*.gtf"

    script:
    """
    stringtie -p ${task.cpus} -o ${group_id}_stringtie_default.gtf ${bam}
    stringtie -f 0.99 -m 120 -a 15 -j 3 -c 3 -s 4.75 -g 50 -p ${task.cpus} -t -o ${group_id}_stringtie_morus.gtf ${bam}
    """
}

process ASSEMBLE_PSICLASS {
    tag "$group_id"
    publishDir "${params.outdir}/transcriptomes/psiclass", mode: 'copy'

    input:
    tuple val(group_id), path(bam)

    output:
    path "*_psiclass.gtf_vote.gtf"

    script:
    """
    psiclass -p ${task.cpus} -b ${bam} -o ${group_id}_psiclass.gtf
    """
}

/* ---------------------------------------------------
 * LONG-READ PROCESSES (minimap2 pathway)
 * --------------------------------------------------- */

process SRA_DOWNLOAD_LONG_READS {
    tag "$group_id:$sra_id"

    input:
    tuple val(group_id), val(sra_id)

    output:
    tuple val(group_id), val(sra_id), path("*.fastq")

    script:
    """
    # Download SRA — long reads are always single-end
    prefetch ${sra_id} --max-size 100G
    fastq-dump ${sra_id}

    # Remove the prefetched .sra (and its folder): it is dead weight once the
    # FASTQ has been extracted and is never used by downstream steps.
    rm -rf ${sra_id}
    """
}

process MINIMAP2_ALIGN {
    tag "$group_id:$sample_id"

    input:
    tuple val(group_id), val(sample_id), path(reads)
    path genome

    output:
    tuple val(group_id), val(sample_id), path("${sample_id}.sam")

    script:
    def reads_list = (reads instanceof List ? reads : [reads]).join(' ')
    """
    minimap2 -a -t ${task.cpus} --secondary=no -x splice ${genome} ${reads_list} -o ${sample_id}.sam
    """
}

process SORT_BAM_LONG_READS {
    tag "$group_id:$sample_id"

    input:
    tuple val(group_id), val(sample_id), path(bam)

    output:
    tuple val(group_id), val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o ${sample_id}.sorted.bam
    """
}

process MERGE_SORT_INDEX_BAM_LONG_READS {
    tag "$group_id"

    input:
    tuple val(group_id), path(bams)

    output:
    tuple val(group_id), path("${group_id}.lr_merged.bam"), path("${group_id}.lr_merged.bam.bai")

    script:
    def bam_list = (bams instanceof List ? bams : [bams]).join(' ')
    """
    samtools merge -f -@ ${task.cpus} ${group_id}.lr_merged.bam ${bam_list}
    samtools index -@ ${task.cpus} ${group_id}.lr_merged.bam
    """
}

process ASSEMBLE_STRINGTIE_LONG_READS {
    tag "$group_id"
    publishDir "${params.outdir}/transcriptomes/stringtie_long_reads", mode: 'copy'

    input:
    tuple val(group_id), path(bam)

    output:
    path "*_stringtie*.gtf"

    script:
    """
    stringtie -p ${task.cpus} -L -o ${group_id}_stringtie_lr_default.gtf ${bam}
    stringtie -f 0.99 -m 120 -a 15 -j 3 -c 3 -s 4.75 -g 50 -p ${task.cpus} -t -L \
        -o ${group_id}_stringtie_lr_morus.gtf ${bam}
    """
}

/* ---------------------------------------------------
 * REPEAT MASKING (EDTA)
 * --------------------------------------------------- */

process RUN_EDTA {
    tag "edta"
    publishDir "${params.outdir}/edta", mode: 'copy'

    input:
    path genome

    output:
    path "${genome}.mod.MAKER.masked", emit: masked_genome

    script:
    """
    EDTA.pl \
        --genome ${genome} \
        --species others \
        --step all \
        --sensitive 1 \
        --anno 1 \
        --threads ${task.cpus}
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
    path "braker/**/augustus.hints.gtf", emit: augustus_gtf
    path "braker/GeneMark-ETP/genemark.gtf", emit: genemark_gtf

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
    aegis_merge.py \
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
    aegis_finalize.py \
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
    genome_ch         = Channel.fromPath(params.genome, checkIfExists: true)

    // Masked genome: run EDTA if not provided, otherwise use the supplied path
    if (params.masked_genome) {
        masked_genome_ch = Channel.fromPath(params.masked_genome, checkIfExists: true)
    } else {
        log.info "No masked genome provided: Running EDTA to generate one."
        masked_genome_ch = RUN_EDTA(genome_ch)
    }
    fasta_dbs_ch      = Channel.fromPath("${params.fasta_databases_dir}/*.{fa,fasta}", checkIfExists: true)
    dedup_proteins_ch = DEDUP_PROTEINS(fasta_dbs_ch)
    diamond_dbs_ch    = MAKE_DIAMOND_DB(dedup_proteins_ch)
    local_ab_initio_ch = params.ab_initio ? Channel.fromPath(params.ab_initio, checkIfExists: true) : Channel.empty()

    // 1. Gather Initial Local Annotations (if any)
    local_transcriptomes_ch = params.transcriptome ? Channel.fromPath(params.transcriptome) : Channel.empty()

    // 2. Gather RAW reads (SRA and/or Local FASTQ)
    //    Every item emitted is a 3-tuple: [group_id, sample_id, reads]
    //    group_id = CSV basename without extension (SRA) or folder basename (FASTQ dir)
    raw_reads_ch = Channel.empty()

    if (params.sra_list) {
        // Parse each CSV file and tag every SRR ID with its source CSV basename as group_id.
        // Reads the file with plain Groovy so the group_id is available inside flatMap.
        sra_group_runs_ch = Channel.fromPath(params.sra_list)
            .flatMap { csv_file ->
                def group_id = csv_file.baseName
                def lines    = csv_file.readLines()
                if (!lines) return []
                def headers  = lines[0].split('\t')*.trim()
                def run_idx  = headers.indexOf('Run')
                lines.drop(1)
                    .findAll { it.trim() }
                    .collect { line -> [group_id, line.split('\t')[run_idx].trim()] }
            }
        // SRA_DOWNLOAD: input [group_id, sra_id] → output [group_id, sra_id, path(*.fastq)]
        raw_reads_ch = raw_reads_ch.mix(SRA_DOWNLOAD(sra_group_runs_ch))
    }

    if (params.fastq_dir) {
        // Build glob patterns for every supplied directory; fromFilePairs handles both PE and SE.
        // The parent folder name of the first file in each pair becomes the group_id.
        def fastq_dirs = params.fastq_dir instanceof List ? params.fastq_dir : [params.fastq_dir]
        def patterns   = fastq_dirs.collect { "${it}/*{1,2}.{fastq,fq}*" }

        Channel.fromFilePairs(patterns, size: -1)
            .map { sample_id, reads ->
                def group_id = reads[0].parent.name
                [group_id, sample_id, reads]
            }
            | set { local_fastq_ch }
        raw_reads_ch = raw_reads_ch.mix(local_fastq_ch)
    }

    // 3. Process the Transcriptomic raw data if it exists
    assembled_transcriptomes_ch = Channel.empty()

    // Create STAR index automatically (cached if already exists)
    star_idx_ch = STAR_INDEX(genome_ch)

    if (params.sra_list || params.fastq_dir) {
        // Trim → align → filter, propagating group_id through every step
        trimmed_ch  = FASTP(raw_reads_ch).trimmed_reads
        aligned_ch  = STAR_ALIGN(trimmed_ch, star_idx_ch.first())
        filtered_ch = PROCESS_BAM(aligned_ch)

        // Collect per-sample filtered BAMs by group_id, then produce one merged,
        // sorted, indexed BAM per input source (fastq_dir folder or SRA CSV file).
        grouped_bams_ch = filtered_ch
            .map   { group_id, sample_id, bam -> [group_id, bam] }
            .groupTuple()
        merged_bams_ch = MERGE_SORT_INDEX_BAM(grouped_bams_ch)

        // Each merged BAM is assembled independently → one GTF set per input source
        merged_bam_only_ch = merged_bams_ch.map { group_id, bam, bai -> [group_id, bam] }
        gtfs_stringtie = ASSEMBLE_STRINGTIE(merged_bam_only_ch)
        gtfs_psiclass  = ASSEMBLE_PSICLASS(merged_bam_only_ch)
        assembled_transcriptomes_ch = gtfs_stringtie.mix(gtfs_psiclass)

        // Collect ALL merged BAMs for the single BRAKER execution
        bams_to_braker_ch = merged_bams_ch.map { group_id, bam, bai -> bam }.collect()
    } else {
        bams_to_braker_ch = Channel.fromPath(params.genome) // Dummy to avoid empty path error
    }

    // --- Long-read processing (minimap2 pathway) ---
    //     Parallel to the short-read path above: SRA/FASTQ → align → sort → merge → StringTie -L
    //     No trimming, no deduplication, no quality filtering, no BRAKER.
    lr_raw_reads_ch  = Channel.empty()
    lr_assembled_ch  = Channel.empty()

    if (params.lr_sra_list) {
        // Parse each CSV file the same way as short-read SRA lists
        lr_sra_group_runs_ch = Channel.fromPath(params.lr_sra_list)
            .flatMap { csv_file ->
                def group_id = csv_file.baseName
                def lines    = csv_file.readLines()
                if (!lines) return []
                def headers  = lines[0].split('\t')*.trim()
                def run_idx  = headers.indexOf('Run')
                lines.drop(1)
                    .findAll { it.trim() }
                    .collect { line -> [group_id, line.split('\t')[run_idx].trim()] }
            }
        lr_raw_reads_ch = lr_raw_reads_ch.mix(SRA_DOWNLOAD_LONG_READS(lr_sra_group_runs_ch))
    }

    if (params.lr_fastq_dir) {
        // Long reads: single files, not paired-end — use a broad glob
        def lr_dirs    = params.lr_fastq_dir instanceof List ? params.lr_fastq_dir : [params.lr_fastq_dir]
        def lr_patterns = lr_dirs.collect { "${it}/*.{fastq,fq,fasta,fa}*" }
        Channel.fromPath(lr_patterns)
            .map { file ->
                def group_id  = file.parent.name
                def sample_id = file.baseName
                [group_id, sample_id, file]
            }
            | set { lr_local_ch }
        lr_raw_reads_ch = lr_raw_reads_ch.mix(lr_local_ch)
    }

    if (params.lr_sra_list || params.lr_fastq_dir) {
        // Align → sort → merge → assemble (no trim, no dedup)
        lr_aligned_ch = MINIMAP2_ALIGN(lr_raw_reads_ch, genome_ch.first())
        lr_sorted_ch  = SORT_BAM_LONG_READS(lr_aligned_ch)

        lr_grouped_ch = lr_sorted_ch
            .map { group_id, sample_id, bam -> [group_id, bam] }
            .groupTuple()
        lr_merged_ch  = MERGE_SORT_INDEX_BAM_LONG_READS(lr_grouped_ch)

        lr_bam_only_ch = lr_merged_ch.map { group_id, bam, bai -> [group_id, bam] }
        lr_assembled_ch = ASSEMBLE_STRINGTIE_LONG_READS(lr_bam_only_ch)
    }

    // 4. Merge All Transcriptomes (Local annotations + Short-read + Long-read)
    all_transcriptomes_ch = local_transcriptomes_ch
        .mix(assembled_transcriptomes_ch)
        .mix(lr_assembled_ch)
        .collect()

    // 5. Run BRAKER once with all merged BAMs → single augustus + genemark output
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

    // .first() turns the single merged-proteins emission into a value channel so it
    // is reused for EVERY diamond database, yielding one *_merged.diamond per DB
    // (Araport11, Viridiplantae, Eudicotyledons) instead of a single result.
    RUN_DIAMOND(
        AEGIS_MERGE_ANNOTATIONS.out.proteins.first(),
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
