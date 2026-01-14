import os
import glob
from datetime import datetime

os.environ["PYTHONNOUSERSITE"] = "1"

##### DIRECTORIES AND PARAMETERS ######
configfile: "config.yaml"
input_dir = config["folder"]["input_dir"]
output_dir = config["folder"]["output_dir"]
sample_tab = config["params"]["sample_tab"]
samples = config["params"]["sample"]
print(samples)
threads = config.get("threads", os.cpu_count())
dorado_pu = config["params"]["dorado_pu"]
filtlong_min_length=config["params"]["filtlong"]["min_length"]
filtlong_keep_percent=config["params"]["filtlong"]["keep_percent"]
racon_rounds=config["params"]["dragonflye"]["racon_rounds"]
db_dir = config["folder"]["db_dir"]
eggnog_db = config["params"]["eggnog"]["db"]
abricate_db_dir = config["params"].get("abricate_db_dir", "")
abricate_db_name = config["params"].get("abricate_db_name", "")
dmr_window = config["params"].get("dmr_window", 50)
dmr_pairs = config["params"].get("dmr_pairs", [])
dorado_fallback_models = config["params"]["dorado_fallback_models"]
model_tier = config["params"]["dorado_model_tier"]
mod_bases = config["params"]["dorado_modified_bases"]

assembler_choice = config["params"].get("assembler", "dragonflye")
genome_size = config["params"].get("genome_size", "5m")

onstart:
    print("Setting executable permissions for scripts...")
    shell("chmod +x scripts/*.py 2>/dev/null || true")
    shell("chmod +x scripts/*.sh 2>/dev/null || true")

def get_assembly_fasta(wildcards):
    if assembler_choice == "dragonflye":
        return f"{output_dir}/assemblies_evaluation/dragonflye/{wildcards.sample}/{wildcards.sample}.fa"
    elif assembler_choice == "canu":
        return f"{output_dir}/assemblies_evaluation/canu/{wildcards.sample}/{wildcards.sample}.fa"
    elif assembler_choice == "raven":
        return f"{output_dir}/assemblies_evaluation/raven/{wildcards.sample}/{wildcards.sample}.fa"
    else:
        raise ValueError(f"Assembler '{assembler_choice}' not supported.")

#### PARSING SAMPLE ####

SAMPLES = {}
with open(sample_tab, 'r') as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) >= 4:
            sample_id = parts[0]
            file_path = parts[2]
            file_type = parts[3].lower()
            
            SAMPLES[sample_id] = {
                "path": file_path,
                "type": file_type
            }

ALL_SAMPLES = list(SAMPLES.keys())

#### INPUT FUNCTIONS ####

def get_pod5_input(wildcards):
    s_info = SAMPLES[wildcards.sample]
    s_type = s_info["type"]
    
    if s_type == "fast5":
        return f"{output_dir}/basecalling/{wildcards.sample}/{wildcards.sample}.pod5"
    elif s_type == "pod5":
        return s_info["path"]
    else:
        return []

def get_bam_input(wildcards):
    s_info = SAMPLES[wildcards.sample]
    s_type = s_info["type"]
    
    if s_type == "bam":
        return s_info["path"]
    else:
        return f"{output_dir}/basecalling/{wildcards.sample}/{wildcards.sample}_basecalling.bam"

def get_fast5_input(wildcards):
    return SAMPLES[wildcards.sample]["path"]

start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(start_time)


#### RULES ####

rule all:
    input:
        expand("{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq", 
               output_dir=output_dir, sample=ALL_SAMPLES),
        
        #Quality check and filtering
        expand("{output_dir}/qc_filtering/nanoplot_before_filtering/{sample}_nanoplot_qc_report_b", 
               output_dir=output_dir, sample=ALL_SAMPLES),
        expand("{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq", 
               output_dir=output_dir, sample=ALL_SAMPLES),
        expand("{output_dir}/qc_filtering/nanoplot_after_filtering/{sample}_nanoplot_qc_report_a", 
               output_dir=output_dir, sample=ALL_SAMPLES),

        # Assembly
        expand("{output_dir}/assemblies_evaluation/{assembler}/{sample}/{sample}.fa", output_dir=output_dir, sample=ALL_SAMPLES, assembler=assembler_choice),
        expand("{output_dir}/assemblies_evaluation/quast/{sample}", output_dir=output_dir, sample=ALL_SAMPLES),

        # Annotation
        expand("{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa", output_dir=output_dir, sample=ALL_SAMPLES),
        
        # EggNog Annotation
        expand("{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations", output_dir=output_dir, sample=ALL_SAMPLES),

        # Abricate
        (expand("{output_dir}/annotation_feat/abricate/{sample}/{sample}_abricate.txt", output_dir=output_dir, sample=ALL_SAMPLES) 
         if config["params"].get("abricate_db_name") else []),

        f"{eggnog_db}",

        # Microbemod
        expand("{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam", output_dir=output_dir, sample=samples),
        expand("{output_dir}/microbemod/MicrobeMod/.installed", output_dir=output_dir),
        expand("{output_dir}/microbemod/results/{sample}/{sample}_methylated_sites.tsv", output_dir=output_dir, sample=samples),
        expand("{output_dir}/microbemod/results/{sample}/annotate_rm/{sample}.rm.genes.tsv", output_dir=output_dir, sample=samples),

        # MeStudio
        expand("{output_dir}/mestudio/{sample}/{sample}_smart.gff", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/{sample}/motiflist.txt", output_dir=output_dir, sample=samples),
        expand("{output_dir}/assemblies_evaluation/{sample}/{sample}_genomic.fasta", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/results/{sample}", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/results/{sample}_checkout.txt", output_dir=output_dir, sample=samples),

        # GSEA
        expand("{output_dir}/mestudio/results/{sample}_checkout_gsea.txt", output_dir=output_dir, sample=samples),

        # DMR analysis
        expand("{output_dir}/DMRs/{control}_vs_{treated}_volcano.png", output_dir=output_dir,
               control=[pair[0] for pair in dmr_pairs], treated=[pair[1] for pair in dmr_pairs]),
        expand("{output_dir}/DMRs/{control}_vs_{treated}_window_results.tsv", output_dir=output_dir,
               control=[pair[0] for pair in dmr_pairs], treated=[pair[1] for pair in dmr_pairs]),

        # Benchmark summary
        expand("{output_dir}/benchmarks/global_benchmark_data.tsv", output_dir=output_dir),
        expand("{output_dir}/benchmarks/global_benchmark_stats_summary.tsv", output_dir=output_dir)

rule convert_fast5_to_pod5:
    message: "Converting .fast5 files to POD5 file before Dorado basecalling."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        get_fast5_input
    output:
        "{output_dir}/basecalling/{sample}/{sample}.pod5"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_fast5_to_pod5_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/benchmarks/{wildcards.sample}
        mkdir -p {output_dir}/basecalling/{wildcards.sample}
        pod5 convert fast5 {input} --output {output}
        """

rule dorado_basecalling:
    message: "Dorado basecalling with Auto-Discovery"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        pod5 = get_pod5_input
    output:
        "{output_dir}/basecalling/{sample}/{sample}_basecalling.bam"
    threads: threads
    params:
        tier = model_tier,
        mods = mod_bases,
        pu = dorado_pu,
        fallback_models=expand(f"{output_dir}/dorado_models/{{model}}", model=dorado_fallback_models),
        fallback_model_names=dorado_fallback_models,
        mdir = "{output_dir}/basecalling/{sample}/dorado_models",
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_basecalling_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/benchmarks/{wildcards.sample}
        mkdir -p {params.mdir}

        export DORADO_MODELS_DIRECTORY="{params.mdir}"

        DORADO_BIN=$(bash scripts/install_dorado.sh | tail -n 1)

        bash scripts/install_dorado.sh

        if ! "$DORADO_BIN" basecaller {params.tier} {input.pod5} --modified-bases {params.mods} \
             -x {params.pu} > {output} 2> {output}.log_tmp; then

            rm -f {output}
            for model in {params.fallback_model_names}
            do
                "$DORADO_BIN" download --model $model --directory {params.mdir}
            done

            "$DORADO_BIN" basecaller {params.mdir}/{params.fallback_model_names[0]} {input.pod5} \
            --modified-bases 4mC_5mC 6mA \
            -x {params.pu} > {output}
        else
            rm -f {output}.log_tmp
        fi
        """


rule convert_bam_to_fastq:
    input:
        get_bam_input
    output:
        "{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_bam_to_fastq_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/benchmarks/{wildcards.sample}
        samtools fastq {input} > {output}
        """

#### QC & FILTERING ####

rule nanoplot_qc_b:
    message: "Quality check (Before Filtering): Executing NanoPlot."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        fastq_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    output:
        qc_report=directory("{output_dir}/qc_filtering/nanoplot_before_filtering/{sample}_nanoplot_qc_report_b")
    threads: threads
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_nanoplot_before_filtering_benchmark.txt"
    shell:
        """
        mkdir -p {output.qc_report}
        NanoPlot --fastq {input.fastq_file} --N50 -o {output.qc_report}
        """

rule filtlong_filter:
    message: "Filtering: Executing Filtlong to remove low-quality reads."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        fastq_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    output:
        filtered_fastq="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_filtlong_benchmark.txt"
    params:
        min_length=filtlong_min_length,
        keep_percent=filtlong_keep_percent
    shell:
        """
        mkdir -p {output_dir}/qc_filtering/filtered/{wildcards.sample}
        filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input.fastq_file} > {output.filtered_fastq}
        """

rule nanoplot_qc_a:
    message: "Quality check (After Filtering): Executing NanoPlot."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        filtered_fastq="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        qc_report_a=directory("{output_dir}/qc_filtering/nanoplot_after_filtering/{sample}_nanoplot_qc_report_a")
    threads: threads
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_nanoplot_after_filtering_benchmark.txt"
    shell:
        """
        mkdir -p {output.qc_report_a}
        NanoPlot --fastq {input.filtered_fastq} --N50 -o {output.qc_report_a}
        """
#### ASSEMBLY ####

rule dragonflye_assembly:
    message: "De novo genome assembly with Dragonflye."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        reads="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        assembly="{output_dir}/assemblies_evaluation/dragonflye/{sample}/{sample}.fa",
        dir=directory("{output_dir}/assemblies_evaluation/dragonflye/{sample}")
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_dragonflye_assembly.txt"
    conda:
        "envs/dragonflye.yaml"
    params:
        racon_rounds=racon_rounds
    priority: 20
    shell:
        """
        dragonflye \
            --outdir {output.dir} \
            --reads {input.reads} \
            --assembler flye \
            --nanohq \
            --racon {params.racon_rounds} \
            --prefix {wildcards.sample} \
            --dnaapler_mode all --force \
            --cpus {threads}
        """


rule canu_assembly:
    message: "De novo genome assembly with Canu."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        reads="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        assembly="{output_dir}/assemblies_evaluation/canu/{sample}/{sample}.fa",
        dir=directory("{output_dir}/assemblies_evaluation/canu/{sample}")
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_canu_assembly_benchmark.txt"
    params:
        gsize=genome_size
    shell:
        """
        mkdir -p {output.dir}
        
        canu \
        -p {wildcards.sample} -d {output.dir} \
        genomeSize={params.gsize} \
        -nanopore {input.reads} \
        useGrid=false maxThreads={threads}

        mv {output.dir}/{wildcards.sample}.contigs.fasta {output.assembly}
        """

rule raven_assembly:
    message: "De novo genome assembly with Raven."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        reads="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        assembly="{output_dir}/assemblies_evaluation/raven/{sample}/{sample}.fa",
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_raven_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/assemblies_evaluation/raven/{wildcards.sample}
        
        raven \
            --threads {threads} \
            {input.reads} > {output.assembly}
        """


#### EVALUATION ####

rule quast:
    message: "Quality assessment using QUAST"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly=get_assembly_fasta 
    output:
        quast_report_dir=directory("{output_dir}/assemblies_evaluation/quast/{sample}")
    threads: threads
    conda:
        "envs/quast.yaml"
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_quast_benchmark.txt"
    shell:
        """
        quast.py {input.assembly} -o {output.quast_report_dir} --threads {threads} > /dev/null 2>&1
        """

#### ANNOTATION ####

rule prokka:
    message: "Annotate genome with Prokka"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly=get_assembly_fasta
    output:
        faa="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa",
        gbk="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gbk",
        gff="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gff"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_prokka_benchmark.txt"
    conda:
        "envs/prokka.yaml"
    params:
        prefix="{sample}_anno"
    shell:
        """
        prokka --outdir {output_dir}/annotation_feat/prokka/{wildcards.sample} \
               --prefix {params.prefix} \
               --force \
               --cpus {threads} \
               {input.assembly}
        """

rule download_eggnog_db:
    output:
        eggnog_db = directory(eggnog_db)
    shell:
        """
        DB_DIR=$(dirname "{output.eggnog_db}")
        mkdir -p "$DB_DIR"
        if [ ! -d "{output.eggnog_db}" ]; then
            git clone https://github.com/eggnogdb/eggnog-mapper.git "{output.eggnog_db}"
        fi

        python "{output.eggnog_db}/download_eggnog_data.py" -y --data_dir "{output.eggnog_db}"
        """

rule emapper:
    message: "Running eggNOG-mapper"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        query = "{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa",
        eggnog_db = eggnog_db
    output:
        #annotations="{output_dir}/annotation_feat/eggnog/{sample}/{sample}.emapper.annotations"
        eggnog_output_dir="{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_emapper_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/annotation_feat/eggnog/{wildcards.sample}

        emapper.py --cpu 20 \
                   --data_dir {input.eggnog_db} \
                   -o out --override -m diamond \
                   -i {input.query} --evalue 0.001 --score 60 --pident 40 \
                   --query_cover 20 --subject_cover 20 \
                   --itype proteins --tax_scope auto \
                   --target_orthologs all \
                   --go_evidence non-electronic \
                   --pfam_realign none --report_orthologs \
                   --decorate_gff yes --excel --output_dir {output_dir}/annotation_feat/eggnog/{wildcards.sample}
        """


rule abricate:
    message: "Screening with Abricate"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        "{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gbk"
    output:
        "{output_dir}/annotation_feat/abricate/{sample}/{sample}_abricate.txt"
    params:
        datadir = abricate_db_dir,
        db = abricate_db_name
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_abricate_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/annotation_feat/abricate/{wildcards.sample}
        if [ -z "{params.datadir}" ]; then
            abricate {input} > {output}
        else
            abricate --datadir {params.datadir} --db {params.db} {input} > {output}
        fi
        """

rule map_to_ref:
    wildcard_constraints: sample="[^/\\.]+"
    message: "Mapping basecalled reads to the assembled genome, including their methylation metadata"
    input:
        bam = get_bam_input,        
        ref = get_assembly_fasta
    output:
        mapped_bam_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam",
        bai = "{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam.bai"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_mapping_to_reference_benchmark.txt"
    shell:
        """
        bash scripts/map_to_ref.sh {input.bam} {input.ref} {output.mapped_bam_file}
        """

rule install_microbemod:
    message: "Installing MicrobeMod"
    wildcard_constraints: sample="[^/\\.]+"
    output:
        touch("{output_dir}/microbemod/MicrobeMod/.installed")
    threads: threads
    shell:
        """
        mkdir -p {output_dir}/microbemod
        cd {output_dir}/microbemod
        git clone https://github.com/cultivarium/MicrobeMod.git
        cd MicrobeMod/MicrobeMod
        python download_db.py
        cd ../
        #cd MicrobeMod
        pip install .
#        pytest
        touch .installed
        """

rule microbemod_call_methylation:
    message: "MicrobeMod: Call methylation"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        bam = "{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam",
        bai = "{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam.bai",
        ref = get_assembly_fasta
    output:
        "{output_dir}/microbemod/results/{sample}/{sample}_methylated_sites.tsv",
        "{output_dir}/microbemod/results/{sample}/{sample}_motifs.tsv"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_call_methylation_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/microbemod/results/{wildcards.sample}
        
        BAM_ABS=$(realpath {input.bam})
        REF_ABS=$(realpath {input.ref})
        
        cd {output_dir}/microbemod/results/{wildcards.sample}
        
        MicrobeMod call_methylation -b "$BAM_ABS" -r "$REF_ABS" -t {threads}
        """

rule microbemod_annotate_rm:
    message: "MicrobeMod: Annotate restriction-modification systems (Auto-fix applied)"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        get_assembly_fasta
    output:
        "{output_dir}/microbemod/results/{sample}/annotate_rm/{sample}.rm.genes.tsv"
    threads: 10
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_rm_annotation_benchmark.txt"
    shell:
        """
        mkdir -p {output_dir}/microbemod/results/{wildcards.sample}/annotate_rm
        
        INPUT_ABS=$(realpath {input})
        cd {output_dir}/microbemod/results/{wildcards.sample}/annotate_rm
        MicrobeMod annotate_rm -f "$INPUT_ABS" -o {wildcards.sample} -t {threads}
        """


rule process_methylated_sites:
    message: "Create smart.gff file from microbemod output"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        ms="{output_dir}/microbemod/results/{sample}/{sample}_methylated_sites.tsv",
    output:
        smart="{output_dir}/mestudio/{sample}/{sample}_smart.gff"
    run:
        import pandas as pd
        ms_df = pd.read_table(input.ms, sep="\t", low_memory=False)
        smart = pd.DataFrame({
            'seqid': ms_df['Contig'],
            'source': 'dorado',
            'type': ms_df['Modification'],
            'start': ms_df['Position'],
            'end': ms_df['Position'],
            'score': '.',
            'strand': ms_df['Strand'],
            'phase': '.',
            'attribute': ms_df.apply(
                lambda row: f"coverage={row['Total_coverage']};context={row['Sequence']};"
                            f"Modified_bases={row['Modified_bases']};Unmodified_bases={row['Unmodified_bases']};"
                            f"Percent_modified={row['Percent_modified']}",
                axis=1
            )
        })

        smart.to_csv(output.smart, sep="\t", index=False, header=True)

rule create_motiflist:
    message: "Creating motiflist.txt from microbemod output"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        motifs="{output_dir}/microbemod/results/{sample}/{sample}_motifs.tsv"
    output:
        motiflist="{output_dir}/mestudio/{sample}/motiflist.txt"
    run:
        import pandas as pd
        motifs_df = pd.read_csv(input.motifs, sep="\t")
        motifs = motifs_df[motifs_df['Motif'] != "No Motif Assigned"]['Motif']
        motifs.to_csv(output.motiflist, index=False, header=False)


rule check_genomic_extension:
    message: "Checking genomic extension for mestudio input"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly = get_assembly_fasta
    output:
        corrected="{output_dir}/assemblies_evaluation/{sample}/{sample}_genomic.fasta"
    shell:
        """
        if [[ "{input.assembly}" != *.fasta && "{input.assembly}" != *.fna && "{input.assembly}" != *.fa ]]; then
             echo "Extension warning managed by copy." 
        fi

        cp {input.assembly} {output.corrected}
        """


rule mestudio:
    message: "Executing MeStudio"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly="{output_dir}/assemblies_evaluation/{sample}/{sample}_genomic.fasta",
        anno="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gff",
        smart="{output_dir}/mestudio/{sample}/{sample}_smart.gff",
        motiflist="{output_dir}/mestudio/{sample}/motiflist.txt"
    output:
        mestudio_results=directory("{output_dir}/mestudio/results/{sample}")
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_mestudio_benchmark.txt"
    shell:
        """
        bash scripts/imestudio -f {input.assembly} -anno {input.anno} -smart {input.smart} -mo {input.motiflist} -o {output.mestudio_results}
        """

rule circular_plots:
    message: "Generating circular density plot for each motif"
    input:
        mscore_dir="{output_dir}/mestudio/results/{sample}"
    output:
        "{output_dir}/mestudio/results/{sample}_checkout.txt"
    threads: threads
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_circular_plots_benchmark.txt"
    params:
        mscore_dir="{output_dir}/mestudio/results/{sample}/mscore"
    shell:
        """
        Rscript scripts/circular_plotter2.R {params.mscore_dir} {output}
        """

rule gsea_analysis:
    input:
        emapper="{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations",
        mscore_dir="{output_dir}/mestudio/results/{sample}"
    output:
        "{output_dir}/mestudio/results/{sample}_checkout_gsea.txt"
    benchmark:
        "{output_dir}/benchmarks/{sample}/{sample}_gsea_benchmark.txt"
    params:
        mscore_dir="{output_dir}/mestudio/results/{sample}/mscore"
    shell:
        """
        Rscript scripts/gsea2.R {params.mscore_dir} {input.emapper} {output}
        """

rule get_contig_lengths:
    message: "Generating contig length file from assembly."
    output:
        contig_lengths="{output_dir}/assemblies_evaluation/{sample}/contig_lengths.tsv"
    run:
        from Bio import SeqIO
        import pandas as pd
        
        if assembler_choice == "dragonflye":
            fa = f"{output_dir}/assemblies_evaluation/dragonflye/{wildcards.sample}/{wildcards.sample}.fa"
        elif assembler_choice == "canu":
            fa = f"{output_dir}/assemblies_evaluation/canu/{wildcards.sample}/{wildcards.sample}.fa"
        elif assembler_choice == "raven":
            fa = f"{output_dir}/assemblies_evaluation/raven/{wildcards.sample}/{wildcards.sample}.fa"
            
        lengths = []
        for record in SeqIO.parse(fa, "fasta"):
            lengths.append({"contig": record.id, "length": len(record.seq)})
            
        pd.DataFrame(lengths).to_csv(output.contig_lengths, sep="\t", index=False)


for control, treated in dmr_pairs:
    rule dmr_analysis:
        message: f"DMR analysis: comparing {control} vs {treated}"
        input:
            control="{output_dir}/microbemod/results/{control}/{control}_methylated_sites.tsv",
            treated="{output_dir}/microbemod/results/{treated}/{treated}_methylated_sites.tsv",
            contig="{output_dir}/assemblies_evaluation/{control}/contig_lengths.tsv"
        output:
            volcano="{output_dir}/DMRs/{control}_vs_{treated}_volcano.png",
            windows="{output_dir}/DMRs/{control}_vs_{treated}_window_results.tsv"
        params:
            window=dmr_window
        benchmark:
            "{output_dir}/benchmarks/{sample}/{sample}_dmr_benchmark.txt"
        shell:
            """
            mkdir -p {output_dir}/DMRs
            python scripts/dmr.py \
                --control {input.control} \
                --treated {input.treated} \
                --contig_lengths {input.contig} \
                --window {params.window} \
                --output {output.volcano} \
                --results {output.windows}
            """

rule global_benchmark_report:
    message: "Generating global benchmark summary for all samples."
    input:
        expand("{output_dir}/mestudio/results/{sample}_checkout_gsea.txt", output_dir=output_dir, sample=samples)
    output:
        raw = "{output_dir}/benchmarks/global_benchmark_data.tsv",
        stats = "{output_dir}/benchmarks/global_benchmark_stats_summary.tsv"
    params:
        bench_dir = directory("{output_dir}/benchmarks")
    script:
        "scripts/summary_benchmark.py"
