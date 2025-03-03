import os
import glob

#Directories and parameters
input_dir = config["folder"]["input_dir"]
output_dir = config["folder"]["output_dir"]
sample_tab = config["params"]["sample_tab"]
samples = config["params"]["sample"]
dorado_models = config["params"]["dorado_models"]
filtlong_min_length=config["params"]["filtlong"]["min_length"]
filtlong_keep_percent=config["params"]["filtlong"]["keep_percent"]
racon_rounds=config["params"]["dragonflye"]["racon_rounds"]
medaka_rounds=config["params"]["dragonflye"]["medaka_rounds"]
medaka_model=config["params"]["dragonflye"]["model"]
db_dir = config["folder"]["db_dir"]
eggnog_db = config["params"]["eggnog"]["db"]
#db_dir=config.get("db_dir", "databases_directory") 
threads = config.get("threads", os.cpu_count())
print(samples)
abricate_db_dir = config.get("abricate_db_dir", "")
abricate_db_name = config.get("abricate_db_name", "")

def sampleInfos(sample_config, basecalling_dir):
    fast5_dict = {}
    pod5_dict = {}
    with open(sample_config, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            sample, index, file, filetype = line.strip().split()
            if filetype == "fast5":
                fast5_dict.setdefault(sample, []).append(file)
            elif filetype == "pod5":
                pod5_dict.setdefault(sample, []).append(file)
    for sample in fast5_dict:
        if sample not in pod5_dict:
            sample_pod5 = os.path.join(basecalling_dir, sample, f"{sample}_output.pod5")
            pod5_dict[sample] = [sample_pod5]
    return fast5_dict, pod5_dict

basecalling_dir = os.path.join(output_dir, "basecalling")
sample_d_fast5, sample_d_pod5 = sampleInfos(sample_tab, basecalling_dir)
#eggnog_db="db_dir/eggnog_db"
#print(sample_d_fast5)
#print(sample_d_pod5)

def getReadFiles(wildcards, d=sample_d_fast5):
    #print(d[wildcards.sample])
    return d[wildcards.sample]

def getReadFilesDict(wildcards, d=sample_d_fast5):
    out = {wildcards.sample: d[wildcards.sample]}
    #print(out)
    return out

def getPodFilesDict(wildcards, d=sample_d_pod5):
    out = {wildcards.sample: d[wildcards.sample]}
    print(out)
    return out

rule all:
    input:
        #Basecalling
        expand("{output_dir}/basecalling/{sample}/{sample}_output.pod5", output_dir=output_dir, sample=sample_d_fast5.keys()),
        expand("{output_dir}/dorado_models/{model}",output_dir=output_dir, model=dorado_models),
        expand("{output_dir}/basecalling/{sample}/{sample}_basecalling.bam",output_dir=output_dir, sample=sample_d_pod5.keys()),
        expand("{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq", output_dir=output_dir, sample=samples),

        #Quality check and filtering
        expand("{output_dir}/qc_filtering/nanoplot_b/{sample}_nanoplot_qc_report_b", output_dir=output_dir, sample=samples),
        expand("{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq", output_dir=output_dir, sample=samples),
        expand("{output_dir}/qc_filtering/nanoplot_a/{sample}_nanoplot_qc_report_a", output_dir=output_dir, sample=samples),

        #Assembly and evaluation
        expand("{output_dir}/assemblies_eva/{sample}/{sample}.fa", output_dir=output_dir, sample=samples),
        expand("{output_dir}/assemblies_eva/quast/{sample}", output_dir=output_dir, sample=samples),

        #Annotation
        expand("{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa", output_dir=output_dir, sample=samples),
        expand("{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations", output_dir=output_dir, sample=samples),
        #expand("{output_dir}/annotation_feat/abricate/{sample}_abricate.txt", output_dir=output_dir, sample=samples) 
        #if config.get("abricate_db_dir", "") and config.get("abricate_db_name", "") else [],
        f"{eggnog_db}",

        #Microbemod
        expand("{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam", output_dir=output_dir, sample=samples),
        expand("{output_dir}/microbemod/MicrobeMod/.installed", output_dir=output_dir),
        expand("{output_dir}/microbemod/results/{sample}/{sample}_methylated_sites.tsv", output_dir=output_dir, sample=samples),

        #MeStudio
        expand("{output_dir}/mestudio/{sample}/{sample}_smart.gff", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/{sample}/motiflist.txt", output_dir=output_dir, sample=samples),
        expand("{output_dir}/assemblies_eva/{sample}/{sample}_genomic.fasta", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/results/{sample}", output_dir=output_dir, sample=samples),
        expand("{output_dir}/mestudio/results/{sample}_checkout.txt", output_dir=output_dir, sample=samples),

        #GSEA
        expand("{output_dir}/mestudio/results/{sample}_checkout_gsea.txt", output_dir=output_dir, sample=samples)

rule convert_fast5_to_pod5:
    message: "Converting .fast5 files to POD5 file before Dorado basecalling."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        unpack(getReadFilesDict)
    output:
        "{output_dir}/basecalling/{sample}/{sample}_output.pod5"
    threads: threads
    conda:
        "envs/basecalling.yaml"
    shell:
        """
        mkdir -p {output_dir}/basecalling/{wildcards.sample}
        pod5 convert fast5 {input} --output {output} > /dev/null 2>&1
        """

rule download_dorado_model:
    message: "Download specific basecalling models"
    output:
        model = directory("{output_dir}/dorado_models/{model}")
    threads: threads
    params:
        model="{model}"
    conda:
        "envs/basecalling.yaml"
    priority: 50
    shell:
        """
        mkdir -p {output_dir}/dorado_models
        dorado download --model {params.model} --directory {output_dir}/dorado_models > /dev/null 2>&1
        """

rule dorado_basecalling:
    message: "Dorado basecalling"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        unpack(getPodFilesDict)
    output:
        "{output_dir}/basecalling/{sample}/{sample}_basecalling.bam"
    threads: threads
    params:
        models=expand(f"{output_dir}/dorado_models/{{model}}", model=dorado_models)
    conda:
        "envs/basecalling.yaml"
    priority: 45
    shell:
        """
        dorado basecaller {params.models[0]} {input} --modified-bases-models {params.models[1]},{params.models[2]} > {output}
        """

rule convert_bam_to_fastq:
    message: "Convert basecalling.bam file to fastq"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        bam_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.bam"
    output:
        fastq_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    threads: threads
    conda:
        "envs/basecalling.yaml"
    priority: 40
    shell:
        """
        samtools fastq {input.bam_file} > {output.fastq_file} 2> /dev/null
        """

rule nanoplot_qc_b:
    message: "Quality check: executing NanoPlot to generate summaries and graphical representations of data statistics."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        fastq_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    output:
        qc_report=directory("{output_dir}/qc_filtering/nanoplot_b/{sample}_nanoplot_qc_report_b")
    threads: threads
    conda:
        "envs/qc_filtering.yaml"
    priority: 35
    shell:
        """
        mkdir -p {output.qc_report}
        NanoPlot --fastq {input.fastq_file} --N50 -o {output.qc_report} > /dev/null 2>&1
        """

rule filtlong_filter:
    message: "Filtering: executing Filtlong to remove low-quality reads."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        fastq_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.fastq"
    output:
        filtered_fastq="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    threads: threads
    params:
        min_length=filtlong_min_length,
        keep_percent=filtlong_keep_percent
    conda:
        "envs/qc_filtering.yaml"
    priority: 30
    shell:
        """
        filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input.fastq_file} > {output.filtered_fastq}
        """

rule nanoplot_qc_a:
    message: "Quality check after filtering: executing NanoPlot."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        filtered_fastq="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        qc_report_a=directory("{output_dir}/qc_filtering/nanoplot_a/{sample}_nanoplot_qc_report_a")
    threads: threads
    conda:
        "envs/qc_filtering.yaml"
    priority: 25
    shell:
        """
        mkdir -p {output.qc_report_a}
        NanoPlot --fastq {input.filtered_fastq} --N50 -o {output.qc_report_a} > /dev/null 2>&1
        """

rule dragonflye_assembly:
    message: "De novo genome assembly with Flye-dragonflye assembler and polishing with Racon and Medaka."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        reads="{output_dir}/qc_filtering/filtered/{sample}/{sample}_basecalling_filtered.fastq"
    output:
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}.fa"
    threads: threads
    params:
        racon_rounds=racon_rounds,
        medaka_rounds=medaka_rounds,
        medaka_model=medaka_model
        #gsize=config["params"]["dragonflye"]["gsize"],
        #prefix="{sample}"
    conda:
        "envs/assemblies_eva.yaml"
    priority: 20
    shell:
        """
        dragonflye \
        --outdir {output_dir}/assemblies_eva/{wildcards.sample} \
        --reads {input.reads} \
        --assembler flye \
        --nanohq \
        --racon {params.racon_rounds} \
        --medaka {params.medaka_rounds} \
        --model {params.medaka_model} \
        --prefix {wildcards.sample} \
        --dnaapler_mode all --force \
        > /dev/null 2>&1
        """

rule quast:
    message: "Quality assessment of the genome assembly using QUAST"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}.fa"
    output:
        quast_report_dir=directory("{output_dir}/assemblies_eva/quast/{sample}")
    threads: threads
    conda:
        "envs/assemblies_eva.yaml"
    shell:
        """
        quast.py {input.assembly} -o {output.quast_report_dir} > /dev/null 2>&1
        """

rule prokka:
    message: "Annotate genome with Prokka"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}.fa"
    output:
        faa="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa",
        gbk="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gbk",
        gff="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gff"
    threads: threads
    params:
        prefix="{sample}_anno"
    conda:
        "envs/annotation_feat.yaml"
    shell:
        """
        prokka --outdir {output_dir}/annotation_feat/prokka/{wildcards.sample} \
               --prefix {params.prefix} \
               --force \
               {input.assembly} > /dev/null 2>&1
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
    message: "Running eggNOG-mapper for functional annotation"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        query = "{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.faa",
        eggnog_db = eggnog_db
    output:
        eggnog_output_dir="{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations"
    threads: threads
    conda:
        "envs/annotation_feat.yaml"
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
    message: "Screening of contigs for antimicrobial resistance or virulence genes with abricate."
    wildcard_constraints: sample="[^/\\.]+"
    input:
        "{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gbk"
    output:
        "{output_dir}/annotation_feat/abricate/{sample}_abricate.txt"
    params:
        datadir = config.get("abricate_db_dir", ""),
        db = config.get("abricate_db_name", "")
    threads: threads
    conda:
        "envs/annotation_feat.yaml"
    shell:
        """
        abricate --nopath --datadir {params.datadir} --db {params.db} {input} > {output} 2> /dev/null
        """

rule map_to_ref:
    wildcard_constraints: sample="[^/\\.]+"
    message: "Mapping basecalled reads to the assembled genome, including their methylation metadata"
    input:
        bam_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.bam",
        ref="{output_dir}/assemblies_eva/{sample}/{sample}.fa"
    output:
        mapped_bam_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam"
    threads: threads
    conda:
        "envs/microbemod.yaml"
    shell:
        """
        bash scripts/map_to_ref.sh {input.bam_file} {input.ref} {output.mapped_bam_file} > /dev/null 2>&1
        """

rule install_microbemod:
    message: "Installing MicrobeMod"
    wildcard_constraints: sample="[^/\\.]+"
    output:
        touch("{output_dir}/microbemod/MicrobeMod/.installed")
    threads: threads
    conda:
        "envs/microbemod.yaml"
    shell:
        """
        mkdir -p {output_dir}/microbemod
        cd {output_dir}/microbemod
        git clone https://github.com/cultivarium/MicrobeMod.git
        cd MicrobeMod/MicrobeMod
        python download_db.py
        cd ../
        pip install .
        pytest
        touch {output}
        """

rule microbemod_call_methylation:
    message: "MicrobeMod call methylation"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        mapped_bam_file="{output_dir}/basecalling/{sample}/{sample}_basecalling.mapped.bam",
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}.fa",
    output:
        microbemod_ms="{output_dir}/microbemod/results/{sample}/{sample}_methylated_sites.tsv",
        microbemod_motifs="{output_dir}/microbemod/results/{sample}/{sample}_motifs.tsv"
    threads: threads
    conda:
        "envs/microbemod.yaml"
    shell:
        """
        mkdir -p {output_dir}/microbemod/results/{wildcards.sample}
        cd {output_dir}/microbemod/results/{wildcards.sample}
        MicrobeMod call_methylation -b {input.mapped_bam_file} -r {input.assembly} -t 10 > /dev/null 2>&1
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
        motifs.to_csv(output.motiflist, index=False, header=False, sep="\n")

rule check_genomic_extension:
    message: "Checking genomic extension for mestudio input"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}.fa"
    output:
        corrected="{output_dir}/assemblies_eva/{sample}/{sample}_genomic.fasta"
    shell:
        """
        if {input.assembly} != *.fasta && {input.assembly} != *.fna; then
            echo "Error: The genomic file {input.assembly} does not have a valid extension. Please add .fasta or .fna." >&2
            exit 1
        fi

        cp {input.assembly} {output.corrected}
        """

rule mestudio:
    message: "Executing MeStudio"
    wildcard_constraints: sample="[^/\\.]+"
    input:
        assembly="{output_dir}/assemblies_eva/{sample}/{sample}_genomic.fasta",
        anno="{output_dir}/annotation_feat/prokka/{sample}/{sample}_anno.gff",
        smart="{output_dir}/mestudio/{sample}/{sample}_smart.gff",
        motiflist="{output_dir}/mestudio/{sample}/motiflist.txt"
    output:
        mestudio_results=directory("{output_dir}/mestudio/results/{sample}")
    threads: threads
    conda:
        "envs/mestudio.yaml"
    shell:
        """
        bash scripts/imestudio -f {input.assembly} -anno {input.anno} -smart {input.smart} -mo {input.motiflist} -o {output.mestudio_results}
        """

rule circular_plots:
    message: "Generating circular density plot for each motif"
    output:
        "{output_dir}/mestudio/results/{sample}_checkout.txt"
    threads: threads
    conda:
        "envs/mestudio.yaml"
    params:
        mscore_dir="{output_dir}/mestudio/results/{sample}/mscore"
    shell:
        """
        Rscript scripts/circular_plotter.R {params.mscore_dir} {output} > /dev/null 2>&1
        """

rule gsea_analysis:
    input:
        "{output_dir}/annotation_feat/eggnog/{sample}/out.emapper.annotations"
    output:
        "{output_dir}/mestudio/results/{sample}_checkout_gsea.txt"
    conda:
        "envs/gsea.yaml"
    params:
        mscore_dir="{output_dir}/mestudio/results/{sample}/mscore"
    shell:
        """
        Rscript scripts/gsea.R {params.mscore_dir} {input} {output}
        """
