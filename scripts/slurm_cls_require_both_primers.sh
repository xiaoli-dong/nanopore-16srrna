#!/bin/bash

#SBATCH --job-name=cls_16srRNA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10-00:00:00
#SBATCH --mem=16G
#SBATCH --partition=vm-cpu
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
#SBATCH --mail-user=xiaoli.dong@albertaprecisionlabs.ca
#SBATCH --mail-type=ALL

# ---------------------------------------------------------------------
# this workflow is developed by xiaoli dong to process 16s rRNA data 
# generated using nanopore technology from provlab south


progs="/nfs/Genomics_DEV/projects/xdong/deve/16srRNA/nanopore/progs"
hostile_minimap2_db="/nfs/APL_Genomics/db/prod/hostile/minimap2_ref/human-t2t-hla-argos985-mycob140.fa.gz"
emu_default_db="${progs}/emu_db/emu_default"
emu_silva_db="${progs}/emu_db/silva"
emu_rdp_db="${progs}/emu_db/rdp"

echo "Current working directory: $(pwd)"
echo "Starting run at: $(date)"
# ---------------------------------------------------------------------

workdir="/nfs/Genomics_DEV/projects/xdong/deve/16srRNA/nanopore/230809_S_N_322/analysis"
results_dir="${workdir}/results_cls"
runid="230809_S_N_322"
minlen=400
maxlen=600

echo ${workdir}
# ---------------------------------------------------------------------
echo "${results_dir}"
[ ! -d "${results_dir}" ] && mkdir -p ${results_dir}

for var in {1..24} {49..72}; do

    prefix="barcode$var"

    if [ "$var" -le 9 ]; then
        prefix="barcode0$var"
    fi

    # raw stats
    [ ! -d "${results_dir}/${prefix}/raw" ] && mkdir -p ${results_dir}/$prefix/raw
    cmd="singularity run ${progs}/seqkit:2.5.1--h9ee0642_0 seqkit stats -a -T -j $SLURM_CPUS_PER_TASK ${workdir}/fastq/${prefix}/${prefix}.fastq.gz --out-file ${results_dir}/${prefix}/raw/${prefix}.seqkit_stats.tsv"
    echo $cmd
    #$cmd

    ###################################### qc ##################
    # raw reads quality assesment

    cmd="singularity run ${progs}/nanoplot:1.41.6--pyhdfd78af_0 NanoPlot --threads $SLURM_CPUS_PER_TASK --fastq ${workdir}/fastq/${prefix}/${prefix}.fastq.gz --outdir ${results_dir}/${prefix}/qc/nanoplot_raw -f json"
    echo $cmd
    #$cmd

    # chop adapters
    [ ! -d "${results_dir}/${prefix}/qc/porechop" ] && mkdir -p ${results_dir}/$prefix/qc/porechop

    cmd="singularity run ${progs}/porechop:0.2.4--py39h1f90b4d_6 porechop --threads $SLURM_CPUS_PER_TASK --input ${workdir}/fastq/${prefix}/${prefix}.fastq.gz --output ${results_dir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz"
    echo $cmd
    echo ""
    #$cmd

    cmd="singularity run ${progs}/seqkit:2.5.1--h9ee0642_0 seqkit stats -a -T -j $SLURM_CPUS_PER_TASK ${results_dir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz --out-file ${results_dir}/${prefix}/qc/porechop/${prefix}.porechop.seqkit_stats.tsv"
    echo $cmd
    #$cmd

    #################### CLS primer
    #fwd: AGA +GTT T+GA TCM TGG CTC AII III AAC GCT
    #Rev CGC GGC TGC TGG CAI IIA ITT RGC

    #Fwd-5'-> 3' - AGAGTTTGATCMTGGCTCANNNNNAACGCT reverse complemnts: AGCGTTNNNNNTGAGCCAKGATCAAACTCT
    #Rev-5'-> 3' - CGC GGC TGC TGG CANNNANTTRGC reverse complemnts: GCYAANTNNNTGCCAGCAGCCGCG
    #get rid of primers
    #in cutadapt, The character I, used to encode the base inosine, is automatically replaced with N within the adapter sequence.

    [ ! -d "${results_dir}/${prefix}/qc/cutadapt" ] && mkdir -p ${results_dir}/$prefix/qc/cutadapt

    cmd="singularity run  ${progs}/cutadapt:4.4--py39hf95cd2a_1 cutadapt -j $SLURM_CPUS_PER_TASK -e 0.2 --overlap 10 --revcomp  -m ${minlen} -M ${maxlen} --max-n 1  -g "AGAGTTTGATCMTGGCTCAIIIIIAACGCT...GCYAAITIIITGCCAGCAGCCGCG" --untrimmed-output=${results_dir}/${prefix}/qc/cutadapt/${prefix}.cutadapt_untrimmed.fastq.gz -o ${results_dir}/${prefix}/qc/cutadapt/${prefix}.cutadapt_trimmed.fastq.gz ${results_dir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz"
    echo $cmd
    echo ""
    #$cmd

    cmd="singularity run ${progs}/seqkit:2.5.1--h9ee0642_0 seqkit stats -a -T -j $SLURM_CPUS_PER_TASK --out-file ${results_dir}/${prefix}/qc/cutadapt/${prefix}.cutadapt_trimmed.seqkit_stats.tsv ${results_dir}/${prefix}/qc/cutadapt/${prefix}.cutadapt_trimmed.fastq.gz"
    echo $cmd
    #$cmd

    # quality filtering
    [ ! -d "${results_dir}/$prefix/qc/chopper" ] && mkdir -p ${results_dir}/$prefix/qc/chopper
    cmd="zcat ${results_dir}/${prefix}/qc/cutadapt/${prefix}.cutadapt_trimmed.fastq.gz \| singularity run ${progs}/chopper:0.6.0--hdcf5f25_0 chopper -q 10 --minlength ${minlen} --maxlength ${maxlen} > ${results_dir}/${prefix}/qc/chopper/${prefix}.chopper.fastq"
    echo $cmd
    echo ""
    #(zcat ${results_dir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz | singularity run ${progs}/chopper:0.6.0--hdcf5f25_0 chopper -q 10 --minlength ${minlen} --maxlength ${maxlen} >${results_dir}/${prefix}/qc/chopper/${prefix}.chopper.fastq)

    cmd="singularity run ${progs}/seqkit:2.5.1--h9ee0642_0 seqkit stats -a -T -j $SLURM_CPUS_PER_TASK --out-file ${results_dir}/${prefix}/qc/chopper/${prefix}.chopper.seqkit_stats.tsv ${results_dir}/${prefix}/qc/chopper/${prefix}.chopper.fastq"
    echo $cmd
    #$cmd

    #dehost
    [ ! -d "${results_dir}/$prefix/qc/dehost" ] && mkdir -p ${results_dir}/$prefix/qc/dehost
    cmd="singularity run ${progs}/hostile:0.1.0--pyhdfd78af_0 hostile clean --threads $SLURM_CPUS_PER_TASK --fastq1 ${results_dir}/${prefix}/qc/chopper/${prefix}.chopper.fastq --aligner minimap2 --index ${hostile_minimap2_db} --out-dir ${results_dir}/${prefix}/qc/dehost"
    echo $cmd
    #$cmd

    cmd="singularity run ${progs}/seqkit:2.5.1--h9ee0642_0 seqkit stats -a -T -j $SLURM_CPUS_PER_TASK --out-file ${results_dir}/${prefix}/qc/dehost/${prefix}.chopper.clean.seqkit_stats.tsv ${results_dir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz"
    echo $cmd
    #$cmd

    cmd="singularity run ${progs}/csvtk:0.28.0--h9ee0642_0 csvtk concat --out-file ${results_dir}/${prefix}/${prefix}_seqkit_stats.tsv ${results_dir}/${prefix}/raw/*.seqkit_stats.tsv ${results_dir}/${prefix}/qc/porechop/*.seqkit_stats.tsv ${results_dir}/${prefix}/qc/cutadapt/*.seqkit_stats.tsv ${results_dir}/${prefix}/qc/chopper/*.seqkit_stats.tsv ${results_dir}/${prefix}/qc/dehost/*.seqkit_stats.tsv"
    echo $cmd
    echo ""
    $cmd
    ################ classify ################
    cmd="singularity run ${progs}/emu:3.4.5--hdfd78af_0 emu abundance --threads $SLURM_CPUS_PER_TASK --db ${emu_default_db} --output-dir ${results_dir}/${prefix}/abundance/emu_default --keep-files --keep-counts --keep-read-assignments --output-unclassified ${results_dir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz"
    echo $cmd
    #$cmd

    cmd="singularity run ${progs}/emu:3.4.5--hdfd78af_0 emu abundance --threads $SLURM_CPUS_PER_TASK --db ${emu_silva_db} --output-dir ${results_dir}/${prefix}/abundance/emu_silva --keep-files --keep-counts --keep-read-assignments --output-unclassified ${results_dir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz"
    echo $cmd
    #$cmd

    cmd="singularity run ${progs}/emu:3.4.5--hdfd78af_0 emu abundance --threads $SLURM_CPUS_PER_TASK --db ${emu_rdp_db} --output-dir ${results_dir}/${prefix}/abundance/emu_rdp --keep-files --keep-counts --keep-read-assignments --output-unclassified ${results_dir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz"
    echo $cmd
    echo ""
    #$cmd

done

cmd="singularity run ${progs}/csvtk:0.28.0--h9ee0642_0 csvtk concat ${results_dir}/*/raw/*.seqkit_stats.tsv --out-file ${results_dir}/${runid}-raw.seqkit_stats.tsv"
echo $cmd
echo ""
#$cmd

for var in porechop cutadapt chopper dehost; do
    cmd="singularity run ${progs}/csvtk:0.28.0--h9ee0642_0 csvtk concat ${results_dir}/*/qc/$var/*.seqkit_stats.tsv --out-file ${results_dir}/${runid}-$var.seqkit_stats.tsv"
    echo $cmd
    echo ""
    #$cmd
done

# silva has no ranking infor
for var in emu_default emu_rdp; do
    [ ! -d "${results_dir}/$var" ] && mkdir -p "${results_dir}/$var"

    cp ${results_dir}/barcode*/abundance/$var/*.fastq_rel-abundance.tsv ${results_dir}/$var/
    for rank in species genus family order class phylum; do
        singularity run ${progs}/emu:3.4.5--hdfd78af_0 emu combine-outputs --split-tables --counts ${results_dir}/$var $rank
    done
done
# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: $(date)"
