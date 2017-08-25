'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import PipelineStages 

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='methylation_pipeline')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = PipelineStages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Run bismark genome preparation on the reference genome 
    pipeline.originate(
        task_func=stages.bismark_genome_prepare,
        name='bismark_genome_prepare',
        output='reference/Bisulfite_Genome')

    # Run FASTQC on the input fastq files 
    pipeline.transform(
        task_func=stages.fastqc,
        name='fastqc',
        input=output_from('original_fastqs'),
        filter=formatter('(?P<path>.+)/(?P<filename>.+).fastq.gz'),
        output='{path[0]}/{filename[0]}_fastqc')

    # Run bismark on the input fastq files 
    (pipeline.transform(
        task_func=stages.bismark,
        name='bismark',
        input=output_from('original_fastqs'),
        filter=formatter('(?P<path>.+)/(?P<filename>.+)_R1_(?P<num>.+).fastq.gz'),
        add_inputs=add_inputs('{path[0]}/{filename[0]}_R2_{num[0]}.fastq.gz'),
        extras=['{path[0]}/bismark_output/'],
        output='{path[0]}/bismark_output/{filename[0]}_R1_{num[0]}_bismark_bt2_pe.sam.gz')).follows('bismark_genome_prepare')

    # Run bismark methylation extractor on the bismark output 
    pipeline.transform(
        task_func=stages.bismark_methylation_extractor,
        name='bismark_methylation_extractor',
        input=output_from('bismark'),
        filter=formatter('(?P<path>.+)/(?P<filename>.+)_R1_(?P<num>.+)_bismark_bt2_pe.sam.gz'),
        extras=['{path[0]}'],
        output='{path[0]}/CpG_context_{filename[0]}_R1_{num[0]}_bismark_bt2_pe.sam.gz.txt')

    # Run methpt on the bismark methylation extractor output 
    pipeline.transform(
        task_func=stages.methpat,
        name='methpat',
        input=output_from('bismark_methylation_extractor'),
        filter=formatter('(?P<path>.+)/CpG_context_(?P<filename>.+)'),
        extras=['{path[0]}', '{filename[0]}'],
        output='{path[0]}/CpG_context_{filename[0]}.methpat.html')


    return pipeline

'''
    # Index the human reference genome with BWA, needed before we can map reads
    pipeline.transform(
        task_func=stages.index_ref_bwa,
        name='index_ref_bwa',
        input=output_from('human_reference_genome'),
        filter=suffix('.fa'),
        output=['.fa.amb', '.fa.ann', '.fa.pac', '.fa.sa', '.fa.bwt'])

    # Align paired end reads in FASTQ to the reference producing a BAM file
    (pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name. 
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        filter=formatter('.+/(?P<sample>[_a-zA-Z0-9]+)_R1.fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        add_inputs=add_inputs('{path[0]}/{sample[0]}_R2.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='{path[0]}/{sample[0]}.bam')
        .follows('index_ref_bwa'))
'''
