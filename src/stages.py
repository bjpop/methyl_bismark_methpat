'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

#from pipeline_base.utils import safe_make_dir, run_java
from pipeline_base.utils import safe_make_dir
from pipeline_base.stages import Stages
from pipeline_base.runner import run_stage
import os.path

class PipelineStages(Stages):
    def __init__(self, *args, **kwargs):
        super(PipelineStages, self).__init__(*args, **kwargs)
        #self.reference = self.get_options('human_reference_genome')


    def original_fastqs(self, output):
        '''Original fastq files'''
        pass

    def human_reference_genome(self, output):
        '''Human reference genome in FASTA format'''
        pass

    def bismark_genome_prepare(self, bisulfite_genome):
        '''Prepare the human genome using bismark'''
        command = "/data/projects/punim0095/methylation/bismark_v0.18.1/bismark_genome_preparation --path_to_bowtie /usr/local/easybuild/software/Bowtie2/2.2.9-intel-2016.u3/bin/ --verbose reference/"
        run_stage(self.state, 'bismark_genome_prepare', command)


    def fastqc(self, fastq_in, dir_out):
        '''Quality check fastq file using fastqc'''
        safe_make_dir(dir_out)
        command = 'fastqc --extract -o {dir} {fastq}'.format(dir=dir_out, fastq=fastq_in)
        run_stage(self.state, 'fastqc', command)

    def bismark(self, inputs, file_out, dir_out):
        '''Run bismark on the fastq files'''
        fastq_read1_in, fastq_read2_in = inputs
        command = '/data/projects/punim0095/methylation/bismark_v0.18.1/bismark -N 0 -q -o {dir_out} --ambiguous --non_directional /data/projects/punim0095/methylation/reference/ --path_to_bowtie /usr/local/easybuild/software/Bowtie2/2.2.9-intel-2016.u3/bin/ -1 {fastq_read1_in} -2 {fastq_read2_in}'.format(dir_out=dir_out, fastq_read1_in=fastq_read1_in, fastq_read2_in=fastq_read2_in)
        self.run('bismark', command)

    def bismark_methylation_extractor(self, bismark_file, file_out, dir_out):
        '''Run bismark methylation extractor on output of bismark'''
        command = '/data/projects/punim0095/methylation/bismark_v0.18.1/bismark_methylation_extractor --comprehensive -o {dir_out} {bismark_file}'.format(dir_out=dir_out, bismark_file=bismark_file) 
        self.run('bismark_methylation_extractor', command)

    def methpat(self, bismark_file, file_out, path_out, suffix_out):
        '''Run methpat on the output of the bismark methylation extractor'''
        CpG_bismark_file = os.path.join(path_out, "CpG_context_" + suffix_out)
        CpG_html = os.path.join(path_out, "CpG_context_" + suffix_out + ".html")
        CpG_methpat_out = os.path.join(path_out, "CpG_context_" + suffix_out + ".methpat.tsv")
        command = 'methpat --html {html_out} --amplicons /data/projects/punim0095/methylation/amplicons.txt {bismark_file} > {methpat_out}'.format(html_out=CpG_html, bismark_file=CpG_bismark_file, methpat_out=CpG_methpat_out) 
        self.run('methpat', command)

#    def index_ref_bwa(self, reference_in, index_outputs):
#        '''Index human genome reference with BWA'''
#        command = "bwa index -a bwtsw {reference}".format(reference=reference_in)
#        self.run('index_ref_bwa', command)
#
#
#    def align_bwa(self, inputs, bam_out, sample_id):
#        '''Align the paired end fastq files to the reference genome using bwa'''
#        fastq_read1_in, fastq_read2_in = inputs
#        cores = self.get_stage_options('align_bwa', 'cores')
#        read_group = '"@RG\\tID:{sample}\\tSM:{sample}\\tPL:Illumina"'.format(sample=sample_id)
#        command = 'bwa mem -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
#                  '| samtools view -b -h -o {bam} -' \
#                  .format(cores=cores,
#                      read_group=read_group,
#                      fastq_read1=fastq_read1_in,
#                      fastq_read2=fastq_read2_in,
#                      reference=self.reference,
#                      bam=bam_out)
#        self.run('align_bwa', command)
