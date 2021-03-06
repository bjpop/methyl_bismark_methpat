'''Run the pipeline as a program.'''

from pipeline import make_pipeline
import pipeline_base.main
import pkg_resources 

PROGRAM_NAME = "methylation_pipeline" 
PROGRAM_INFO = pkg_resources.require(PROGRAM_NAME)[0]
PROGRAM_VERSION = PROGRAM_INFO.version

def main():
    pipeline_base.main.main(PROGRAM_NAME, PROGRAM_VERSION, make_pipeline)
