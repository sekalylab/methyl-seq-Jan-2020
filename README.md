# methylSeq data processing pipeline

master shell script:      script-methyl-seq.sh  
Slurm batch running script:     batch-methyl-seq.slurm  
  
one methyl-seq bam file with ~50GB will take ~12 hours  
  
## prerequisite
g++ -o sam-overlap-rp sam-overlap-rp.C  
g++ -o methyl-seq-level.strandOnly methyl-seq-level.strandOnly.C  