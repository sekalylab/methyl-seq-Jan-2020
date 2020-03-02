#!/bin/bash
# @author Linyong Mao

bamFile="$1" # bf
sampleId=$(echo $bamFile | \
	       sed -r 's|.+/||' | \
	       sed -r 's/.bam//' | \
	       sed 's/ /_/g') # pref

# extract reads properly paired
properpairsFile=$(echo $bamFile | sed -r 's/.bam/.properpairs.sam/g')
samtools view $bamFile | \
    awk 'BEGIN {FS = "\t";} $2 == 99 || $2 == 147 || $2 == 83 || $2 == 163' > \
	$properpairsFile
wc -l $properpairsFile

### frag > 0, then read upsteam its mate; frag < 0, then read downsteam its mate
forwardReadUpstreamFile=$(echo $properpairsFile | \
			      sed -r 's/.properpairs.sam/.fowReadUpstream.sam/g')
awk 'BEGIN {FS = "\t";} ($9 > 0 && $4 <= $8) || ($9 < 0 && $4 >= $8)' $properpairsFile > \
    $forwardReadUpstreamFile
wc -l  $forwardReadUpstreamFile
# clean up
rm $properpairsFile

### overlapped portion of read pair 
overlapFile=$(echo $forwardReadUpstreamFile | \
		  sed -r 's/.fowReadUpstream.sam/.overlap.sam/g')
overlapOutput=$(echo $forwardReadUpstreamFile | \
		    sed -r 's/.fowReadUpstream.sam/.overlapOutput.sam/g')
sam-overlap-rp $forwardReadUpstreamFile $overlapFile > $overlapOutput
wc -l $overlapOutput
# clean up
rm $forwardReadUpstreamFile
rm $overlapOutput

# split Watson and Crick reads
watsonFile=$(echo $overlapFile | \
                    sed -r 's/.overlap.sam/.watson.sam/g')
crickFile=$(echo $overlapFile | \
			 sed -r 's/.overlap.sam/.crick.sam/g')
awk 'BEGIN {FS = "\t";} $2 == 99 || $2 == 147' $overlapFile > $watsonFile
awk 'BEGIN {FS = "\t";} $2 == 83 || $2 == 163' $overlapFile > $crickFile

wc -l $watsonFile
wc -l $crickFile
# clean up
rm $overlapFile

# convert to BAM
watsonBamFile=$(echo $watsonFile | \
                    sed -r 's/.sam/.bam/g')
samtools view \
	 -b \
	 -S \
	 -t /home/lxm416/hg19/hg19.fa.fai \
	 -o $watsonBamFile \
	 $watsonFile
rm $watsonFile
watsonSortedFile=$(echo $watsonBamFile | \
		       sed -r 's/.bam/.sorted.bam/g')

samtools sort -o $watsonSortedFile $watsonBamFile
rm $watsonBamFile

# pileup
watsonPileFile=$(echo $watsonSortedFile | \
                                sed -r 's/.sorted.bam/.pile.txt/g')
samtools mpileup \
	 -f /home/lxm416/hg19/hg19.fa \
	 -o $watsonPileFile \
	 $watsonSortedFile
rm $watsonSortedFile
watsonCPGfile=$(echo $watsonPileFile | \
		    sed -r 's/.pile.txt/.cpg.txt/g')

methyl-seq-level.strandOnly \
    /home/lxm416/hg19/hg19.CpG.sites.list \
    $watsonCPGfile \
    $watsonPileFile \
    $sampleId \
    5 \
    reverse
rm $watsonPileFile

# convert to BAM
crickBamFile=$(echo $crickFile | \
                    sed -r 's/.sam/.bam/g')
samtools view \
	 -b \
	 -S \
	 -t /home/lxm416/hg19/hg19.fa.fai \
	 -o $crickBamFile \
	 $crickFile
rm $crickFile
crickSortedFile=$(echo $crickBamFile | \
		       sed -r 's/.bam/.sorted.bam/g')

samtools sort -o $crickSortedFile $crickBamFile
rm $crickBamFile

# pileup
crickPileFile=$(echo $crickSortedFile | \
                                sed -r 's/.sorted.bam/.pile.txt/g')
samtools mpileup \
	 -f /home/lxm416/hg19/hg19.fa \
	 -o $crickPileFile \
	 $crickSortedFile
rm $crickSortedFile
crickCPGfile=$(echo $crickPileFile | \
		    sed -r 's/.pile.txt/.cpg.txt/g')

methyl-seq-level.strandOnly \
    /home/lxm416/hg19/hg19.CpG.sites.list \
    $crickCPGfile \
    $crickPileFile \
    $sampleId \
    5 \
    reverse
rm $crickPileFile
