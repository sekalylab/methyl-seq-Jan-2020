# methyl-seq-Jan-2020
# methyl-seq data processing

#

bF="$1"
pref=`echo "$bF" | sed 's/.bam//' | sed 's/ /_/g'`

samtools view "$bF" | awk 'BEGIN {FS = "\t";} $2 == 99 || $2 == 147 || $2 == 83 || $2 == 163'  > temp-556.sam
wc -l temp-556.sam
### frag > 0, then read upsteam its mate; frag < 0, then read downsteam its mate
awk 'BEGIN {FS = "\t";} ($9 > 0 && $4 <= $8) || ($9 < 0 && $4 >= $8)' temp-556.sam  > temp.fowReadUpstream.sam 
rm temp-556.sam
wc -l temp.fowReadUpstream.sam
### overlapped portion of read pair 
sam-overlap-rp  temp.fowReadUpstream.sam out-556 > temp-overlap-556
wc -l temp-overlap-556
rm temp.fowReadUpstream.sam
rm temp-overlap-556

awk 'BEGIN {FS = "\t";} $2 == 99 || $2 == 147' out-556 > out-556-watson
awk 'BEGIN {FS = "\t";} $2 == 83 || $2 == 163' out-556 > out-556-crick
rm out-556
wc -l out-556-watson out-556-crick

samtools view -b -S -t /home/lxm416/hg19/hg19.fa.fai -o out-watson.bam out-556-watson
rm out-556-watson
samtools sort out-watson.bam out-watson.sorted 
rm out-watson.bam 
samtools mpileup -Q 0 -q 20 -d 1000000 -f /home/lxm416/hg19/hg19.fa  out-watson.sorted.bam   > out.watson.pile
rm out-watson.sorted.bam
methyl-seq-level.strandOnly  /home/lxm416/hg19/hg19.CpG.sites.list $pref.watson.CpG.txt out.watson.pile $pref 10 forward

samtools view -b -S -t /home/lxm416/hg19/hg19.fa.fai -o out-crick.bam out-556-crick
rm out-556-crick
samtools sort out-crick.bam out-crick.sorted
rm out-crick.bam
samtools mpileup -Q 0 -q 20 -d 1000000 -f /home/lxm416/hg19/hg19.fa  out-crick.sorted.bam   > out.crick.pile
rm out-crick.sorted.bam
methyl-seq-level.strandOnly  /home/lxm416/hg19/hg19.CpG.sites.list $pref.crick.CpG.txt out.crick.pile $pref 10 forward
rm out.watson.pile out.crick.pile

#
#
