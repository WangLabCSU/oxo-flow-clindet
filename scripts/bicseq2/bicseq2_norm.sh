#!/bin/bash
cleanup() {
        pkill -P $$
        kill 0
}

for sig in INT QUIT HUP TERM; do
        trap "
            cleanup
            trap - $sig EXIT
            kill -s $sig "'"$$"' "$sig"
done
declare -a pids;

output=$1
reference=$2
map_file=$3
bam=$4
mod_samtools=$5
threads=$6

read_length=150
fragment_size=550


samtools stats -@ $threads $bam > $bam.bam.stats & pids[0]=$!;

read_length=`cat $bam.bam.stats | grep ^SN | cut -f 2- | grep "average length" | awk -F "\t" '{split($2,f,"."); print f[1]}'`
fragment_size=`cat $bam.bam.stats | grep ^SN | cut -f 2- | grep "insert size average" | awk -F "\t" '{split($2,f,"."); print f[1]}'`

mkdir -p $output
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > $output/bicseq2.configFile;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "$chr\t$reference/$chr.fa\t${map_file}/hg19.50mer.CRC.chr$chr.txt\tWES/b37/results/cnv/paired/bicseq2/CGGA_653/samtools/$chr.seq\tWES/b37/results/cnv/paired/bicseq2/CGGA_653/samtools/$chr.norm.bin" >> $output/bicseq2.configFile;
done

mod_samtools view -U BWA,WES/b37/results/cnv/paired/bicseq2/samtools/,N,N $bam 

mkdir -p $output/bicseq2.tmp;
NBICseq-norm.pl -l $read_length -s $fragment_size --tmp $output/bicseq2.tmp --fig $output/norm.pdf  $output/bicseq2.configFile $output/bicseq2.NB_parameters