#!/usr/bin/perl

$bam_list = $ARGV[0];
$tsv_path = $ARGV[1];

#print "$bam_list\n";

#open(F,"<SPARK.bam.list");
open(F,"<$bam_list");
@bams = <F>;

chomp @bams;

foreach $b (@bams){

	$b =~ /.*\/(.*?)\.bam/;
	$sample = $1;

	$lookup{$sample} = $b;
	#print "$sample\n";
	
	
}

close(F);

#@tsvs = glob("/mnt/home/mlek/mlek/readviz/SPARK_WGS/chr1_tsvs/*.bgz");
@tsvs = glob("$tsv_path/*.bgz");

foreach $t (@tsvs){
	$t =~ /.*\/(.*?)\.tsv.bgz/;
	$sample = $1;

	if(defined($lookup{$sample})){
		if(-e "$lookup{$sample}.bai"){
			print "$sample\t$t\t$lookup{$sample}\t$lookup{$sample}.bai\n";
		}
		else{
			print "$sample\t$lookup{$sample}\t$bai2\tINDEX_NOT_FOUND\n";
		}

	}
	else{
		#print "$sample NOT_FOUND\n";
	}
}
