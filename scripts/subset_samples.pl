#!/usr/bin/perl


open(F,"<$ARGV[1]");

@samples = <F>;
chomp @samples;

foreach $s (@samples){
	$lookup{$s} = 1;
}

close(F);


open(G,"<$ARGV[0]");
@lines = <G>;
chomp @lines;


foreach $line (@lines){

	@data = split(/\t/,$line);
	
	if(defined($lookup{$data[0]})){
		print "$line\n";
	}

}

close(G);
