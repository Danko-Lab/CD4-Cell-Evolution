#!/usr/bin/perl
#	Outputs the upstream promoter sequence froma bed file of transcripts.
#	usage: perl this.pl INPUTFILE UPSTREAM DOWNSTREAM MIN_SIZE
open(FILE, "< $ARGV[0]");
$up = $ARGV[1];
$dn = $ARGV[2];
$sz = $ARGV[3];

sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

while(<FILE>) {
	@SPL = split("\t");
	$chrom = $SPL[0];
	$start = $SPL[1];
	$end   = $SPL[2];
	$name  = $SPL[3];
	$score = $SPL[4];
	$str   = trim($SPL[5]);
	if($str eq "+") {
		if(($end-$start)>$sz) {
			print "$chrom\t".($start-$up)."\t".($start+$dn)."\t$name\t$score\t$str\t$start\t$end\n";
		}
	}
	else {
	if($str eq "-") {
                if(($end-$start)>$sz) {
 	               print "$chrom\t".($end-$dn)."\t".($end+$up)."\t$name\t$score\t$str\t$start\t$end\n";
		}
	}
	else {
		print "$str";
		exit;
	}
	}
}
