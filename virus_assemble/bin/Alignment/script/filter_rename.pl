my ($filter,$info) = @ARGV;
my (%hash,@report);
open INFO,$info or die $!;
while(<INFO>){
	chomp;
	my ($bms_id,$report_id) = (split/\t/)[0,-2];
	$hash{$bms_id} = $report_id;
	push @report,$report_id;
}
close INFO;

my (%info,@subtitle);
open FL,$filter or die $!;
chomp(my $title = <FL>);
my @title = split/\t/,$title;
while(<FL>){
	chomp;
	my @tmp=split/\t/;
	my $subtitle = $tmp[0];
	push @subtitle,$subtitle;
	for (my $i=1;$i<=$#tmp;$i++){
		if (exists $hash{$title[$i]}){
			$info{$hash{$title[$i]}}{$subtitle} = $tmp[$i];
		}else{
			$info{$title[$i]}{$subtitle} = $tmp[$i];
		}
	}
}
close FL;
my $header = join("\t",@report);
print "Sample\t$header\n";
for(@subtitle){
	my @values;
	for (my $i=0;$i<=$#report;$i++){
		push @values,$info{$report[$i]}{$_}; 
	}
	my $value = join("\t",@values);
	print "$_\t$value\n";
}

