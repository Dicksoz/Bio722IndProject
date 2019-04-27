#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
#use Parallel::ForkManager;

unless(@ARGV){
    die "Usage: ".basename($0)."HeaderMap.tab exactFilteredSAM\n";
}

#Load up map
my %OrgDict;
my %OrgCount;
my $MapFile = shift(@ARGV);
open(my $fh, $MapFile) or die "Could not open $MapFile: $!\n";
while(my $line = <$fh>){
    chomp($line);
    my @Taxa = split(/\t/,$line);
    my $head = shift(@Taxa);
    my $org = join('|',@Taxa);
    $org =~ s/ //g;
    $OrgDict{$head} = $org;
    $OrgCount{$org} = 0;
}
close($fh);


#Count Reads for each organism
my %ReadRecord;
while(my $line = <>){
    next if($line =~ /^@/);
    chomp($line);
    #  print $line,"\n";
    my @fieldList = split(/\t/,$line);
    #print "@fieldList\n";
    if(exists($ReadRecord{$fieldList[0]})){ #Only count one if both read pairs mapped
        next;
    } else {
        $ReadRecord{$fieldList[0]} = 1;
        my $org = $OrgDict{$fieldList[2]};
        $OrgCount{$org}++;
    }
}

#Output Read Counts by Organism
foreach my $org (sort(keys(%OrgCount))){
    print "$org\t$OrgCount{$org}\n";
}
