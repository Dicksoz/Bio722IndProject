#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

#D I M S are the letters in the BWA CIGAR strings (Deletion, Insertion,
#Match/Mismatch, and Soft Clip)

sub MultiOverlapCalc($$@);


my %opts;
getopts('m:o:hSs',\%opts);
my $MismatchProp = exists($opts{m}) ? $opts{m} : 0;
if($MismatchProp < 0 || $MismatchProp > 1){
    die "Mismatch[-m] must be between 0 and 1\n";
}
my $OverlapProp = exists($opts{o}) ? $opts{o} : 0;
if($OverlapProp < 0 || $OverlapProp > 1){
    die "Overlap[-o] must be between 0 and 1\n";
}
my $bHeader = exists($opts{h}) ? 1 : 0;
my $bSoft = exists($opts{S}) ? 1 : 0;
my $bOutside = exists($opts{s}) ? 1: 0;
if($bOutside && !$bSoft){
    print STDERR "[WARNING] -s flag has no effect without -S flag\n";
    $bOutside = 0;
}

unless(@ARGV){
    die "Usage: perl [-options] $0 samFile\n".
        "samFile is a sam formated file\n".
        "\t-m float\tThe maximum allowable proportions of mismatches (Default: 0)\n".
        "\t\tRequires that the lines of the SAM file contain the the NM field\n".
        "\t-o float\tThe proportion of soft-clipped read which must overlap the region of interest\n".
        "\t\tRequires that the region have an field with the format ';interval=s1-e1,s2-e2,...;'\n".
        "\t\t(Default: 0 - Disabled)\n".
        "\t-h\tPrint header Lines\n".
        "\t-S\tExclude soft clipped reads\n".
        "\t-s\tUsed with -S: Keep soft clipped reads where the soft clipped portion is outside the region\n".
        "\t\tRequres samFile to have a header\n";
}


my $OrganismCountDict;
my %RegionLengthDict;

while(my $line = <>){
    if($line =~ /^@/){
        print $line if($bHeader);
        if($bOutside){
            if($line =~ /\@SQ/){
                if($line =~ /^\@SQ\tSN:(.+)\tLN:([0-9]+)$/){
                    $RegionLengthDict{$1} = $2;
                } else {
                    print STDERR "[WARNING] Could not get LN field for \@SQ at $.\n";
                }
            }
        }
        next;
    }
    my @SAMfieldList = split(/\t/,$line);
    my $bClipped = $SAMfieldList[5] =~ /S/;
    my $nMismatches = (split(/:/,$SAMfieldList[11]))[2];
    $SAMfieldList[5] =~ /([0-9]+)M/;
    my $scLength = $1;
    my $mapPos = $SAMfieldList[3];
    if($bSoft && $bClipped) {#CIGAR String indicates soft clipped
        if($bOutside){
            if($SAMfieldList[5] =~ /^[0-9]+S/ && $mapPos > 1){ #Left Soft clipping not outside region
                next;
            }
            if($SAMfieldList[5] =~ /[0-9]+S$/){ #Right Soft Clipping
                #Check if soft clipped read ends with the region
                next unless($mapPos + $scLength - 1 == $RegionLengthDict{$SAMfieldList[2]});
            }
        } else {
            next;
        }   
    }
    unless(defined($scLength)){
        print STDERR "[WARNING] Could not find NM field for $SAMfieldList[0]: Skipping\n";
        next;
    }
    if($nMismatches / $scLength > $MismatchProp){#Too many mismatches
        next;
    }
    my @regionFieldList = split(/;|=/,$SAMfieldList[2]);
    my $Organism = shift(@regionFieldList);
    my %regionFieldDict = @regionFieldList;
    if($OverlapProp){ #If an overlap requirement is requested
        unless(exists($regionFieldDict{interval})){
         print STDERR "[WARNING] Could not find interval field in region name for $SAMfieldList[0]: Skipping\n";
         next;
        }
        my @intervalList = split(/,/,$regionFieldDict{interval});
        my $scOverlap = MultiOverlapCalc($mapPos,$scLength,@intervalList);
        if($scOverlap / $scLength < $OverlapProp){ #Not enough overlap
            next;
        }
    }
    print $line;
}


sub MultiOverlapCalc($$@){
    my ($pos,$len,@intervalList) = @_;
    my @olList;
    my $end = $pos+$len-1;
    my $overlap = 0;
     #print "@intervalList\n";
    #Restrict to the intervals which the read overlaps
    while(my $int = shift(@intervalList)){
        my($s,$e) = split('-',$int);
        unless($end < $s or $pos > $e){
            push(@olList,$int)
        } 
    }
     #print "@olList\n";
    if(@olList){
        #Sort selceted intervals by start position
        @olList = sort {(split('-',$a))[0] <=> (split('-',$b))[0]} (@olList);
         #print "@olList\n";
        #Collapse all intervals into intervals intervals
        my($s0,$e0);
        while(my $int = shift(@olList)){
            my($s,$e) = split('-',$int);
            if(defined($s0)){
                if($s < $e0){
                    $e0 = $e if ($e >= $e0);;
                } else {
                    push(@intervalList,join('-',($s0,$e0)));
                    $s0 = $s;
                    $e0 = $e;
                }
            } else {
                $s0 = $s;
                $e0 = $e;
            }
        }
        push(@intervalList,join('-',($s0,$e0)));
    }
     #print "@intervalList\n";
    #Calculate amount of overlap  across all independant intervals
    if(@intervalList){
        foreach my $int (@intervalList){
            my($s,$e) = split('-',$int);
            if($end >= $s and $pos <= $e){
                 #print $overlap,"\n";
                $overlap += $end - $s + 1 - (($pos > $s) ? ($pos - $s) : 0) - (($end > $e) ? ($end - $e) : 0);
            } else {
                print STDERR "Independant Interval Not overlapping read\n";
            }
        }
    }
    return $overlap;
 }
