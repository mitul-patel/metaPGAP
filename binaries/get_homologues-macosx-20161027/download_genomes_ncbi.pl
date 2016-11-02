#!/usr/bin/perl -w
# 2013-6 Bruno Contreras-Moreira, Pablo Vinuesa
# download_genomes_ncbi.pl

# Attempts to download a list of genomes from the NCBI.
# Queries 1) e-utils, 2) the FTP site and 3) the WGS repository.

# Uses wget and zcat binaries, which should be installed on most systems.

use strict;
use File::Fetch;

my $VERSION  = 2.0;

# http://www.ncbi.nlm.nih.gov/books/NBK25500/
my $EFETCHEXE = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text&id=';
my $NCBIHOST  = 'ftp://ftp.ncbi.nlm.nih.gov';
my $WGSURL    = 'http://www.ncbi.nlm.nih.gov/Traces/wgs/?download='; 

my $COMPLETE_GENOMES_DIR = '/genomes/all/'; # includes also refseq folders Jan2016
my $LOCAL_GENOMES_DIR    = './';            # where to store downloaded files

my $EXTENSION    = 'gbk';
my $WGSEXTENSION = 'gbff';

my $SLEEPTIME    = 1; # seconds between queries to avoid overloading remote servers

$File::Fetch::WARN = 0; # no warnings by default
$File::Fetch::DEBUG = 0;

#######################################################

my ($size,$localsize,$file,$localfile,$listfile,%wanted_genome,%new_name,$newname,$corr_filename);
my ($genome,$dir,$n_of_files,$single_file,$n_of_downloaded_files,$n_of_wanted_files);

if(!$ARGV[0])
{
  die "# usage: $0 v.$VERSION: <genome_list.txt>\n\n".
    "File <genome_list.txt> must contain names of genomes to be downloaded, one per row, as explained in the manual.\n".
    "Only the first column is required. Example:\n\n".
    "NC_010159                   Yersinia_pestis_Angola      # known accession\n".
    "GCA_000016445.1_ASM1644v1   Yersinia_pestis_Pestoides_F # complete genome\n".
    "AAYU01                      Yersinia_pestis_B42003004   # WGS genome, possibly a draft\n\n"
}
else{ $listfile = $ARGV[0] }

# 1) parse list file ###################################
$n_of_wanted_files = 0;
open(LIST,$listfile) || die "# cannot read $listfile\n";
while(<LIST>)
{
  next if(/^#/ || /^\s+/);
  ($genome,$newname) = split;
  $wanted_genome{$genome}{'order'} = $n_of_wanted_files;
  $wanted_genome{$genome}{'wanted'} = 1;
  $wanted_genome{$genome}{'downloaded'} = 0;
  $wanted_genome{$genome}{'appears'}++;
  if($wanted_genome{$genome}{'appears'} > 1)
  {
    die "# $0: error, genome $genome seems to appear twice in input list\n";
  }
  elsif($newname){ $new_name{$genome} = $newname }
  $n_of_wanted_files++;
}
close(LIST);
print "# number of genomes in '$listfile' = $n_of_wanted_files\n";

# 2) look for parsed genomes ###########################
$n_of_downloaded_files = 0;
chdir($LOCAL_GENOMES_DIR);


# 2.1) check FTP server of complete genomes #
print "\n# checking complete genomes...\n";
foreach $dir (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
{  
  if($new_name{$dir})
  {
    $single_file = "_".$new_name{$dir} . ".$EXTENSION";
  }
  else{ $single_file = "_".$dir . ".$EXTENSION"; }

  if(-s $single_file)
  {
    print "> $single_file (done)\n";
    $wanted_genome{$dir}{'downloaded'}=1;
    next;
  }
  
  eval
  {
    #GCA_000016445.1_ASM1644v1_genomic.gbff.gz
    $file = "$dir\_genomic.$WGSEXTENSION.gz"; 
    
    if(!-s $file)
    {
      #system("wget -q $NCBIHOST$COMPLETE_GENOMES_DIR$dir/$file -O $file");
      my ($ff,$content); 
      $ff = File::Fetch->new(uri =>"$NCBIHOST$COMPLETE_GENOMES_DIR$dir/$file");   
      if($ff->fetch( to => \$content ))
      {
        open(GBKFILE,">$file") || die "# cannot append to $file\n";
        binmode(GBKFILE);
        print GBKFILE $content;
        close(GBKFILE);
      }
      #else { warn "# cannot fetch $NCBIHOST$COMPLETE_GENOMES_DIR$dir/$file\n" }
    }  
    
    if(!-s $file)
    {
      unlink($file) if(-e $file); 
    }
    else
    {
      system("zcat $file > $single_file");   
      if(-s $single_file)
      {
        print ">$single_file\n";
        if(gbk_contains_CDSs($single_file))
        {
          $wanted_genome{$dir}{'downloaded'} = 1;
          $n_of_downloaded_files++;
        }
        else
        {
          $wanted_genome{$dir}{'downloaded'} = 1;
          print "# file $single_file does not contain annotated CDSs\n";
        }
      }
    }      
  };
  
  sleep($SLEEPTIME);
}

if($n_of_downloaded_files == $n_of_wanted_files)
{ 
  print "\n# number of handled genomes = $n_of_downloaded_files (out of $n_of_wanted_files)\n";
  die "\n# bye!\n"; 
}

# 2.2) check WGS repository #
print "\n# checking bacterial WGS genomes...\n";
foreach $dir (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
{
  next if($wanted_genome{$dir}{'downloaded'});

  if($new_name{$dir})
  {
    $single_file = "_".$new_name{$dir} . ".$EXTENSION";
  }
  else{ $single_file = "_".$dir . ".$EXTENSION"; }

  if(-s $single_file)
  {
    print "> $single_file (done)\n";
    $wanted_genome{$dir}{'downloaded'}=1;
    next;
  }
  
  eval
  {
    $file = "$dir.1.$WGSEXTENSION.gz";
  
    # wget hack
    system("wget -q $WGSURL.$file -O $file") if(!-s $file);
    
    ##File::Fetch fails as it creates tmp file with invalid name on my system
    #my ($ff,$content); 
    #$ff = File::Fetch->new(uri=>"$WGSURL$file");   
    #if($ff->fetch( to => \$content )){
    #  open(GBKFILE,">$file") || die "# cannot append to $file\n";
    #  binmode(GBKFILE);
    #  print GBKFILE $content;
    #  close(GBKFILE);
    #}else { warn "# cannot fetch $WGSURL$file\n" }
    
    if(!-s $file)
    {
      unlink($file);
    }
    else
    {
      my $outputOK = 1;
      open(OUTFILE,$file);
      while(<OUTFILE>)
      {
        if(/Cannot find/){ $outputOK = 0; last }
      }
      close(OUTFILE);
    
      if($outputOK)
      {
        system("zcat $file > $single_file");  
        if(-s $single_file)
        {
          print ">$single_file\n";
          if(gbk_contains_CDSs($single_file))
          {
            $wanted_genome{$genome}{'downloaded'} = 1;
            $n_of_downloaded_files++;
          }
          else
          {
            $wanted_genome{$genome}{'downloaded'} = 1;
            print "# file $single_file does not contain annotated CDSs\n";
          }
        }
      }
      else
      {
        unlink($file);
      }        
    }      
  };

  sleep($SLEEPTIME);
}

if($n_of_downloaded_files == $n_of_wanted_files)
{ 
  print "\n# number of handled genomes = $n_of_downloaded_files (out of $n_of_wanted_files)\n";
  die "\n# bye!\n"; 
}


# 2.3) finally query e-utils  in case name is an acession #
print "\n# checking GenBank accessions with NCBI e-utils efetch ...\n";
foreach $genome (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
{
  next if($wanted_genome{$genome}{'downloaded'});

  $file = $new_name{$genome} || $genome;
  $file = "_$file.$EXTENSION"; #print "$file $wanted_genome{$genome}{'order'}\n"; #next;
  if(-s $file)
  {
    print ">$file (done)\n";
    $wanted_genome{$genome}{'downloaded'} = 1;
    $n_of_downloaded_files++;
    next;
  }

  eval
  {
    foreach my $g (split(/,/,$genome)) # supports comma-sep accessions
    {
      my ($ff,$gbkcontent);
      $ff = File::Fetch->new(uri =>$EFETCHEXE.$g);
      if($ff->fetch( to => \$gbkcontent ))
      {
        open(GBKFILE,">>$file") || die "# cannot append to $file\n";
        print GBKFILE $gbkcontent;
        close(GBKFILE);
      }
    }

    if(-s $file)
    {
      print ">$file\n";
      if(gbk_contains_CDSs($file))
      {
        $wanted_genome{$genome}{'downloaded'} = 1;
        $n_of_downloaded_files++;
      }
      else
      {
        $wanted_genome{$genome}{'downloaded'} = 1;
        $corr_filename = $file;
        $corr_filename =~ s/^_//;
        rename($file,$corr_filename);
        print "# file $corr_filename does not contain annotated CDSs\n";
      }
    }
  }; 

  sleep($SLEEPTIME);
} 

print "\n# number of handled genomes = $n_of_downloaded_files (out of $n_of_wanted_files)\n";





sub gbk_contains_CDSs
{
  my ($gbkfile) = @_;
  my $n_of_CDSs = 0;
  open(GBK,$gbkfile) || die "# check_gbk_content: cannot read $gbkfile\n";
  while(<GBK>){ if(/^     CDS /){ $n_of_CDSs++ } }
  close(GBK);
  return $n_of_CDSs;
}

