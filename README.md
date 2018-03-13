# motifSearch.pl
## version 2.1

  What is it?
  -----------

  Searches for perl regular expressions in sequence files.  This script outputs a tab-delimited file with coordinates.

## INSTALLATION

To install this module type the following:

    perl Makefile.PL
    make
    make install

And optionally (to remove unnecessary files):

    make clean

## RUNNING

To get the usage:

    motifSearch.pl

To get a detailed usage:

    motifSearch.pl --extended

To get help:

    motifSearch.pl --help

Example run:

    motifSearch.pl -s 'AUAU.{1,100}AUAU.{25,50}(?i:ggg)' -e -i rna.fa --verbose

## DEPENDENCIES

This module comes with a pre-release version of a perl module called "CommandLineInterface".  CommandLineInterface requires these other modules and libraries:

  Getopt::Long
  File::Glob

## COPYRIGHT AND LICENCE

See LICENSE
