#!/usr/bin/perl -s
#
# A script to help writing the skeleton of a SusyEventSelector
#
# Original author: F. Ronga (ETH) - May 23, 2008
#
# $Id: $
#

use strict;

$| = 0; # No buffering

if ( $main::h || $main::help ) {
  chomp(my $cmd = `basename $0`);
  print <<EOF
*
* $cmd
*
* A simple script to help write the skeleton of a SusyEventSelector.
*
* The script should be run from the directory of the package where 
* the selector is to be created. E.g., SusyAnalysis/AnalysisSkeleton.
* It will create two files:
* - interface/<Name>.h with the skeleton interface
* - plugins/<Name>.cc  with the skeleton code
* where <Name> is the name of the selector: it can either be given
* as an option on the command line or during the creation process.
* (Missing directories will be created if they do not exist yet.)
*
* The CMSSW environment has to be set before running this script
* (just run "cmsenv" to set it).
*
* Usage: $cmd [-name=<selector name>]
*
EOF
    ;
  exit(0);
}

### Check environment
my $cmssw_base = $ENV{"CMSSW_BASE"};
if ( !$cmssw_base ) {
  die "\n*** Couldn't get CMSSW environment: please go to your local area\n"
    ."*** and simply run \"cmsenv\"\n";
}

### Get name of selector
my $selectorName;
if ( $main::name ) {
  # Was given on input
  $selectorName = $main::name;
} else {
  # Prompt user
  print "Please give the name of the selector (without extension): ";
  chomp($selectorName  = <STDIN>);
}

# Check that we have a name...
if ( length($selectorName) < 2 ) {
  die "Selector name ($selectorName) too short: please use at least two letters...\n";
}

# Complete name
if ( $selectorName !~ /Selector/ ) {
  print "*** Selector \"$selectorName\" does not contain \"Selector\" in the name.\n";
  print "--> Shall I append it (".$selectorName."Selector) [y*/n] ";
  my $answ = <STDIN>;
  if ( $answ !~ /^n/i ) {
    $selectorName .= "Selector";
  }
}

# Print some information
chomp(my $cwd = `pwd`);
print <<INFO
 Will create files:
 - $cwd/interface/$selectorName.h
 - $cwd/plugins/$selectorName.cc
INFO
  ;


### Check existence of directories
foreach my $dir ( "interface", "plugins" ) {
  if ( opendir(DIR,$dir) ) { closedir(DIR); } # Exists
  else
    {
      print "Shall I create directory $cwd/$dir? [y*/n] ";
      my $answ = <STDIN>;
      if ( $answ !~ /^n/i ) {
        mkdir( $dir ) or die "Couldn't create $cwd/$dir: $!";
      } else {
        die "Please create directory $cwd/$dir first.\n";
      }
    }
}

### Create files
(my $root = $cwd) =~ s/.*CMSSW_.*?\/src\///; # Local directory relative to base
(my $defineRoot = $root) =~ s@/@\_@g;
my $userName = $ENV{'USER'};
chomp(my $date = `date`);

# Header
my $file = "interface/$selectorName.h";
open(HEADER,">$file") or die "Couldn't create $file: $!";
print HEADER <<EOF
#ifndef $defineRoot\_$selectorName\_h_
#define $defineRoot\_$selectorName\_h_
///
/// $selectorName
///
/// \\author $userName - $date
///
/// \$Id: \$
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
class $selectorName : public SusyEventSelector {
public:
  $selectorName (const edm::ParameterSet&);
  virtual bool select(const edm::Event&) const;
  virtual ~$selectorName () {}

 private:
  // Define parameters here
};
#endif
EOF
  ;
close(HEADER);

# Source code
$file = "plugins/$selectorName.cc";
open(SRC,">$file") or die "Couldn't create $file: $!";
print SRC <<EOF
#include "$root/interface/$selectorName.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//________________________________________________________________________________________
$selectorName\::$selectorName (const edm::ParameterSet& pset ) :
  SusyEventSelector(pset)
{

  // Define variables to store
  // defineVariable("myVariable");

}

//________________________________________________________________________________________
bool
$selectorName\::select (const edm::Event& event) const
{

  // Reset cached variables
  resetVariables();

  // Compute variables

  // Store variables
  // setVariable("myVariable",myValue);

  // Return selection
  // return (myValue > minValue);

}

//__________________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, $selectorName, "$selectorName");

EOF
  ;
close(SRC);

print "\n  All done!\n";
