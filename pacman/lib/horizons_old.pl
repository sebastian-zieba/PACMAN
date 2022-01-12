#!/usr/bin/perl
# 
# Access JOP's Horizons web interface.
#
# This script is an adapted version of a sample
# perl script written 2006-Mar-02 by Alain B. Chamberlin (JPL/Caltech).
# Unfortunately, I am unable to locate the script's url anymore.
#
# Use of this script is at your own risk.
#

use strict;
use LWP::UserAgent;

my @data = ( #_{
	"COMMAND= -48",
	"CENTER= 500\@0",
	"MAKE_EPHEM= YES",
	"TABLE_TYPE= VECTORS",
	"START_TIME= $ARGV[0]",
	"STOP_TIME= $ARGV[1]",
	"STEP_SIZE= 5m",
	"OUT_UNITS= KM-S",
	"REF_PLANE= FRAME",
	"REF_SYSTEM= J2000",
	"VECT_CORR= NONE",
	"VEC_LABELS= YES",
	"VEC_DELTA_T= NO",
	"CSV_FORMAT= NO",
	"OBJ_DATA= YES",
	"VEC_TABLE= 3",
); #_}

for (@data) { #_{

# Remove any spaces surrouning '=' for compactness
  s/ *= */=/;

# Escape special URL characters
  s/ /%20/g;
  s/\&/%26/g;
  s/;/%3B/g;
  s/\?/%3F/g;

} #_}

# Create URL
my $url = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&';
$url .= join('&', @data);

# print "URL: $url\n";
print "URL: $url\n";
my $ua  = new LWP::UserAgent;
my $req = new HTTP::Request;

$req->method("GET");
$req->url($url);
my $res = $ua->request($req);
die "$url\n" . $res->status_line . "\n" if ( $res->is_error );

print $res->content;

