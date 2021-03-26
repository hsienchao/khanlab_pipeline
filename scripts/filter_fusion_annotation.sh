annotation=$1
perl -lane '{s/,/\t/g;@s=split(/\t/);print if ($s[11] eq "SpanningReads" || $s[11] > 2)}' $annotation
