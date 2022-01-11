#!/bin/bash
[ "$#" -lt 2 ] && echo usage: $0 PeptideLevel.csv[.gz] TransitionLevel.csv && exit 1

case $1 in
  *.gz) gunzip -c $1;;
  *) cat $1;;
esac | awk -F"," '
NR==1 {
  for (f=1; f<=NF; f++) {
    cols[$f] = f
  }
  print "ProteinName,FullPeptideName,Charge,RTsubgroup,iRT,culture,Intensity,aggr_Fragment_Annotation,aggr_Peak_Area"
  next
} {
  split($cols["aggr_Fragment_Annotation"],farr,";")
  split($cols["aggr_Peak_Area"],parr,";")
  for (f=1; f<=length(farr); f++) {
     print $cols["ProteinName"] "," $cols["FullPeptideName"] "," $cols["Charge"] ","  $cols["RTsubgroup"] ","  $cols["iRT"] "," $cols["culture"] ","  $cols["Intensity"] "," farr[f] "," parr[f]
  }
}' > $2
