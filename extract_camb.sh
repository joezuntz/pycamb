#! /bin/sh

if [ "x$1" == x ]; then
  echo trying do get CAMB.tar.gz from camb.info
  cat<< EOF
Licence
You are licensed to use this software free of charge until January 2013 on condition that:
    * Any publication using results of the code must be submitted to arXiv at the same time as, or before, submitting to a journal. arXiv must be updated with a version equivalent to that accepted by the journal on journal acceptance.
    * If you identify any bugs you report them as soon as confirmed 

  If you do not agree with the Licence, contact camb.info, and abort from here with Ctrl-C.
  Otherwise, press <RETURN> to continue.
EOF
  read
  wget http://camb.info/CAMB.tar.gz
  FILE=CAMB.tar.gz
else
  FILE=$1
fi

if [ ! -f $FILE ]; then
  echo $FILE is not a file.
  exit 1
fi

tar -xzvf $FILE

trap "rm -rf tmp" EXIT
