#! /bin/sh

if [ "x$1" == x ]; then
 FILE=CAMB.tar.gz
else
  FILE=$1
fi

if [ ! -f $FILE ]; then
  echo "Please visit the CAMB website (opening now, hopefully)"
  echo "to register for and download CAMB into this directory."
  python -c 'import webbrowser;webbrowser.open("http://camb.info/CAMBsubmit.html")'
  exit 1
fi

tar -xzvf $FILE

trap "rm -rf tmp" EXIT
