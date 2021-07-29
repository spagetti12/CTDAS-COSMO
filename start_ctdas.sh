# CarbonTracker Data Assimilation Shell (CTDAS) Copyright (C) 2017 Wouter Peters. 
# Users are recommended to contact the developers (wouter.peters@wur.nl) to receive
# updates of the code. See also: http://www.carbontracker.eu. 
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation, 
# version 3. This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this 
# program. If not, see <http://www.gnu.org/licenses/>. 

#!/bin/bash
set -e

cat > heredocfile.txt <<_EOF_
    where <projectdir> is a base folder for the project
    and <projectname> is a name to use for the project.

    !! A folder projectdir/projectname will be created !!
_EOF_

while getopts "h" opt; do
  case $opt in
    h) cat heredocfile.txt
      exit 1
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
      exit 1
    ;;
    *) cat heredocfile.txt
      exit 1
    ;;
  esac
done

EXPECTED_ARGS=2
E_BADARGS=666

if [ $# -ne $EXPECTED_ARGS ]
then
  echo ""
  echo "Usage: `basename $0` projectdir projectname"
  cat heredocfile.txt
  exit $E_BADARGS
fi

echo "New project to be started in folder $1"
echo "               ...........with name $2"

rootdir=$1/$2
rundir=$1/$2/exec
sedrundir=$1/$2/exec

if [ -d "$rootdir" ]; then
    echo "Directory already exists, please remove before running $0"
    exit 1
fi

mkdir -p ${rundir}
rsync -au --cvs-exclude * ${rundir}/
cd ${rundir}

echo "Creating jb file, py file, and rc-file"
sed -e "s/template/$2/g" template.jb > $2.jb
sed -e "s/template/$2/g" template.py > $2.py
sed -e "s,template,${rootdir},g" template.rc > $2.rc
rm -f template.py
rm -f template.jb
rm -f template.rc
rm -f start_ctdas.sh

chmod u+x $2.jb

echo ""
echo "************* NOW USE ****************"
ls -lrta $2.*
echo "**************************************"
echo ""
cd ${rundir}
pwd

