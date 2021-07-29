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
    where <rootdir> is a base folder for the project
    where <projectsource> is a source folder for the project, to clone
    and <projectclone> is a name to use for the cloned project.

    !! A folder rootdir/projectclone will be created !!
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

EXPECTED_ARGS=3
E_BADARGS=666

if [ $# -ne $EXPECTED_ARGS ]
then
  echo ""
  echo "Usage: `basename $0` rootdir projectsource projectclone"
  cat heredocfile.txt
  exit $E_BADARGS
fi

echo "New project to be started in folder $1"
echo "               ...........with name $3"
echo "             ...........cloned from $1/$2"

sourcedir=$1/$2/exec
rundir=$1/$3/exec
sedrundir=$1/$3/exec

if [ -d "$rundir" ]; then
    echo "Directory already exists, please remove before running $0"
    exit 1
fi

mkdir -p ${rundir}
rsync -ru --cvs-exclude --exclude=*nc ${sourcedir}/* ${rundir}/
rsync -ru --cvs-exclude ${sourcedir}/da/analysis/*nc ${rundir}/da/analysis/
cd ${rundir}

echo "Modifying jb file, py file, and rc-file"
sed -e "s/$2/$3/g" ${sourcedir}/$2.jb > $3.jb
sed -e "s/$2/$3/g" ${sourcedir}/$2.py > $3.py
sed -e "s/$2/$3/g" ${sourcedir}/$2.rc > $3.rc
rm -f clone_ctdas.sh
rm -f $2.jb $2.rc $2.py
make clean

chmod u+x $3.jb

echo ""
echo "************* NOW USE ****************"
ls -lrta $3.*
echo "**************************************"
echo ""
cd ${rundir}
pwd

