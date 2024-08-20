#!/bin/bash
if [ $# == 0 ] ; then
	echo "[usage:] sh qdel_job.sh *db"
fi

all_db=$*
for i in $all_db
do
	if [ "${i##*.}" = "db"  ] ; then
		sqlite3 $i 'select JOBID from JOB;' | xargs -i echo qdel {}
		sqlite3 $i 'select JOBID from JOB;' | xargs -i  qdel {}
	fi
done 

