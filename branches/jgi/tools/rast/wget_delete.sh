#!/bin/bash

if [ $# -ne 3 ]; then
	echo "Usage: `basename $0` <user> <pass> <job_id|job_id_file>"
	exit
fi 
user=$1
pass=$2
jobid=$3

cookie=`mktemp`
url=http://rast.nmpdr.org/rast.cgi\?page=JobDelete\&confirm=%20Delete%20this%20job%20

function share {
	wget --load-cookie $cookie $url\&job=$1
}

wget --save-cookie $cookie http://rast.nmpdr.org/rast.cgi\?login=$user\&password=$pass\&action=perform_login
if [ -f $jobid ]; then
	for id in `cut -f 2 $jobid`; do
		share $id
	done
else
	share $jobid
fi

rm $cookie
