#!/bin/bash

if [ $# -ne 4 ]; then
	echo "Usage: `basename $0` <user> <pass> <email_to_share> <job_id|job_id_file>"
	exit
fi 
user=$1
pass=$2
email=$3
jobid=$4

cookie=`mktemp`
url=http://rast.nmpdr.org/rast.cgi\?page=JobShare\&email=$email\&action=share_job\&share_job=%20Share%20job%20with%20this%20user%20or%20group%20

function share {
	wget --load-cookie $cookie $url\&job=$1
}

wget --save-cookie $cookie http://rast.nmpdr.org/rast.cgi\?login=$user\&password=$pass\&action=perform_login
if [ -f $jobid ]; then
	for id in `cat $jobid`; do
		share $id
	done
else
	share $jobid
fi

rm $cookie
