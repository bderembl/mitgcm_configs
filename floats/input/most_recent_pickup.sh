#!/bin/bash

# script from https://gist.github.com/rabernat/90cd307e9cdf8c774a42fd6ed2264fe9

pmeta=$( ls -t pickup.*.meta |head -1 )
plabel=$( echo $pmeta | sed 's/pickup.\(.*\).meta/\1/' )

# strip from data
sed -i.bak -n -e '/pickupSuff/I!p' data

if [ "${plabel:0:1}" == "c" ]
then
	# we have a partial checkpoint
	iter=$( grep timeStepNumber $pmeta | sed -e 's/.*\[ \(.*\) \];/\1/' )
	extra="\n pickupSuff='$plabel',"
else
	iter=$plabel
fi

# update iteration number
sed -i.bak -e "s/nIter0=\(.*\)/nIter0=$iter,$extra/I" data

echo "Updated data: pickup $iter ($plabel)"
