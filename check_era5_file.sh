#!/bin/bash -f

filesize_max1=`ls -lS ERA5*| head -n 1| awk '{print $5}'`
filesize_max2=`ls -lS ERA5*| head -n 2|tail -n 1| awk '{print $5}'`

# ------------

[ $filesize_max1 -ne $filesize_max2 ] && exit

while true
do

    filesize_min=`ls -lS ERA5*| tail -n 1| awk '{print $5}'`

    [ $filesize_max1 -eq $filesize_min ] && exit

    filename_min=`ls -lS ERA5*| tail -n 1| awk '{print $9}'`

    echo $filesize_max1
    echo $filesize_min
    echo $filename_min

    rm $filename_min

    # exit

done