#!/bin/bash
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
# area bounds:

west=77.9
east=79.5
south=40.3
north=41.4
region=$west"/"$east"/"$south"/"$north

#time interval:
t1=20240110
t2=20240205

dir0="$PWD"
trackname=DES136
track_number=`echo $trackname|cut -c 4-`
track_dir=`echo $trackname|cut -c -3`
cd $dir0
#echo "get_s1_data.csh R$west/$east/$south/$north 20180101 20190101 $track_number Y"
echo "get_s1_data.csh R$region $t1 $t2 $track_number ASF Y"
#exit 1
get_s1_data.csh R$region $t1 $t2 $track_number ASF Y
echo "get_s1_data.csh R$region $t1 $t2 $track_number ASF Y"
echo "make_s1_bounds.csh $region $track_dir"
 
#exit 1
nnew=`cat $dir0/$trackname/data.new |wc -l` 
 
if [ "$nnew" -gt 0 ]; then
   echo "Job started on "`date`
   cd $dir0/$trackname
   if [ ! -e frames.ll ]; then
    imgtmp=`tail -1 raw/downloaded.list`
    cp $dir0/$trackname/raw/$imgtmp".SAFE"/preview/map-overlay.kml .
    echo "make_s1_bounds.csh $region $track_dir" > make_frame
    make_s1_bounds.csh $region $track_dir
   fi
fi

if [ "$nnew" -gt 0 ]; then 
   echo "Job started on "`date` 
   cd $dir0/$trackname 
   if [ -e frames.ll ] && [ "$nnew" -gt 1 ]; then 
     echo "Creating new frames ..." 
     echo "Downloading new orbits ..." 
#     check_s1_res_orbit.csh 1 
     make_s1a_frame.csh data_new.list frames.ll 
 
     echo "Start the processing ..." 
exit 1 
     cd $dir0/$trackname/FRAME_1 
     ndiff=`comm -23 dates.merge.new dates.merge.old|wc -l` 
     if [ "$ndiff" -gt 0 ]; then 
       batch_s1.csh  batch.config F.list 
       echo "Now start the post-processing "`date` 
       cd $dir0/$trackname/FRAME_1 
# 
####   The following line is for doing the time series analysis 
#       $dir0/$trackname/FRAME_1/process_ts.csh   
# 
#####       ###################################################### 
       echo "post-processing finished on " `date` 
       echo "Job finished on " `date` 
       cd $dir0 
      echo "That's all folks ..." 
     else 
       echo "Now new good data from Framing ...." 
       exit 1 
     fi 
   else 
     echo "No frame definition file found .." 
     exit 
   fi 
 
else 
   echo "No new data. Good Bye!!!" 
fi 
 
