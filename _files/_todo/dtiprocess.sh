#! /bin/sh
#  Dewarping of B0 inhomogeneity distortion, correction of eddy current distortion and motion, and tensor calculation using FSL
#  Takuya Hayashi, Center for Life Science Technologies, RIKEN (takuya.hayashi@riken.jp, takuya.hayashi@gmail.com), Kobe, JP

VERSION="1.8" # released on Thu Aug 1 09:11:15 JST 2013

#---------variables and functions---------#
CMD=`echo $0 | sed -e 's/^\(.*\)\/\([^\/]*\)/\2/'`
inparg="$@";
CWD=`pwd`
dim=0;		           # reslice dti data in a specified dimension
ud=y-;		           # default undistortion direction
SL=10;		           # default signal loss threshold
outdir=dti;		    # default output dirname
tmpdir=dti_preprocess;   # default dir to save intermediate files
reportdir=report;	    # report dir
LF=dti_preprocess.log;   # default log filename
reg=1; 		    # registration between dti and fieldmap
refvol=0;	           # reference volume for motion correcion and eddy current distortion
test=0;                  # test or full calc
FORCEDIR=0;              # Use the actual directory name given (i.e. do not add + to make a new directory)
FLIRT=0;                 # use flirt (1) or mcflirt (0)
MCFLIRT=3;               # number of step in mcflirt
bvec_ecc=0;              # correction of bvector file using rotation parameter
FC=1                     # fieldmap correcion on (1) or off (0)

usage_exit() {
      cat <<EOF

  DTI calculation with correction for motion, eddy current distortion, B0 inhomogeneity distortion

  Version `echo "$VERSION" | awk '{print $1}'`

  Usage:
  
    $CMD -k <img> -b <bvals.txt> -r <bvecs.txt> [option]
      : corrects motion & eddy current distortion and calculates dti (w/o correction of B0 inhomogentity distortion)     
  
    $CMD -k <img> -t <num> -e <num> -f <img> -m <img> -b <bvals.txt> -r <bvecs.txt> [option]
      : corrects motion & eddy current distortion and B0 inhomogeneity and calculates dti

    $CMD -k <img> -t <num> -e <num> -f <img> -m <img> [option]
      : corrects B0 inhomogeneity disortion only for reference volume
  
    -k <img>    : DTI 4D data
    -t <num>    : DTI dwell time (ms)
    -e <num>    : DTI TE (ms)
    -f <img>    : B0 fieldmap image (radian/sec)
    -m <img>    : B0 fieldmap magnitude image
    -b <bvals.txt> : a text file containing a list of b-values
    -r <bvecs.txt> : a text file containing a list of b-vectors
  
    Option: 
    -u <x, x-, y, y-, z, or z->  : unwarp direction (default: y-)
    -s <num> : %signal loss threshold for B0 unwarping (default: 10)
    -v <num> : reslice resolution in mm (default: no reslicing)
    -K <img> : mask file for dti image (by default, this is calculated by bet)
    -M <img> : mask file for fieldmap magnitude image (by default, this is calculated by bet)
    -d <dir> : directory to save outputs (default: ./dti)
    -F       : Use the actual directory name given (i.e. do not add + to make a new directory)
    -R <num> : reference volume (default : 0th volume)
    -4   : 4-step registration in mcflirt (slow but more accurate, default is 3-step)
    -l   : use flirt instead of mcflirt (much slower but more accurate)
    -n   : do not register between fieldmap and dti data
    -B   : correct b-vector file using rotation parameters (output subdirectory, dti_bvecmc)
    -U   : do not mask DTI images

EOF
    exit 1;
}

Exit () {
    cd $CWD
    exit 1;
}

make_absolute (){
    file=$1
    if [ "`echo ${file} | head -c 1`" = "/" ]; then abs=${file}; else abs=${PWD}/${file}; fi
    echo $abs 
}

test_varimg (){
    var=$1
    if [ "x$var" = "x" ]; then test=0; else  test=`imtest $1`; fi
    echo $test
}

test_varfile (){
    var=$1
    if [ "x$var" = "x" ]; then test=0 ; elif [ ! -f $var ]; then test=0; else test=1; fi
    echo $test
}

transposematrix (){
    nl=`cat $1 | awk 'BEGIN {N=0}{N++} END{print N}'`
    nw=`cat $1 | awk 'NR==1 {print NF}'`
    if [ -e $2 ] ; then rm -f $2; fi
    touch $2
    i=1
    while [ $i -le $nw ] ; do
   	j=1
   	while [ $j -le $nl ] ; do
      	v=`cat $1 | awk 'NR=='$j' {printf "%f ",  $'$i'}'`
      	echo -n "$v " >> $2
      	j=`expr $j + 1`
   	done
    echo "" >> $2
    i=`expr $i + 1`
    done
}

threecolumnmeansd () {
    mean=`cat $1 | awk 'BEGIN {x=0;y=0;z=0;N=0};{x=x+$1;y=y+$2;z=z+$3;N=N+1}END {printf("%f, %f, %f",x/N,y/N,z/N)}'`
    xm=`echo $mean | awk '{print $1}'`; ym=`echo $mean | awk '{print $2}'`; zm=`echo $mean | awk '{print $3}'`
    sd=`cat $1 | awk 'BEGIN {x=0;y=0;z=0;N=0};{x=x+($1-"$xm")^2;y=y+($2-"$ym")^2;z=z+($3-"$zm")^2;N=N+1}END {printf("%f, %f, %f",sqrt(x/N),sqrt(y/N),sqrt(z/N))}'`
    echo $mean $sd
}


round () {
   LANG=C printf "%1.7f" $1
}
	
T () {
    E=0; if [ "$1" = "-e" ] ; then E=1; shift ; fi; cmd="$*"; echo $* | tee -a $LF; if [ "$E" != "1" ] ; then $cmd 2>&1 | tee -a $LF; fi; echo  | tee -a $LF
}

#------------- Preparatory process --------------------#
[ "$6" = "" ] && usage_exit

while getopts k:t:e:f:m:b:r:u:s:v:K:M:d:FR:l4nBU OPT
 do
 case "$OPT" in 
   "k" ) dti="$OPTARG";;
   "t" ) esp="$OPTARG";;
   "e" ) te="$OPTARG";;
   "f" ) dph="$OPTARG";;
   "m" ) mag="$OPTARG";;
   "b" ) bval="$OPTARG";;
   "r" ) bvec="$OPTARG";;
   "u" ) ud="$OPTARG";;
   "s" ) SL="$OPTARG";;
   "v" ) dim="$OPTARG";;
   "K" ) maskd="$OPTARG";;
   "M" ) maskf="$OPTARG";;   
   "d" ) outdir="$OPTARG";;
   "F" ) FORCEDIR=1;;
   "R" ) refvol="$OPTARG";;
   "l" ) FLIRT=1;;
   "4" ) MCFLIRT=4;;
   "n" ) reg=0;;
   "B" ) bvec_ecc=1;;
   "U" ) unmaskd=1;;
    * )  usage_exit;;
 esac
done;

if [ `test_varimg $dti` -eq 0 ]; then
   echo "ERROR: cannot find image for dti 4D data"
   exit 1;
else
   dti=`make_absolute $dti`
   dtidim4=`fslval $dti dim4`
fi

if [ "$dph" = "" ] && [ "$mag" = "" ] ; then FC=0
else
 if [ `test_varimg $dph` -eq 0 ]; then echo "ERROR: cannot find image for B0 fieldmap"; exit 1; else dph=`make_absolute $dph`; fi
 if [ `test_varimg $mag` -eq 0 ]; then echo "ERROR: cannot find image for B0 fieldmap magnitude"; exit 1; else mag=`make_absolute $mag`; fi
fi

if [ "x$maskf" != "x" ]; then
   if [ `test_varimg $maskf` -eq 0 ]; then echo "ERROR: cannot find image: $maskf"; exit 1; else maskf=`make_absolute $maskf`; fi; 
fi

if [ "x$maskd" != "x" ]; then
   if [ `test_varimg $maskd` -eq 0 ]; then echo "ERROR: cannot find image: $maskd"; exit 1; else maskd=`make_absolute $maskd`; fi
fi

if [ "$bvec" = "" ] && [ "$bval" = "" ] ;  then
    test=1
else
    if [ `test_varfile $bvec` -eq 0 ]; then echo "ERROR: no bvecs file specified"; exit 1;fi
    bvec=`make_absolute $bvec`; bvecl=`cat $bvec | awk 'END{print NR}'`; bvecw=`cat $bvec | wc -w`	
    if [ $bvecl != 3 ]; then echo "ERROR: bvecs file contains $bvecl lines, it should be 3 lines, each for x, y, z"; exit 1;fi
    if [ "$bvecw" != "`expr 3 \* $dtidim4`" ]; then echo "ERROR: bvecs file contains $bvecw words, it should be 3 x $dtidim4 words"; exit 1;fi
    if [ `test_varfile $bval` -eq 0 ]; then echo "ERROR: no bvals file specified"; exit 1;fi
    bval=`make_absolute $bval`; bvall=`cat $bval | awk 'END{print NR}'`; bvalw=`cat $bval | wc -w`
    if [ $bvall != 1 ]; then echo "ERROR: bvals file contains $bvall lines, it should be 1 lines"; exit 1;fi
    if [ $bvalw != $dtidim4 ]; then echo "ERROR: bvalc file contains $bvalw words, it should be $dtidim4 words"; exit 1;fi 
fi

if [ "$SL" -gt 100 ] || [ "$SL" -lt 0 ]; then echo "Error: signal loss threshold should be from 0 to 100";exit 1;fi

if [ "$FC" = "1" ] ; then
 if [ "x${esp}" = "x" ]; then echo "Cannot find option -t (echo space time)...exit."; exit 1;fi
 if [ "x${te}" = "x" ]; then echo "Cannot find option -e (echo time)...exit."; exit 1;fi
fi

outdir=`make_absolute $outdir`;

if [ $FORCEDIR = 0 ] ; then while [ -e $outdir ]; do outdir=${outdir}+; done; mkdir -p $outdir; elif [ $FORCEDIR = 1 ] ; then mkdir -p $outdir; elif [ ! -w $outdir ] ; then echo "ERROR: cannot create $outdir"; exit 1; fi;


#-----------Start calculation------------#
trap Exit 2
cd $outdir
mkdir -p $tmpdir
if [ -e $reportdir ]; then /bin/rm -Rf $reportdir;fi
mkdir -p $reportdir
mkdir -p .files; /bin/cp ${FSLDIR}/doc/fsl.css .files;
if [ -e report.html ]; then rm report.html;fi
touch report.html
echo "<HTML><HEAD><link REL=\"stylesheet\" TYPE=\"text/css\" href=\".files/fsl.css\"><TITLE>DTI_PREPROCESS REPORT</TITLE></HEAD><BODY BGCOLOR="#aaaaaa">" >> report.html
echo "<hr><B>DTI_PREPROCESS REPORT</B><BR>" >> report.html
echo "<I>DTI calculation with correction for motion and distortion from eddy current and B0 inhomogeneity</I><BR>" >> report.html
echo "Version $VERSION <BR>"  >> report.html
echo "Output directory:$outdir<BR>"  >> report.html
echo "<HTML><HEAD><link REL=\"stylesheet\" TYPE=\"text/css\" href=\"../.files/fsl.css\"><TITLE>DTI_PREPROCESS REPORT</TITLE></HEAD><BODY BGCOLOR="#aaaaaa">" >> $reportdir/fieldmap.html
echo "<hr><B>DTI_PREPROCESS REPORT</B><BR>"  >> $reportdir/fieldmap.html
echo "<I>DTI calculation with correction for motion and distortion from eddy current and B0 inhomogeneity</I><BR>" >>  $reportdir/fieldmap.html
echo "Version $VERSION <BR>" >>  $reportdir/fieldmap.html
echo "Output directory:$outdir<BR>" >> $reportdir/fieldmap.html

cp $reportdir/fieldmap.html $reportdir/stats.html
cp $reportdir/fieldmap.html $reportdir/motion.html
cp $reportdir/fieldmap.html $reportdir/bvecbval.html
cp $reportdir/fieldmap.html $reportdir/input.html

echo "<a href=./$reportdir/input.html>Inputs<a> - <a href=./$reportdir/fieldmap.html>Fieldmap<a> - <a href=./$reportdir/motion.html>Motion&distortion</a> - <a href=./$reportdir/bvecbval.html>  bvec&bval</a> - <a href=./$reportdir/stats.html> b0-statistics</a> - <a href=./report.html><u>Results</u></a> - <a href=./dti_preprocess.log>log</a><BR><hr>" >> report.html
echo "..not yet run" >> report.html
echo "<a href=./input.html>Inputs<a> - <a href=./fieldmap.html>Fieldmap<a> - <a href=./motion.html>Motion&distortion</a> - <a href=./bvecbval.html>bvec&bval</a> - <a href=./stats.html><u>b0-statistics</u></a> - <a href=../report.html>Results</a> - <a href=../dti_preprocess.log>log</a><BR><hr>" >> $reportdir/stats.html
echo "..not yet run" >> $reportdir/stats.html
echo "<a href=./input.html>Inputs<a> - <a href=./fieldmap.html>Fieldmap<a> - <a href=./motion.html>Motion&distortion</a> - <a href=./bvecbval.html><u>bvec&bval</u></a> - <a href=./stats.html>b0-statistics</a> - <a href=../report.html>Results</a> - <a href=../dti_preprocess.log>log</a><BR><hr>" >> $reportdir/bvecbval.html
echo "..not yet run" >> $reportdir/bvecbval.html
echo "<a href=./input.html>Inputs<a> - <a href=./fieldmap.html>Fieldmap<a> - <a href=./motion.html><u>Motion&distortion</u></a> - <a href=./bvecbval.html>bvec&bval</a> - <a href=./stats.html>b0-statistics</a> - <a href=../report.html>Results</a> - <a href=../dti_preprocess.log>log</a><BR><hr>" >> $reportdir/motion.html
echo "..not yet run" >> $reportdir/motion.html
echo "<a href=./input.html>Inputs<a> - <a href=./fieldmap.html><u>Fieldmap</u></a> - <a href=./motion.html>Motion&distortion</a> - <a href=./bvecbval.html>bvec&bval</a> - <a href=./stats.html>b0-statistics</a> - <a href=../report.html>Results</a> - <a href=../dti_preprocess.log>log</a><BR><hr>" >> $reportdir/fieldmap.html
echo "..not yet run" >> $reportdir/fieldmap.html
echo "<a href=./input.html><u>Inputs</u></a> - <a href=./fieldmap.html>Fieldmap<a> - <a href=./motion.html>Motion&distortion</a> - <a href=./bvecbval.html>bvec&bval</a> - <a href=./stats.html>b0-statistics</a> - <a href=../report.html>Results</a> - <a href=../dti_preprocess.log>log</a><BR><hr>" >> $reportdir/input.html
echo "<b>DWI 4D data :</b> $dti <BR>" >> $reportdir/input.html
echo "<b>Dwell time [msec] :</b> $esp  <BR>" >> $reportdir/input.html
echo "<b>TE [msec] :</b> $te  <BR>" >> $reportdir/input.html
echo "<b>Fieldmap :</b> $dph <BR>" >> $reportdir/input.html
echo "<b>Fieldmap magnitude :</b> $mag <BR>" >> $reportdir/input.html
echo "<b>bvecs :</b> $bvec <BR>" >> $reportdir/input.html
echo "<b>bvals :</b> $bval <BR><BR>" >> $reportdir/input.html
echo "<b>Signal loss threshold :</b> $SL [%]<BR>" >> $reportdir/input.html
echo "<b>Unwarp direction :</b> $ud <BR>" >> $reportdir/input.html


T -e "${CMD} ${inparg}"
T -e "Started at `date`"

for i in firefox mozilla netscape iexplorer opera konqueror; do sh $i $outdir/report.html & if [ $? = 0 ] ; then status=0; break; else status=1; fi; done 2> /dev/null
if [ $status = 1 ]; then T -e "You can view the results by pointing your browser at $outdir/report.html"; fi

T -e "Logfile is $outdir/$LF"
T -e "dti_preprocess : $0"
T -e "Version: $VERSION"
T -e "FSLDIR = ${FSLDIR}"
fslversion=`cat ${FSLDIR}/etc/fslversion`
T -e "FSL version = $fslversion"
if [ `echo $fslversion | head -c 1` -lt 4 ]; then
   echo "Error: dti_preprocess requires FSL 4.0 or later version."
   exit 1;
fi
T -e "Current dir = $CWD"
T -e "Output dir = $outdir"
T -e "Directory space:"
T -e "`df -h $outdir`"

if [ ${dim} != 0 ]; then
 T flirt -in $dti -ref $dti -applyisoxfm $dim -o $tmpdir/dti
 dti=$tmpdir/dti
 if [ "x${maskd}" != "x" ]; then
  T flirt -in $maskd -ref $maskd -applyisoxfm $dim -o $tmpdir/maskd
  maskd=${tmpdir}/maskd
 fi
fi

T fslroi $dti $tmpdir/ED_D_example_dti $refvol 1

if [ "x${maskd}" = "x" ]; then
 T bet $tmpdir/ED_D_example_dti $tmpdir/ED_D_example_dti_brain -f 0.2 -m
else
 T fslmaths $tmpdir/ED_D_example_dti -mas $maskd $tmpdir/ED_D_example_dti_brain
 T fslmaths $maskd $tmpdir/ED_D_example_dti_brain_mask
fi

#------------- Undistortion using fieldmap----------------#
if [ "$FC" = 0 ] ; then
 T fslmaths $tmpdir/ED_D_example_dti_brain nodif_brain
 T fslmaths $tmpdir/ED_D_example_dti_brain_mask nodif_brain_mask
else 
 T fslmaths $dph $tmpdir/FM_UD_fmap
 T fslmaths $mag $tmpdir/FM_UD_fmap_mag
 if [ "x${maskf}" = "x" ]; then T bet $tmpdir/FM_UD_fmap_mag $tmpdir/FM_UD_fmap_mag_brain -m -f 0.3
 else T fslmaths $tmpdir/FM_UD_fmap_mag -mas $maskf $tmpdir/FM_UD_fmap_mag_brain; T fslmaths $maskf -bin $tmpdir/FM_UD_fmap_mag_brain_mask
 fi
 sl=`echo "scale=3; 1 - $SL/100" | bc`;
 te=`echo "scale=3; $te * 0.001" | bc`;
 T fslmaths $tmpdir/FM_UD_fmap_mag_brain_mask -ero  $tmpdir/FM_UD_fmap_mag_brain_mask
 T fslstats $tmpdir/FM_UD_fmap -k $tmpdir/FM_UD_fmap_mag_brain_mask -P 50
 v=`fslstats $tmpdir/FM_UD_fmap -k $tmpdir/FM_UD_fmap_mag_brain_mask -P 50`
 T fslmaths $tmpdir/FM_UD_fmap -sub $v  $tmpdir/FM_UD_fmap
 T sigloss -i $tmpdir/FM_UD_fmap --te=$te -m $tmpdir/FM_UD_fmap_mag_brain_mask -s $tmpdir/FM_UD_fmap_sigloss
 T fslmaths $tmpdir/FM_UD_fmap_sigloss -mul $tmpdir/FM_UD_fmap_mag_brain $tmpdir/FM_UD_fmap_mag_brain_siglossed -odt float
 esp=`echo "scale=5; $esp/1000" | bc`
 T fugue -i $tmpdir/FM_UD_fmap_mag_brain_siglossed --loadfmap=$tmpdir/FM_UD_fmap --mask=$tmpdir/FM_UD_fmap_mag_brain_mask --dwell=$esp -w $tmpdir/FM_D_fmap_mag_brain_siglossed --nokspace --unwarpdir=$ud
 T fugue -i $tmpdir/FM_UD_fmap_sigloss --loadfmap=$tmpdir/FM_UD_fmap --mask=$tmpdir/FM_UD_fmap_mag_brain_mask --dwell=$esp -w $tmpdir/FM_D_fmap_sigloss --nokspace --unwarpdir=$ud
 T fslmaths $tmpdir/FM_D_fmap_sigloss -thr $sl $tmpdir/FM_D_fmap_sigloss
 
 if [ "$reg" != "0" ] ; then
  T flirt -in $tmpdir/ED_D_example_dti_brain -ref $tmpdir/FM_D_fmap_mag_brain_siglossed -omat $tmpdir/ED_2_FM.mat -o $tmpdir/grot -dof 6 -refweight $tmpdir/FM_D_fmap_sigloss
  T convert_xfm -omat $tmpdir/FM_2_ED.mat -inverse $tmpdir/ED_2_FM.mat
 else
  echo "1 0 0 0" >> $tmpdir/ED_2_FM.mat 
  echo "0 1 0 0" >> $tmpdir/ED_2_FM.mat 
  echo "0 0 1 0" >> $tmpdir/ED_2_FM.mat 
  echo "0 0 0 1" >> $tmpdir/ED_2_FM.mat 
  cp $tmpdir/ED_2_FM.mat $tmpdir/FM_2_ED.mat
 fi;
 
 T flirt -in $tmpdir/FM_UD_fmap -ref $tmpdir/ED_D_example_dti -init $tmpdir/FM_2_ED.mat -applyxfm -out $tmpdir/ED_UD_fmap
 T flirt -in $tmpdir/FM_UD_fmap_mag -ref $tmpdir/ED_D_example_dti -init $tmpdir/FM_2_ED.mat -applyxfm -out $tmpdir/ED_UD_fmap_mag
 T flirt -in $tmpdir/FM_UD_fmap_mag_brain -ref $tmpdir/ED_D_example_dti -init $tmpdir/FM_2_ED.mat -applyxfm -out $tmpdir/ED_UD_fmap_mag_brain
 T flirt -in $tmpdir/FM_UD_fmap_mag_brain_mask -ref $tmpdir/ED_D_example_dti -init $tmpdir/FM_2_ED.mat -applyxfm -out $tmpdir/ED_UD_fmap_mag_brain_mask
 T flirt -in $tmpdir/FM_UD_fmap_sigloss -ref $tmpdir/ED_D_example_dti -init $tmpdir/FM_2_ED.mat -applyxfm -out $tmpdir/ED_UD_fmap_sigloss
 T fslmaths $tmpdir/FM_UD_fmap_mag_brain_mask -thr 0.5 -bin $tmpdir/FM_UD_fmap_mag_brain_mask -odt float
 T fslmaths $tmpdir/ED_UD_fmap_sigloss -thr $sl $tmpdir/ED_UD_fmap_sigloss -odt float

 T fugue --loadfmap=$tmpdir/ED_UD_fmap --dwell=$esp -i $tmpdir/ED_D_example_dti -u $tmpdir/ED_UD_example_dti --unwarpdir=$ud --saveshift=$tmpdir/ED_UD_shift #--mask=$tmpdir/ED_UD_fmap_mag_mask
 T convertwarp -s $tmpdir/ED_UD_shift -o $tmpdir/ED_UD_warp -r $tmpdir/ED_D_example_dti --shiftdir=$ud
 T fslmaths $tmpdir/ED_UD_example_dti $outdir/nodif
 T applywarp -i $tmpdir/ED_D_example_dti_brain -o $outdir/nodif_brain -w $tmpdir/ED_UD_warp -r $tmpdir/ED_D_example_dti --abs #--mask=$tmpdir/ED_UD_fmap_mag_mask
 T fslmaths $outdir/nodif_brain -bin  $outdir/nodif_brain_mask
 T fslmaths $tmpdir/ED_UD_fmap_sigloss nodif_sigloss

 cat $reportdir/fieldmap.html | head -n 6 > $reportdir/tmpfieldmap.html
 mv $reportdir/tmpfieldmap.html $reportdir/fieldmap.html
 T fslstats $tmpdir/FM_UD_fmap -R
 T fslmaths $tmpdir/FM_UD_fmap -mas $tmpdir/FM_UD_fmap_mag_brain_mask $tmpdir/FM_UD_fmap_brain
 v=`fslstats $tmpdir/FM_UD_fmap -R | awk '{print $1}'`
 V=`fslstats $tmpdir/FM_UD_fmap_brain -P 5 -P 95`
 T fslmaths $tmpdir/FM_UD_fmap -sub $v -add 10 -mas $tmpdir/FM_UD_fmap_mag_brain_mask grot
 T fslstats grot -P 5 -P 95
 v=`fslstats grot -P 5 -P 95`
 T overlay 0 0 $tmpdir/FM_UD_fmap_mag -a grot $v grot
 T slicer grot -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/fmap+mag.gif
 T /bin/cp ${FSLDIR}/etc/luts/ramp.gif $reportdir/ramp.gif
 T /bin/cp ${FSLDIR}/etc/luts/ramp2.gif $reportdir/ramp2.gif
 echo "<B>Fieldmap overlaid on magnitude image </B><BR>`echo $V | awk '{print $1}'`<IMG src=./ramp.gif width=106 height=14 border=0 align=middle >`echo $V | awk '{print $2}'` [rad/sec]<BR>"  >> $reportdir/fieldmap.html
 echo "<IMG src= "./fmap+mag.gif" width="1200" height="100" border="0"><BR><BR>" >> $reportdir/fieldmap.html

 T fslstats $tmpdir/ED_UD_shift -R -P 1 -P 99
 O=`fslstats $tmpdir/ED_UD_shift -R -P 1 -P 99`
 p=`echo $O | awk '{print $3}'`
 q=`echo $O | awk '{print $4}'`
 p=`echo "scale=1; $p * -1" | bc`
 T fslmaths $tmpdir/ED_UD_shift -mul -1 grot
 T overlay 1 0 $tmpdir/ED_UD_fmap_mag_brain -a $tmpdir/ED_UD_shift 0 $q grot 0 $p grot
 T slicer grot -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_UD_shift+mag.png
 echo "<B>Unwarping shift map in voxels </B><BR>-$p <IMG src=./ramp2.gif width=106 height=14 border=0 align=middle> 0 <IMG src=./ramp.gif width=106 height=14 border=0 align= middle >$q<BR>"  >> $reportdir/fieldmap.html
 echo "<IMG src="./ED_UD_shift+mag.png" width="1200" height="100" border="0"><BR><BR>" >> $reportdir/fieldmap.html 

 T flirt -in $tmpdir/FM_D_fmap_mag_brain_siglossed -ref $tmpdir/ED_D_example_dti_brain -applyxfm -init $tmpdir/FM_2_ED.mat -o grot
 T slicer grot -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/FM_D_fmap_mag_brain_siglossed.gif
 T slicer $tmpdir/ED_D_example_dti_brain  -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65  sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_D_example_dti_brain.gif
 T whirlgif -o $reportdir/ED_FM_D_movie2.gif -loop -time 50 $reportdir/ED_D_example_dti_brain.gif $reportdir/FM_D_fmap_mag_brain_siglossed.gif
 echo "<B>Registration of brain images between original b=0 and estimated from fieldmap </B><BR> "  >> $reportdir/fieldmap.html
 echo "<IMG src= "./ED_FM_D_movie2.gif" width="1200" height="100" border="0"><BR><BR>" >> $reportdir/fieldmap.html

 T slicer $tmpdir/ED_D_example_dti  -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_D_example_dti.gif
 T slicer $tmpdir/ED_UD_example_dti   -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_UD_example_dti.gif
 T whirlgif -o $reportdir/ED_UD_movie2.gif -loop -time 50 $reportdir/ED_D_example_dti.gif $reportdir/ED_UD_example_dti.gif

 cat report.html | head -n 6 > tmpreport.html; mv tmpreport.html report.html
 echo "<B>b=0 image</B><BR>"  >> report.html
 echo "Corrected b=0 image: nodif<BR>"  >> report.html
 echo "Corrected b=0 brain image: nodif_brain<BR>" >>  report.html
 echo "Signal loss image: nodif_sigloss<BR><BR>"  >> report.html
 echo "<B>Uncorrected and corrected b=0 images</B><BR>"  >> report.html
 echo "<IMG src="./${reportdir}/ED_UD_movie2.gif" width="1200" height="100" border="0"><BR><BR>" >>  report.html
 T fslstats $tmpdir/ED_UD_fmap_mag_brain.nii.gz -P 20 -P 90
 v=`fslstats $tmpdir/ED_UD_fmap_mag_brain.nii.gz -P 20 -P 90`
 T slicer $tmpdir/ED_UD_fmap_mag -s 3 -i $v -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_UD_fmap_mag.gif
 T whirlgif -o $reportdir/ED_UD_movie3.gif -loop -time 50  $reportdir/ED_D_example_dti.gif $reportdir/ED_UD_example_dti.gif $reportdir/ED_UD_fmap_mag.gif 
 echo "<B>Uncorrected, corrected b=0 images and a fieldmap magnitude image </B><BR>"  >> report.html
 echo "<IMG src=./${reportdir}/ED_UD_movie3.gif width=1200 height=100 border=0><BR><BR>" >> report.html

 T overlay 1 0 nodif -a nodif_sigloss $sl 1 grot
 T slicer grot -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/nodif+sigloss.png
 echo "<B>Corrected b0 image and signal loss estimated from fieldmap</B><BR>"  >> report.html
 echo "<IMG src="./${reportdir}/nodif+sigloss.png" width="1200" height="100" border="0"><BR><BR>" >> report.html
fi


#------------ Correction for eddy current distortion and motion---------------#
if [ $test != 1 ]; then
 if [ $FC != 0 ] ; then
 echo "Finished fieldmap distortion correction and now performing corrention for motion and eddy-current distortion."
 echo "Check fieldmap distortion correction by pointing browser at $reportdir/fieldmap.html"
 echo ""
 fi
 
 if [ $FLIRT = 0 ]; then
   if [ -e $tmpdir/dti_ecc.mat ]; then \rm -Rf $tmpdir/dti_ecc.mat;fi
   T mcflirt -in $dti -o ${tmpdir}/dti_ecc -report -stages $MCFLIRT -mats -dof 12 -rmsabs -rmsrel -reffile $tmpdir/ED_D_example_dti; 
 else
  mkdir -p $tmpdir/dti_ecc.mat
  i=0
  while [ $i -lt $dtidim4 ]; do
  j=`zeropad $i 4`
  T fslroi $dti $tmpdir/vol_$j $i 1
  T flirt -in $tmpdir/vol_$j -ref $tmpdir/ED_D_example_dti -nosearch -paddingsize 1 -omat $tmpdir/dti_ecc.mat/MAT_$j
  i=`expr $i + 1`
  done  
 fi

 transposematrix $bval $tmpdir/ts_bval
 transposematrix $bvec $tmpdir/ts_bvec
 for i in translation.par rotation.par scale.par skew.par ts_bvec_ecc; do if [ -e ${tmpdir}/${i} ] ; then rm -Rf ${tmpdir}/${i}; fi; touch $tmpdir/${i};done
 zero=`round 0`; pi=$(echo "scale=20; 4*a(1)" | bc -l)
 
 i=0
 while [ $i -lt $dtidim4 ]; do
  j=`zeropad $i 4`
  v=`avscale --allparams $tmpdir/dti_ecc.mat/MAT_$j | head -13 | tail -7 | cut -d "=" -f 2 | grep [0-9]`
  echo $v | awk '{print $4,$5,$6}' >> $tmpdir/translation.par
  echo $v | awk '{printf "%f %f %f\n",$1*180/3.141592653,$2*180/3.141592653,$3*180/3.141592653}' >> $tmpdir/rotation.par
  echo $v | awk '{print $7,$8,$9}' >> $tmpdir/scale.par
  echo $v | awk '{print $10,$11,$12}' >> $tmpdir/skew.par
  if [ "$bvec_ecc" != 0 ]; then
    k=`expr $i + 1`;
    xR=`echo $v | awk '{print $1}'`; yR=`echo $v | awk '{print $2}'`; zR=`echo $v | awk '{print $3}'`;
    xUF=`cat ${tmpdir}/ts_bvec | awk 'NR=='$k' {print $1}'`; yUF=`cat ${tmpdir}/ts_bvec | awk 'NR=='$k' {print $2}'`; zUF=`cat ${tmpdir}/ts_bvec | awk 'NR=='$k' {print $3}'`
    x=`round $xUF`; y=`round $yUF`; z=`round $zUF`
    if [ $x = "$zero" ] && [ $y = "$zero" ] && [ $z = "$zero" ] ; then echo "$x $y $z" >> $tmpdir/ts_bvec_ecc; else
        #-----------Rotation around x-------------#
        MagnitudeX=`echo "sqrt ( $y ^ 2 + $z ^2 )" | bc -l` ;
        if [ $z != "$zero" ]; then DirectionX=`echo "a ( $y / $z )" | bc -l`; else DirectionX=`echo "$pi / 2" | bc -l`; fi;
        FactorX=0;
        ReverseX=`echo "$DirectionX + $FactorX" | bc -l` ;
        AxUF="$x";AyUF=`echo "$MagnitudeX * ( s ($ReverseX) )" | bc -l` ; AzUF=`echo "$MagnitudeX * ( c ($ReverseX) )" | bc -l` ;
        Ax=`round $AxUF`;Ay=`round $AyUF`; Az=`round $AzUF`;
        Rx=`round $x`; Ry=`round $y`; Rz=`round $z`;
        if [ $Rx = $Ax ] && [ $Ry = $Ay ] && [ $Rz = $Az ] ; then FactorX=$FactorX ; else FactorX=$pi; fi;
        NewDirectionX=`echo "$DirectionX + $FactorX + $xR" | bc -l` ;
        x=$x; y=`echo "$MagnitudeX * ( s ($NewDirectionX) )" | bc -l` ; z=`echo "$MagnitudeX * ( c ($NewDirectionX) )" | bc -l` ;
        #-----------Rotation around y-------------#
        MagnitudeY=`echo "sqrt ( $x ^ 2 + $z ^2 )" | bc -l` ;
        if [ $x != "$zero" ]; then DirectionY=`echo "a ( $z / $x )" | bc -l` ; else DirectionY=`echo "$pi / 2" | bc -l`; fi;
        FactorY=0;
        ReverseY=`echo "$DirectionY + $FactorY" | bc -l` ;
        AxUF=`echo "$MagnitudeY * ( c ($ReverseY) )" | bc -l`; AyUF=$y; AzUF=`echo "$MagnitudeY * ( s ($ReverseY) )" | bc -l` ;
        Ax=`round $AxUF`; Ay=`round $AyUF`; Az=`round $AzUF`;
        Rx=`round $x`; Ry=`round $y`; Rz=`round $z`;
        if [ $Rx = $Ax ] && [ $Ry = $Ay ] && [ $Rz = $Az ] ; then FactorY=$FactorY ; else FactorY=$pi; fi;
        NewDirectionY=`echo "$DirectionY + $FactorY + $yR" | bc -l` ;
        x=`echo "$MagnitudeY * ( c ($NewDirectionY) )" | bc -l`; y=$y; z=`echo "$MagnitudeY * ( s ($NewDirectionY) )" | bc -l` ;
        #-----------Rotation around z-------------#
        MagnitudeZ=`echo "sqrt ( $x ^ 2 + $y ^2 )" | bc -l` ;
        if [ $x != "$zero" ] ; then DirectionZ=`echo "a ( $y / $x )" | bc -l` ; else DirectionZ=`echo "$pi / 2" | bc -l`; fi;
        FactorZ=0;
        ReverseZ=`echo "$DirectionZ + $FactorZ" | bc -l` ;
        AxUF=`echo "$MagnitudeZ * ( c ($ReverseZ) )" | bc -l`;AyUF=`echo "$MagnitudeZ * ( s ($ReverseZ) )" | bc -l` ; AzUF=$z;
        Ax=`round $AxUF`;Ay=`round $AyUF`;Az=`round $AzUF`;
        Rx=`round $x`; Ry=`round $y`; Rz=`round $z`;
        if [ $Rx = $Ax ] && [ $Ry = $Ay ] && [ $Rz = $Az ] ; then FactorZ=$FactorZ ; else FactorZ=$pi; fi;
        NewDirectionZ=`echo "$DirectionZ + $FactorZ + $zR" | bc -l`;
        x=`echo "$MagnitudeZ * ( c ($NewDirectionZ) )" | bc -l`; y=`echo "$MagnitudeZ * ( s ($NewDirectionZ) )" | bc -l`; z=$z;
        echo "`round $x` `round $y` `round $z`" >> $tmpdir/ts_bvec_ecc
    fi;
  fi;
  i=`expr $i + 1`
 done
 T fsl_tsplot -i $tmpdir/translation.par -o $reportdir/translation.png -t "Translation" -y [mm] -x Volume -a x,y,z
 T fsl_tsplot -i $tmpdir/rotation.par -o $reportdir/rotation.png -t "Rotation" -y [degree] -x Volume -a x,y,z

 cat $reportdir/motion.html | head -n 6 > $reportdir/tmpmotion.html
 mv $reportdir/tmpmotion.html $reportdir/motion.html
 echo "<B>Estimated motion</B><BR>" >> $reportdir/motion.html
 meandisl=`cat $tmpdir/translation.par | awk 'BEGIN {x=0;N=0};{x=x+($1^2+$2^2+$3^2)^0.5;N=N+1}END {printf("%f",x/N)}'`
 echo "Mean dislocation : $meandisl [mm] <BR>"  >> $reportdir/motion.html
 meansdt=`threecolumnmeansd $tmpdir/translation.par`
 meansdr=`threecolumnmeansd $tmpdir/rotation.par`
 echo "Mean translation: (x y z)=(`echo $meansdt | awk '{print $1,$2,$3}'`) [mm]<BR>"  >> $reportdir/motion.html
 echo "Mean rotation: (x y z)=(`echo $meansdr | awk '{print $1,$2,$3}'`) [degree]<BR><BR>"  >> $reportdir/motion.html
 echo "<IMG src=../$reportdir/translation.png><BR><BR>" >> $reportdir/motion.html
 echo "<IMG src=../$reportdir/rotation.png><BR><BR>" >> $reportdir/motion.html
 T fsl_tsplot -i $tmpdir/scale.par -o $reportdir/scale.png -t "Scale"  -x Volume -a x,y,z
 T fsl_tsplot -i $tmpdir/skew.par -o $reportdir/skew.png -t "Skew" -x Volume -a x,y,z
 meansdS=`threecolumnmeansd $tmpdir/scale.par`
 meansds=`threecolumnmeansd $tmpdir/skew.par`
 echo "<B>Estimated distortion</B><BR>" >> $reportdir/motion.html
 echo "Mean scale: (x y z)=(`echo $meansdS | awk '{print $1,$2,$3}'`) <BR>"  >> $reportdir/motion.html
 echo "Mean skew: (x y z)=(`echo $meansds | awk '{print $1,$2,$3}'`) <BR><BR>"  >> $reportdir/motion.html
 echo "<IMG src=../$reportdir/scale.png><BR><BR>" >> $reportdir/motion.html
 echo "<IMG src=../$reportdir/skew.png><BR><BR>" >> $reportdir/motion.html
 
 #bvecs&bvals
 T fslstats $tmpdir/ED_D_example_dti -x
 v=`fslstats $tmpdir/ED_D_example_dti -x`
 T tsplot $outdir -f $dti -C $v $tmpdir/ts_maxvoxel.txt
 T fsl_tsplot -i $tmpdir/ts_maxvoxel.txt,$tmpdir/ts_bval -o $reportdir/ts_signal-bval.png -a Signal,b-value -t ts-signal-bval-plot -x Volume
 cat $reportdir/bvecbval.html | head -n 6 > $reportdir/tmpbvecbval.html
 mv $reportdir/tmpbvecbval.html $reportdir/bvecbval.html
 echo "<B>Time seris plot of signal and b-value </B><BR>"  >> $reportdir/bvecbval.html
 echo "Time series signal was taken from voxel at $v<BR>" >> $reportdir/bvecbval.html 
 echo "<IMG src=./ts_signal-bval.png><BR><BR>" >> $reportdir/bvecbval.html
 T fsl_tsplot -i $tmpdir/ts_bvec -o $reportdir/ts_bvec.png -a x,y,z -t ts-bvec-plot -x Volume
 cat ${tmpdir}/ts_bvec | awk '{print -$1,-$2,-$3"\n"$1,$2,$3""}' > tmp_bvec
 statsbvec=`threecolumnmeansd tmp_bvec`
 /bin/rm tmp_bvec 
 echo "<B>Statistics of b-vector</B><BR>"  >> $reportdir/bvecbval.html
 echo "Center of gravity: (x,y,z)=(`echo $statsbvec | awk '{print $1,$2,$3}'`)<BR>" >> $reportdir/bvecbval.html
 echo "Standard deviation: (x,y,z)=(`echo $statsbvec | awk '{print $4,$5,$6}'`)<BR><BR>" >> $reportdir/bvecbval.html
 echo "<B>Time seris plot of b-vector </B><BR>"  >> $reportdir/bvecbval.html
 echo "<IMG src=./ts_bvec.png><BR><BR>" >> $reportdir/bvecbval.html
 if [ "$bvec_ecc" != 0 ]; then 
  T fsl_tsplot -i $tmpdir/ts_bvec_ecc -o $reportdir/ts_bvec_ecc.png -a x,y,z -t ts-motion-correced-bvec-plot -x Volume
  cat ${tmpdir}/ts_bvec_ecc | awk '{print -$1,-$2,-$3"\n"$1,$2,$3""}' > tmp_bvec
  statsbvec=`threecolumnmeansd tmp_bvec`
  /bin/rm tmp_bvec
  echo "<B>Statistics of motion-corrected b-vector</B><BR>"  >> $reportdir/bvecbval.html
  echo "Center of gravity: (x,y,z)=(`echo $statsbvec | awk '{print $1,$2,$3}'`)<BR>" >> $reportdir/bvecbval.html
  echo "Standard deviation: (x,y,z)=(`echo $statsbvec | awk '{print $4,$5,$6}'`)<BR><BR>" >> $reportdir/bvecbval.html
  echo "<B>Time seris plot of motion-corrected b-vector </B><BR>"  >> $reportdir/bvecbval.html
  echo "<IMG src=./ts_bvec_ecc.png><BR><BR>" >> $reportdir/bvecbval.html
 fi;

 #Reslice images & calculate DTI
 if [ $FC != 0 ] ; then
  i=0
  while [ $i -lt $dtidim4 ]; do
  j=`zeropad $i 4`
  if [ $FLIRT = 0 ]; then T fslroi $dti $tmpdir/vol_$j $i 1; fi
   T applywarp -i $tmpdir/vol_$j -o $tmpdir/vol_$j --premat=$tmpdir/dti_ecc.mat/MAT_$j -w $tmpdir/ED_UD_warp -r $tmpdir/ED_D_example_dti --abs #--mask=$tmpdir/ED_UD_fmap_mag_brain_mask
  i=`expr $i + 1`
  done;
 T fslmerge -t data `imglob -oneperimage $tmpdir/vol_????.*`
 T imrm $tmpdir/vol_????*
 else 
  if [ $FLIRT = 0 ] ; then
   T fslmaths $tmpdir/dti_ecc data
  else
   i=0
   while [ $i -lt $dtidim4 ]; do
   j=`zeropad $i 4`
   flirt -in $tmpdir/vol_$j -o $tmpdir/vol_$j -applyxfm -init $tmpdir/dti_ecc.mat/MAT_$j -ref $tmpdir/ED_D_example_dti
   i=`expr $i + 1`
   done
  T fslmerge -t data `imglob -oneperimage $tmpdir/vol_????.*`
  T imrm $tmpdir/vol_????*
  fi
 fi
 T cp -f $bvec bvecs
 T cp -f $bval bvals

 if [ "$unmaskd" != 1 ] ; then
  T dtifit -k $dti -o $tmpdir/dti_D -b bvals -r bvecs -m $tmpdir/ED_D_example_dti_brain_mask --sse
 else
  T fslmaths $tmpdir/ED_D_example_dti_brain_mask -add 1 -bin ${tmpdir}/mask
  T dtifit -k $dti -o $tmpdir/dti_D -b bvals -r bvecs -m $tmpdir/mask --sse
 fi
 T slicer $tmpdir/dti_D_FA -s 3 -i 0 1.0 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ${reportdir}/dti_FA_D.gif

 if [ "$unmaskd" != 1 ] ; then
  T dtifit -k data -o dti -b bvals -r bvecs -m nodif_brain_mask --sse
 else
  T dtifit -k data -o dti -b bvals -r bvecs -m $tmpdir/mask --sse
 fi

 T slicer dti_FA -s 3 -i 0 1.0 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/dti_FA_UD.gif
 T whirlgif -o $reportdir/dti_FA_movie2.gif -loop -time 50 $reportdir/dti_FA_D.gif $reportdir/dti_FA_UD.gif

 T slicer $tmpdir/ED_UD_fmap_mag -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/ED_UD_fmap_mag.gif 

 if [ `cat report.html | wc -l` -lt 8 ] ;then cat report.html | head -n 6 > tmpreport.html; mv tmpreport.html report.html; fi

 echo "<B>DTI image</B><BR>"  >> report.html
 echo "Motion&distortion-corrected DTI images: dti_FA, dti_V1-3, dti_L1-3,dti_S0, dti_MD, dti_M0<BR>"  >> report.html
 if [ "$bvec_ecc" != 0 ]; then echo "bvec-corrected DTI images, as well as bvecs are saved in the directory: dti_bvecmc<BR>"  >> report.html; fi; echo "<BR>" >> report.html
 echo "<B>Uncorrected and motion&distortion-corrected FA images </B><BR>"  >> report.html
 echo "<IMG src=./$reportdir/dti_FA_movie2.gif width=1200 height=100 border=0><BR><BR>" >> report.html
 if [ $FC != 0 ] ; then
  T whirlgif -o $reportdir/mag+dti_FA_movie3.gif -loop -time 50 $reportdir/dti_FA_D.gif  $reportdir/dti_FA_UD.gif $reportdir/ED_UD_fmap_mag.gif 
  echo "<B>Uncorrected, motion&distortion-corrected FA images and a fieldmap magnitude image </B><BR>"  >> report.html
  echo "<IMG src=./$reportdir/mag+dti_FA_movie3.gif width=1200 height=100 border=0><BR><BR>" >> report.html
 fi

 if [ "$bvec_ecc" != 0 ]; then
  mkdir -p dti_bvecmc
  transposematrix $tmpdir/ts_bvec_ecc dti_bvecmc/bvecs;
  T /bin/cp bvals dti_bvecmc/
  if [ "$unmaskd" != 1 ] ; then
   T dtifit -k data -o dti_bvecmc/dti -b dti_bvecmc/bvals -r dti_bvecmc/bvecs -m nodif_brain_mask --sse
  else
   T dtifit -k data -o dti_bvecmc/dti -b dti_bvecmc/bvals -r dti_bvecmc/bvecs -m ${tmpdir}/mask --sse
  fi
  T slicer dti_bvecmc/dti_FA -s 3 -i 0 1.0 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
  T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/dti_bvecsecc_FA_UD.gif
  T whirlgif -o $reportdir/dti_FA_movie3.gif -loop -time 50 $reportdir/dti_FA_D.gif $reportdir/dti_FA_UD.gif $reportdir/dti_bvecsecc_FA_UD.gif
  echo "<B>Uncorrected, motion&distortion-corrected and bvecs-corrected FA images </B><BR>"  >> report.html
  echo "<IMG src=./$reportdir/dti_FA_movie3.gif width=1200 height=100 border=0><BR><BR>" >> report.html
 fi

 #Statistic of b=0 images
 for i in `cat $tmpdir/ts_bval  | awk '{if ($1==0) print NR-1}'`; do 
  j=`zeropad $i 4`
  T fslroi ${dti} ${tmpdir}/nodif_${j}_orig $i 1
  if [ $FLIRT != 1 ] ; then T fslroi ${tmpdir}/dti_ecc ${tmpdir}/nodif_${j}_ecc $i 1; else T flirt -in ${tmpdir}/nodif_${j}_orig -ref $tmpdir/ED_D_example_dti -applyxfm -init $tmpdir/dti_ecc.mat/MAT_$j -o ${tmpdir}/nodif_${j}_ecc; fi;
  T fslroi data ${tmpdir}/nodif_${j}_eccUD $i 1
 done
 T fslmerge -t ${tmpdir}/nodifmerge_orig `imglob -onerperimage ${tmpdir}/nodif_????_orig.nii.gz`
 T fslmaths ${tmpdir}/nodifmerge_orig -Tmean ${tmpdir}/nodifmean_orig -odt float
 T fslmaths ${tmpdir}/nodifmerge_orig -Tstd ${tmpdir}/nodifstd_orig -odt float
 T fslmaths ${tmpdir}/nodifstd_orig -div ${tmpdir}/nodifmean_orig -mas $tmpdir/ED_D_example_dti_brain_mask -mul 100 ${tmpdir}/nodifcov_orig -odt float
 v=`fslstats $tmpdir/nodifcov_orig -P 10 -P 90`
 T slicer $tmpdir/nodifcov_orig -s 3 -i $v -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/nodifcov_orig.gif
 cat $reportdir/stats.html | head -n 6 > $reportdir/tmpstats.html
 mv $reportdir/tmpstats.html $reportdir/stats.html
 echo "<B>COV of b=0 images (original)</B><BR>"  >> $reportdir/stats.html
 cov=`fslstats $tmpdir/nodifcov_orig -k $tmpdir/ED_D_example_dti_brain_mask -m`
 echo "Averaged COV in brain region = $cov[%]<BR>" >> $reportdir/stats.html
 echo "<IMG src=./nodifcov_orig.gif  width=1200 height=100 border=0><BR><BR>" >> $reportdir/stats.html

 T fslmerge -t ${tmpdir}/nodifmerge_ecc `imglob -onerperimage ${tmpdir}/nodif_????_ecc.nii.gz`
 T fslmaths ${tmpdir}/nodifmerge_ecc -Tmean ${tmpdir}/nodifmean_ecc -odt float
 T fslmaths ${tmpdir}/nodifmerge_ecc -Tstd ${tmpdir}/nodifstd_ecc -odt float
 T fslmaths ${tmpdir}/nodifstd_ecc -div ${tmpdir}/nodifmean_ecc -mas $tmpdir/ED_D_example_dti_brain_mask -mul 100 ${tmpdir}/nodifcov_ecc -odt float
 T slicer $tmpdir/nodifcov_ecc -s 3 -i $v -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
 T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/nodifcov_ecc.gif
 echo "<B>COV of b=0 images (corrected for eddy-current distortion and motion)</B><BR>"  >> $reportdir/stats.html
 cov=`fslstats $tmpdir/nodifcov_ecc -k $tmpdir/ED_D_example_dti_brain_mask -m`
 echo "Averaged COV in brain region = $cov[%]<BR>" >> $reportdir/stats.html
 echo "<IMG src=./nodifcov_ecc.gif  width=1200 height=100 border=0><BR><BR>" >> $reportdir/stats.html

 if [ $FC != 0 ] ; then
  T fslmerge -t ${tmpdir}/nodifmerge_eccUD `imglob -onerperimage ${tmpdir}/nodif_????_eccUD.nii.gz`
  T fslmaths ${tmpdir}/nodifmerge_eccUD -Tmean ${tmpdir}/nodifmean_eccUD -odt float
  T fslmaths ${tmpdir}/nodifmerge_eccUD -Tstd ${tmpdir}/nodifstd_eccUD -odt float
  T fslmaths ${tmpdir}/nodifstd_eccUD -div ${tmpdir}/nodifmean_eccUD -mas $tmpdir/ED_D_example_dti_brain_mask -mul 100 ${tmpdir}/nodifcov_eccUD -odt float
  T slicer $tmpdir/nodifcov_eccUD -s 3 -i $v -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
  T pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png $reportdir/nodifcov_eccUD.gif
  echo "<B>COV of b=0 images (corrected for eddy-current distortion, motion, B0 inmogeneity distortion)</B><BR>"  >> $reportdir/stats.html
  cov=`fslstats $tmpdir/nodifcov_eccUD -k nodif_brain_mask -m`
  echo "Averaged COV in brain region = $cov[%]<BR>" >> $reportdir/stats.html
  echo "<IMG src=./nodifcov_eccUD.gif width=1200 height=100 border=0><BR><BR>" >> $reportdir/stats.html
 fi
fi

#---------------Finish-------------------#
T /bin/rm -f sla.png slb.png slc.png sld.png sle.png slf.png slg.png slh.png sli.png slj.png slk.png sll.png grot*
T imrm $tmpdir/grot $tmpdir/dti_ecc $tmpdir/nodif_* $tmpdir/nodifmerge* $tmpdir/nodifstd*
echo "<BR></P><hr><I>Created by dti_preprocess written by Takuya Hayashi (takuya.hayashi@gmail.com)</I></BODY></HTML>" >> report.html
T -e "Finished $CMD at `date`"
cd $CWD
exit 0;

