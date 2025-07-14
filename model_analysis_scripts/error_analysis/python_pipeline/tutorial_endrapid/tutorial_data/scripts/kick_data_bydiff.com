#! /bin/tcsh -f
#
#	change FP by RMS (FP - FC)				-James Holton 6-24-12
#
#


set mtzfile = refined.mtz
set F = ""
set SIGF = ""
set FC = ""
set seed = `date +%N | awk '{print $1/1000}'`


set tempfile = ${CCP4_SCR}/kick_data$$


foreach arg ( $* )

    if( "$arg" =~ *.mtz ) set mtzfile  = "$arg"
    if(( "$arg" =~ *[0-9] )&&( "$arg" =~ [1-9]* )) set steps = "$arg"

    if( "$arg" =~ seed=* ) then
        set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
        if("$test" != "") set seed = $test
    endif
    if( "$arg" =~ F=* ) then
        set test = `echo $arg | awk -F "[=]" '{print $2}'`
        if("$test" != "") set user_F = $test
    endif
    if( "$arg" =~ FC=* ) then
        set test = `echo $arg | awk -F "[=]" '{print $2}'`
        if("$test" != "") set user_FC = $test
    endif
end

# examine MTZ file
echo "go" | mtzdump hklin $mtzfile |\
awk '/OVERALL FILE STATISTICS/,/No. of reflections used/' |\
awk 'NF>5 && $(NF-1) ~ /^[FJDGKQLMPWABYIRUV]$/' |\
cat >! ${tempfile}mtzdmp

# use completeness, or F/sigF to pick default F
cat ${tempfile}mtzdmp |\
awk '$(NF-1) == "F"{F=$NF; meanF=$8; reso=$(NF-2); comp=substr($0,32)+0; \
  getline; if($(NF-1) != "Q") next; \
  S=$NF; if($8) meanF /= $8; print F, S, reso, comp, meanF;}' |\
sort -k3n,4 -k4nr,5 -k5nr >! ${tempfile}F

cat ${tempfile}mtzdmp |\
awk '$(NF-1) == "F"{F=$NF; meanF=$8; reso=$(NF-2); comp=substr($0,32)+0; \
  getline; if($(NF-1) != "P") next; \
  PHI=$NF; if($8) meanF /= $8; print F, PHI, reso, comp, meanF;}' |\
cat >! ${tempfile}FC

# and extract all dataset types/labels
cat ${tempfile}mtzdmp |\
awk 'NF>2{print $(NF-1), $NF, " "}' |\
cat >! ${tempfile}cards

#clean up
rm -f ${tempfile}mtzdmp

if("$F" == "" || "$SIGF" == "") then
    # pick F with best resolution, or F/sigma
    if($?user_F) then
	set F = `awk -v F=$user_F '$1==F{print}' ${tempfile}F`
	if("$F" == "") then
	    echo "WARNING: $user_F not found in $mtzfile"
	    cat ${tempfile}cards
	endif
    endif
    if($#F < 2) set F = `head -1 ${tempfile}F`
    if($#F > 2) then
	set SIGF = $F[2]
	set F    = $F[1]
	echo "selected F=$F SIGF=$SIGF "
    endif
    rm -f ${tempfile}F
endif
if("$F" == "") then
    set BAD = "no Fobs in $mtzfile "
    goto exit
endif
if("$FC" == "") then
    # pick FC with best resolution?  Or just first one in the file.
    if($?user_FC) then
	set FC = `awk -v FC=$user_FC '$1==FC{print}' ${tempfile}FC`
	if("$FC" == "") set FC = `awk -v FC=$user_FC '$2==FC{print $2}' ${tempfile}cards`
	if("$FC" == "") then
	    echo "WARNING: $user_FC not found in $mtzfile"
	    cat ${tempfile}cards
	endif
    endif
    if("$FC" == "") set FC = `head -1 ${tempfile}FC`
    if($#FC > 1) then
	set FC    = $FC[1]
	echo "selected FC=$FC "
    endif
    rm -f ${tempfile}FC
endif
if("$FC" == "") then
    set BAD = "no Fcalc in $mtzfile "
    goto exit
else
    echo "using FC=$FC"
endif
# capture remaining columns for imort later...
set otherstuff = `awk -v F=$F -v SIGF=$SIGF -v FC=$FC '$2!=F && $2!=SIGF && $2!=FC{++n;print "E"n"="$2}' ${tempfile}cards`
rm -f ${tempfile}cards
#echo "$otherstuff"


# extract only F/SIGF so sftools does not get pissed off
cad hklin1 $mtzfile hklout F_FC.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF E3=$FC
EOF

echo "using seed = $seed"
sftools << EOF >&! ${tempfile}sftools.log
read F_FC.mtz
calc seed $seed
calc col delta = col $F col $FC -
calc col noise = col delta RAN_G *
calc F col Fnew  = col $F col noise +
select col Fnew < 0
calc col Fnew = 0
select all
calc F col $F = col Fnew
delete  col Fnew noise delta
write new.mtz 
y
exit
y
EOF
if($status) then
    set BAD = "sftools failed"
    goto exit
endif

# one last clean-up
mv new.mtz cleanme.mtz
sftools << EOF > /dev/null
read cleanme.mtz
absent col $SIGF if col $F ABSENT
write new.mtz
EOF
cad hklin1 new.mtz hklout cleaned.mtz << EOF > /dev/null
labin file 1 all
EOF


# now put other stuff back in
if("$otherstuff" != "") then
    cad hklin1 cleaned.mtz hklin2 $mtzfile hklout kicked.mtz << EOF > /dev/null
labin file 1 all
labin file 2 $otherstuff
EOF
else
    cad hklin1 cleaned.mtz hklout kicked.mtz << EOF > /dev/null
labin file 1 all
EOF
endif
if($status) then
    set BAD = "output mtz corrupted."
    goto exit
endif
rm -f F_FC.mtz
rm -f new.mtz
rm -f cleaned.mtz

echo "kicked.mtz contains $F from $mtzfile modified by rms ( $F - $FC )"

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

exit


######################################################################3
#
#	notes and stuff
#

../kick_data_bydiff.com ${ID}.updated_refine_001_f_model.mtz seed=1

cad hklin1 ${ID}.updated_refine_001_f_model.mtz hklin2 kicked.mtz hklout test.mtz << EOF
labin file 1 E1=FOBS E2=FMODEL
labin file 2 E1=FOBS 
labou file 2 E1=Fnew
EOF

mtzdmp test.mtz -1 |\
awk 'substr($0,1,13)==sprintf(" %4d%4d%4d",$1,$2,$3){print $4,$5,$6}' |\
cat >! plotme.txt





