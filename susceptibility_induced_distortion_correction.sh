#unwarping data... 

ww=CrM_aOnly
mkdir applyB
cd applyB

run=(1 2 3 4 5 6 7 8 9 10)
for mm in "${run[@]}"
do
3dvolreg -verbose -zpad 1 -base CrM4.nii[0] \
	-1Dfile dfile.r${mm}.1D -prefix volreg.r${mm} \
	-cubic \
	-1Dmatrix_save mat.r${mm}.vr.aff12.1D \
	CrM${mm}.nii
done

fslroi CrM4.nii AP.nii 0 -1 0 -1 0 -1 0 1; fslroi phase_s1.nii PA.nii 0 -1 0 -1 0 -1 0 1
fslmerge -a mergedAP-PA AP.nii PA.nii

#make slice number even (75 to 74) - see fsl topup.wiki!! 
3dZcutup -keep 0 73 -overwrite -prefix mergedAP-PA.nii mergedAP-PA.nii; 

printf "0 -1 0 0.026\n0 1 0 0.026\n" > acqparams.txt
cat acqparams.txt;

topup --imain=mergedAP-PA.nii --datain=acqparams.txt --config=b02b0.cnf --out=topup_base --iout=b0_unwarped --fout=fieldmap_Hz 


#apply B0 to all 
run=( 1 2 3 4 5 6 7 8 9 10 ) # 
for mm in "${run[@]}"
do
3dZcutup -keep 0 73 -overwrite -prefix volreg.CrM${mm}.nii volreg.r${mm}+orig.
	echo 'Run No.' ${mm}
	applytopup --imain=volreg.CrM${mm}.nii --topup=topup_base --method=jac --datain=acqparams.txt --inindex=1 --out=top.CrM${mm}
echo 'Run no.'${mm} 'finished. Removing tmp_CrM'${mm}' files..'
3dcopy top.CrM${mm}.nii.gz top.CrM${mm}
rm -rf top.CrM${mm}.nii.gz
done


cp ../pb01.${ww}.r04.volreg+orig.* .
cp ../pb01.${ww}.r11.volreg+orig.* .
3dcopy pb01.CrM_main.r04.volreg+orig.HEAD CrM_main.r04.nii
3dcopy pb01.CrM_main.r11.volreg+orig.HEAD CrM_main.rphase.nii
#cut volume 0-1 
fslroi CrM_main.r04.nii b0_CrM.nii 0 -1 0 -1 0 -1 0 1; fslroi CrM_main.rphase.nii b0_phase.nii 0 -1 0 -1 0 -1 0 1
#make slice number even (75 to 74) - see fsl topup.wiki!! 
3dZcutup -keep 0 73 -overwrite -prefix b0_CrM.nii b0_CrM.nii; 
3dZcutup -keep 0 73 -overwrite -prefix b0_phase.nii b0_phase.nii; 
rm -f b0_*.nii.gz

#merge (A>P first. See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/Faq) 
fslmerge -a merged b0_CrM.nii b0_phase.nii
 # To compute the etl (echo train length), run:
   # fslhd merged.nii | grep desc
printf "0 -1 0 0.026\n0 1 0 0.026\n" > acqparams.txt
cat acqparams.txt;

#calculate... 
#acqparams.txt contains how I acquired fMRI data 
#(by column, x y z (roughly, TR/slice N. which is total readout time for each slice) )
topup --imain=merged.nii --datain=acqparams.txt --config=b02b0.cnf --out=topup_base --iout=b0_unwarped --fout=fieldmap_Hz 

#apply B0 to all 
run=( 01 02 03 04 05 06 07 08 09 10 ) # 
for mm in "${run[@]}"
do
3dZcutup -keep 0 73 -overwrite -prefix tmp_CrM${mm}.nii ../pb01.${ww}.r${mm}.volreg+orig.HEAD
	echo 'Run No.' ${mm}
	applytopup --imain=tmp_CrM${mm}.nii --topup=topup_base --method=jac --datain=acqparams.txt --inindex=1 --out=top.CrM${mm}
echo 'Run no.'${mm} 'finished. Removing tmp_CrM'${mm}' files..'
3dcopy top.CrM${mm}.nii.gz top.CrM_main.r${mm}
rm -rf tmp_CrM${mm}.nii; top.CrM${mm}.nii.gz
done

folder=(CrM_main)  #CrM_aOnly

for ww in "${folder[@]}"
do

if [ "$ww" == "CrM_main" ]; then
run=( 01 02 03 04 05 06 07 08 09 10 )
elif [ "$ww" == "CrM_aOnly" ]; then
run=( 01 02 03 04 05 06 07 08 09 10 )
fi
echo $ww
echo "${run[*]}"

#gedit /home/jiyeongha/Script/phaseCorr.sh &

	##------------------------02 Scale & remove skull in EPI -----------

	for r in "${run[@]}"
	do
	3dClipLevel top.${ww}.r${r}+orig.HEAD >> clip.txt
	3dTstat -mean -prefix r.${r}.base top.${ww}.r${r}+orig.HEAD'[0..$]'
	done
	more clip.txt # check the smallest clip value across all runs
	clip=$(cut -f1 -d"," clip.txt | sort -n | head -1) # assign the clip value

	for r in "${run[@]}"
	do
	3dcalc -a top.${ww}.r${r}+orig. -b r.${r}.base+orig. \
	       -expr "(100 * a/b) * step(b-$clip)" -prefix top.pb02.${ww}.r${r}.scaled #remove the space and scale 
	done
	# (tip) when calling variables(like $clip) within 3dcalc -expr, be sure to use "", and not ''

	##------------------------ 03 Detrending & 04 highpass filter & 05 add mean (100)-----
	for i in "${run[@]}"
	do

	3dDetrend -polort 1 -prefix top.pb03.${ww}.r${i}.sc_dt top.pb02.${ww}.r${i}.scaled+orig.
	3dBandpass -prefix top.pb04.${ww}.r${i}.sc_dt_hp 0.01 99999 top.pb03.${ww}.r${i}.sc_dt+orig
	#add mean of scaled data (pb02.~)
	3dTstat -mean -prefix r.${i}.sc_base top.pb02.${ww}.r${i}.scaled+orig.HEAD'[0..$]'
	3dcalc  -a top.pb04.${ww}.r${i}.sc_dt_hp+orig.HEAD -b r.${i}.sc_base+orig.HEAD \
	        -expr 'a+b' -prefix top.pb05.${ww}.r${i}.sc_dt_hp_am
	done
	##Finished. Let's combine all data-------------------------------------
	3dTcat -prefix top.${SN}${ww}FNL.nii \
	top.pb05.${ww}.*.HEAD 
	mv top.${SN}${ww}FNL.nii ../../../Decoding_${SN}/
done














#apply B0 to all 
#I'll use eddy instead of applytopup -> impossible since we don't have bvec & bval (only DWI data can have them)
#ref: https://lcni.uoregon.edu/kb-articles/preprocessing-of-diffusion-data-with-multiple-phase-encoding-#directions


#make a mask that separate brain from non-brain
fslmaths b0_unwarped -Tmean b0_unwarped
bet b0_unwarped b0_unwarped_brain -m

indx=""
for ((i=1; i<=170; i+=1)); do indx="$indx 1"; done
echo $indx > index.txt

rcnt=1  
run=( 1 2 3 4 5 6 7 8 9 10 ) # 
for mm in "${run[@]}"
do
3dZcutup -keep 0 73 -overwrite -prefix tmp_CrM${mm}.nii ../CrM${mm}.nii
fslsplit tmp_CrM${mm}.nii
	for fileName in ./vol*.nii.gz
	do
	echo $fileName
	f=${fileName:4:3}
	echo $f
	applytopup --imain=$fileName,b0_phase.nii --topup=topup_base --datain=acqparams.txt --inindex=1,2 --out=top_r${rcnt}_$f
	done
3dTcat -prefix top_CrM$(($rcnt+1)).nii top_r$rcnt_*
rm -rf vol*.nii.gz; rm -rf tmp_CrM${mm}.nii
rcnt=$(($rcnt+1))
done














i=0;
for i in {0..168};
do
   printf "0 -1 0 0.026\n" >> acqparams.txt
done 


run=(1)  #2 3 4 5 6 7 8 9 10

rcnt=1 
for rawfile in "${run[@]}"
do
echo "split" 
echo $rcnt
fslsplit ../CrM${rawfile}.nii
	for fileName in ./vol*.nii.gz
	do
	echo $fileName
	f=${fileName:4:3}
	echo $f
3dZcutup -keep 0 73 -overwrite -prefix $fileName $fileName; 
	applytopup --imain=b0_phase.nii,$fileName --topup=topup_base --datain=acqparams.txt --inindex=1,2 --out=top_r${rcnt}_$f
	done
3dTcat -prefix top_CrM$rcnt.nii top_r$rcnt_*
rm -rf vol*.nii.gz
rcnt=$(($rcnt+1))
done
