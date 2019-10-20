# 1 comp number
rm *.nii.gz
for jj in 1
do
   /home/jingleil/bin/Map_reference2brain /home/jingleil/ldrive/Lab_ChristineG/JingleiL/EEG-fMRI/MNI152_3mm_mask.nii.gz P${jj}.txt $1 ./P${jj}test
	for ((i=1;i<=$1;i++)) do
	   echo $i
	   flirt -applyisoxfm 2 -in ./P${jj}test_map$i.nii -ref ./P${jj}test_map$i.nii -out ./P${jj}test_map${i}_2mm.nii
           fdr -i ./P${jj}test_map${i}_2mm.nii -m /home/jingleil/ldrive/Lab_ChristineG/JingleiL/EEG-fMRI/MNI152_2mm_mask_c.nii.gz -q 0.05 -a P${jj}test_map${i}_2mm_a.nii
           fslmaths ./P${jj}test_map${i}_2mm_a.nii -ptoz ./P${jj}test_map${i}_2mm_z.nii
	   fslmaths ./P${jj}test_map${i}_2mm_z.nii -mul -1 ./P${jj}test_map${i}_2mm_minusz.nii
	   
	   overlay 1 0 /home/jingleil/ldrive/Lab_ChristineG/JingleiL/EEG-fMRI/MNI152_T1_2mm_brain.nii -a ./P${jj}test_map${i}_2mm_z.nii 2.0 5 ./P${jj}test_map${i}_2mm_minusz.nii 2.0 5 ./temp.nii.gz
	   fslroi ./temp.nii.gz ./temp.nii.gz 0 90 0 108 4 72
	   slicer ./temp.nii.gz -S 2 860 ./P${jj}test_map$i.png
	done
done

