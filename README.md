# M.S.-Thesis
Public repository containing associated MATLAB scripts and relevant data files used in the completion of the thesis titled "Needle Optimization for Use in HDR GYN Brachytherapy" (Ryan Sanford, Duke University, 2025).

These scripts were created solely by Ryan Sanford in the completion of his M.S. Thesis. 

Copyright (c) 2025 Ryan Sanford

List of Major Projects and Associated Scripts
	<br /> &emsp; DICOM File Conversion 
		<br /> &emsp; &emsp; convertRTStructToSeg.m
	<br /> &emsp; Needle Path Optimization 
    	<br /> &emsp; &emsp; comittee_sphere.m
    	<br /> &emsp; &emsp; comittee_anatomy.m
    	<br /> &emsp; &emsp; optimize_dwell_times.m
    	<br /> &emsp; &emsp; dose_calc.m
    	<br /> &emsp; &emsp; dose_eval.m
    	<br /> &emsp; &emsp; greedy_needle_selection.m 
    	<br /> &emsp; &emsp; organ_mesh.m
    	<br /> &emsp; &emsp; createBinaryMask.m
    	<br /> &emsp; &emsp; randTri.m
    	<br /> &emsp; &emsp; uniformSample.m 
    	<br /> &emsp; &emsp; checkIntersect.m
    	<br /> &emsp; &emsp; generatePointsAlongLine.m
    	<br /> &emsp; &emsp; isPointInSphere.m
     	<br /> &emsp; &emsp; Consensus_data.m

Detailed descriptions of each of the major projects can be found below.


## DICOM File Conversion 
Description: Code used to convert a DICOM RT Structure (RT) File into a DICOM Surface Segmentation (SEG) File. User provides an image set and the RT file. The script extracts the relevant information from the RT file needed to 
<br> input: 
	<br> &emsp; *RT filename*.dcm - RT file exported from treatment planning system; linked to volume MR image set 
 	<br> &emsp; *image filepath* - path to folder containing associated volume MR image set that the RT file was created on 
<br> 

## Needle Path Optimization 
Description: Collection of scripts used in the two wrappers - *committee_sphere.m* and *committee_anatomy.m* - to facillitate optimized selection of needle paths for use in GYN HDR brachytherapy. *committee_sphere.m* requires 

