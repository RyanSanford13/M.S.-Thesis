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
Associated Data Files
	<br> &emsp; 192ir-hdr_gammamed_plus.xls

Detailed descriptions of each of the major projects can be found below.


## DICOM File Conversion 
Description: Code used to convert a DICOM RT Structure (RT) File into a DICOM Surface Segmentation (SEG) File. User provides an image set and the RT file. The script extracts the relevant information from the RT file needed to process the contour information and generate a triangulated mesh. Using the DICOM library in MATLAB, the mesh is written to the SEG file. 
<br> _convertRTStructToSeg.m_
<br> inputs: 
	<br> &emsp; *RT filename*.dcm - RT file exported from treatment planning system; linked to volume MR image set 
 	<br> &emsp; *image filepath* - path to folder containing associated volume MR image set that the RT file was created on 
<br> output:
	<br> &emsp; *SEG filename.dcm* - converted SEG file.


## Needle Path Optimization 
Description: Collection of scripts used in the two wrappers - *committee_sphere.m* and *committee_anatomy.m* - to facillitate optimized selection of needle paths for use in needle path planning in HDR GYN brachytherapy. User employs scripts to prepare data for optimization process and the wrapper scripts make use of support functions to accomplish repeated processes. 
<br>
<br> **Data preparation**
<br> _organMesh.m_
<br> Description: Creates triangulated meshes of anatomy from provided DICOM RT file. Anatomy serves as base geometry that can be optimized over for needle path optimization
<br> inputs: 
	<br> &emsp; *RT filename.dcm* - RT file exported from treatment planning system; Contains anatomy contours; linked to volume  MR image set
 	<br> &emsp; *image filepath* - path to folder containing associated volume MR image set that the RT file was created on
  <br> output:
  	<br> &emsp; *organ_name.mat* - structure containing the organ mesh (faces and vertices) and the binary mask used in the triangulation process
<br> functions:
	<br> &emsp; *createBinaryMask.m* 
 	<br> &emsp; Description: Creates a 3D binary mask by concatenating 2D binary masks based on contour data from the provided MR volume image <br> &emsp; set and RT file
<br>
<br> _Consensus_data.m_
<br> Description: Processes HDR source consensus data for use in TG-43 dose calculations in the needle path optimization process
<br> inputs: 
	<br> &emsp; *192ir-hdr_gammamed_plus.xls* - spreadsheet containing HDR source consensus data 
 <br> output:
 	<br> &emsp; *consensus data varaible name.mat* - structure containing the consensus data for the associated variable
<br>
<br> **Needle Path Optimization**
<br> _committee_sphere.m_
<br> Description: Wrapper script to perform needle path optimization over a sphere. 
<br> functions:
	<br> &emsp; *randTri.m*
 	<br> &emsp; &emsp; Description: Selects a triangle from a provided mesh based on probability proportional to area; takes the triangulated entry and/or  <br> &emsp; &emsp; target zone(s) and selects a triangle based on a probability distribution defined by the area of the triangles within the mesh 
  	<br> &emsp; *uniformSample.m*
   	<br> &emsp; &emsp; Description: Randomly selects a point within a provided polygon; takes the triangle chosen from the randTri() function and randomly <br> &emsp; &emsp; selects a point within the boundaries of the triangle
    	<br> &emsp; *checkIntersect.m*
     	<br> &emsp; &emsp; Description: determines if a line intersects a trinagulated mesh; determines if the needle path, represented as a line, intersects the <br> &emsp; &emsp; sphere target (sphere triangulated mesh)
	<br> &emsp; *generatePointsAlongLine.m*
 	<br> &emsp; &emsp; Description: generates points along a line a uniform distances; generates the dwell locations along the path of the needle path, 
  	<br> &emsp; *isPointInSphere.m*
   	<br> &emsp; &emsp; Description: determines if a point lies withing a sphere mesh; determines if a dwell position is located within the spherical target
    	<br> &emsp; *greedy_needle_selection.m*
     	<br> &emsp; &emsp; Description: See entry below; Optimization script used to determine the optimal needle paths over a target

<br> _committee_anatomy.m_
<br> Description: Wrapper script to perform needle path optimization over a provided anatomical geometry 
<br> inputs:
	<br> &emsp; uterus.mat
 	<br> &emsp; cerv_hrctv.mat
  	<br> &emsp; rectum.mat
   	<br> &emsp; vagina.mat
    	<br> &emsp; hrctv_marg.mat 
<br> functions: 
	<br> &emsp; _randTri.m_
 	<br> &emsp; &emsp; Description: Selects a triangle from a provided mesh based on probability proportional to area; takes the triangulated entry and/or  <br> &emsp; &emsp; target zone(s) and selects a triangle based on a probability distribution defined by the area of the triangles within the mesh 
	<br> &emsp; _uniformSample.m_
   	<br> &emsp; &emsp; Description: Randomly selects a point within a provided polygon; takes the triangle chosen from the randTri() function and randomly <br> &emsp; &emsp; selects a point within the boundaries of the triangle



