function segfile = convertRTStructToSeg(rtStructFile, imageFilePath)
% Reads and processes provided DICOM RT file and converts + writes to a DICOM SEG file 
% note: template seg file provided with code and needed to run script 

% rtStructFile - struct.dcm; rtStruct file that provides the ROI information for conversion
% imageFilePath - string; path on local machine to the associated folder containing all MR images for the rtStruct file 
% use: segFle = convertRTStructToSeg(rtStructFile, ImageFilePath)  
    %%
    % Confirm correct number of input arguments
    if nargin ~= 2
        disp("Incorrect number of Arguments: " + nargin);
        disp("To use: convertRTStructToSeg(<rt_struct_filename>, <string_image_file_path>)");
        return
    end

    % Switch dicom dictionary to include Eigen private tags
    dicomdict('set', 'new-dicom-dict.txt');
    dicomdict('get');

    % Initialize headers for manipulation 
    seg_info = dicominfo('SEG1');
    new_headr = seg_info;

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Data Formatting %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    % read RT file and extract relevant information 

    rtStructInfo = dicominfo(rtStructFile);
    % Confirm file is an rtStruct

    if(strcmp(rtStructInfo.SOPClassUID,'1.2.840.10008.5.1.4.1.1.481.3')) 
        disp('Starting the conversion...');
        
        % Process depends on homogenous coordinate systems and matrix / vector algebra
        % For a specifc contour set, collect points from each contour ROI into a single '3D contour'
        % Find the number of ROIs that need to be converted 
        nROI = length(fieldnames(rtStructInfo.ROIContourSequence));
        disp(nROI);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GET VOLUME INFO FOR IMAGE AND STRUCTURE BEFORE LOOPING THROUGH CONTOURS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Load MR files
        filesArray = string({dir(imageFilePath + "/MR*.dcm").name});
        numFiles = length(filesArray);

        % Sort MR files by number once read in
        numArray = [];
        for a = 1:numFiles
            number = extractAfter(filesArray(a), " ");
            number = str2double(regexp(number, '\d+', 'match'));
            numArray(a) = number;
        end

        % Sort extracted file numbers
        [~, idx] = sort(numArray);
        filesArray = filesArray(idx);
        
        % Find lowest slice position
        % Defines translation and is constant for any given volume image set 
        % Initialize position
        ppz = 0;
        patientPos = 0;
        for m = 1 :length(filesArray)
            tempInfo = dicominfo(imageFilePath + "/" + filesArray(m) + "");
            
            if(abs(tempInfo.SliceLocation) > abs(ppz))
                ppz = tempInfo.SliceLocation;
                patientPos = tempInfo.ImagePositionPatient;
            end
        end

        % Below info doesnt constant for entire image set 
        % Choose first MR Image to get MR info
        mrInfo = dicominfo(imageFilePath + "/" +filesArray(1));
        image_height = double(mrInfo.Height);
        image_width = double(mrInfo.Width);
        
        % Define dicom 'unit' vectors
        xd = [mrInfo.ImageOrientationPatient(1), mrInfo.ImageOrientationPatient(2), mrInfo.ImageOrientationPatient(3)];
        yd = [mrInfo.ImageOrientationPatient(4), mrInfo.ImageOrientationPatient(5), mrInfo.ImageOrientationPatient(6)];
        zd = cross(xd, yd);

        % Get scaling factors per slice
        dx = mrInfo.PixelSpacing(1,1);
        dy = mrInfo.PixelSpacing(2,1);
        dz = mrInfo.SliceThickness;

        % Create homogeneous transformation matrix
        % This matrix converts from index to dicom coordinates 
        in2di_rotationMat = [([xd.*dx, patientPos(1,1)]); ([yd.*dy, patientPos(2,1)]); ([zd.*dz, patientPos(3,1)]); [0 0 0 1]];
        
        % Define blank 3d mask
        mask3d = zeros(image_width, image_height, numFiles);

        % Loop through ROIs
        for i = 1:nROI
            

            roiIndex = sprintf('Item_%d', i);
            % Each ROI is called Item_#
            % roiIdx is a string, "Item_i", where i varies in loop
            % For Anatomical Phantom:
                % Item_1 = Uterus ROI
                % Item_2 = CTV-HR ROI
            % For Deformable Registration Phantom:
                % Item_1 = Uterus ROI
                % Item_2 - Item_7 = BB ROI's (for registration)
            roi = getfield(rtStructInfo.ROIContourSequence, roiIndex);
            roi_name = getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName;
            disp(roi_name);

            if ~strcmp(roi_name, 'HRCTV') && ~strcmp(roi_name, 'Cervix + HRCTV')
                continue;
            end

            nConts = length(fieldnames(roi.ContourSequence));
            disp(nConts);

            % Loop through the contours in the ROI 
            for j = 1:nConts
                
                % Each Contour is called Item_#
                % itemIdx is a string, "Item_j", where j varies in loop 
                % next_itemIdx is for the next contour layer 
                itemIdx = sprintf('Item_%d', j);
    
                % Reshape data such that each point is a column vector in the larger point matrix
                layer = getfield(roi.ContourSequence, itemIdx);
                di_points = reshape(layer.ContourData, 3, []);
    
                % Sort the layer you are working with based on x-y polar angle 
                for k = 1:length(di_points)
                    ang = atan2(di_points(2,k), di_points(1,k));
                    di_points(4,k) = ang;
                end
    
                di_points = sort(di_points, 4, 'ascend');
                
                % Once sorted, set the fourth point (your homogeneous dimmension) to a filler variable, in this case 1
                di_points(4,:) = 1;
    
                % Find matching MR file
                for l = 1:numFiles
                    % Get current file 
                    tempInfo = dicominfo(imageFilePath + "/" +filesArray(l) + "");
    
                    % Find matching instance ID's
                    if(strcmp(string(tempInfo.SOPInstanceUID), string(getfield(roi.ContourSequence, itemIdx).ContourImageSequence.Item_1.ReferencedSOPInstanceUID)))
                        % Enter this portion if theres a match
                        disp('Matchng MR file found...');
                        
                        % Convert dicom points to index
                        % in2di converts index to dicom -> to go from dicom to index multiply by the inverse
                        % in2di_rotationMat\ = inv(in2di_rotationMat)
                        % Find slice location based on rotation matrix and ImagePositionPatient 
                        imgIdx = in2di_rotationMat\[tempInfo.ImagePositionPatient; 1];
                        
                        in_points = in2di_rotationMat\di_points;
                        
                        % Create binary mask of contour
                        bw = poly2mask(in_points(1,:), in_points(2,:), image_height, image_width);
                        
                        % Add to mask3d at the file location
                        % Image was scanned in reverse, therefore add slices in reverse order
                        % 96 - l + 1, since 1 based indexing system, if zero based index then 96 - l
                        % 96 is num slices
                        mask3d(:,:,round(abs(imgIdx(3))+1)) = bw;
                        disp("Match at slice" + l + " added to mask...")

                        % Qualitative check to confirm the conversion us proceeding correctly 
                        figure(1)
                        subplot(121)
                            imshow(mask3d(:,:,round(abs(imgIdx(3)) + 1)))
                            hold on
                            plot(in_points(1,:), in_points(2,:), 'b', 'LineWidth',2)
                            hold off
                        subplot(122)
                            imshow(dicomread(tempInfo), [])
                            hold on 
                            plot(in_points(1,:), in_points(2,:), 'b', 'LineWidth',2)
                            hold off
                    else
                        continue;
                    end 
                end   
            end

            % Convert binary mask into surface mesh
            % Extract triangle correspondences and points for seg file 
            [faces, verts] = extractIsosurface(mask3d, .5);
            
            % Faces = Triangle Correspondences || Verts = Coordinate data 
            % Format vertices to be used with the homogeneous matrix 
            verts = verts.';
            verts(4,:) = 1;
    
            % verts in index coordinates, convert to dicom 
            verts_di = in2di_rotationMat * verts;
            
            % display dicom and index verts
            % Qualitative check to confirm generated meshes are the same 
            figure
            subplot(121)
                title('Index Verts');
                trisurf(faces, verts(1,:), verts(2,:), verts(3,:), image_height, 'facealpha', 0.3)
            subplot(122)
                title('Dicom Verts')
                trisurf(faces, verts_di(1,:), verts_di(2,:), verts_di(3,:), image_height, 'facealpha', 0.3)
            
            % volumeViewer(mask3d)

            % get ROI name for seg file assignment 
            roi_name = getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName;
            roi_in_header = struct();
            
            % Add mesh information to the SEG header 
            if strcmp(roi_name, 'Uterus') || strcmp(roi_name, 'uterus')
                    % loop rois
                    roi_in_header = getfield(new_headr.SurfaceSequence, 'Item_1');
            
                    % Change triangle correspondences 
                    faces = uint16(reshape(faces', [], 1));
                    % C++ indexing begins at 0
                    roi_in_header.SurfaceMeshPrimitivesSequence.Item_1.TrianglePointIndexList = faces-1;
            
                    % Change coordinate point data
                    num_points = length(verts_di(1,:));
                    verts_di = reshape(verts_di(1:3, :), 1, []);
                    roi_in_header.SurfacePointsSequence.Item_1.PointCoordinatesData = verts_di;
                    roi_in_header.SurfacePointsSequence.Item_1.NumberOfSurfacePoints = num_points;

                    new_headr.SurfaceSequence.Item_1 = roi_in_header;

            elseif strcmp(roi_name, 'HRCTV') || strcmp(roi_name, 'hrctv') || strcmp(roi_name, 'Cervix + HRCTV + margin')
                    % loop rois
                    roi_in_header = getfield(new_headr.SurfaceSequence, 'Item_2');
            
                    % Change triangle correspondences 
                    faces = uint16(reshape(faces', [], 1));
                    % C++ indexing begins at 0
                    roi_in_header.SurfaceMeshPrimitivesSequence.Item_1.TrianglePointIndexList = faces-1;
            
                    % Change coordinate point data
                    % Volume came from aria hardcoded rn, discuss with RV on how to change 
                    volume = 16.69;
                    num_points = length(verts_di(1,:));
                    verts_di = reshape(verts_di(1:3, :), 1, []);
                    roi_in_header.SurfacePointsSequence.Item_1.PointCoordinatesData = verts_di;
                    roi_in_header.SurfacePointsSequence.Item_1.NumberOfSurfacePoints = num_points;
                    roi_in_header.WsVolume = volume;

                    new_headr.SurfaceSequence.Item_2 = roi_in_header;

            else
                    % if i > 2 && strcmp(roi_name, 'Uterus') ~= 1 && strcmp(roi_name, 'HRCTV') ~=1
                    %     item_number = sprintf('Item_%d', i);
                    %     % loop rois
                    %     roi_in_header = getfield(new_headr.SurfaceSequence, item_number);
                    % 
                    %     % Change triangle correspondences 
                    %     faces = uint16(reshape(faces, [], 1));
                    %     % C++ indexing begins at 0
                    %     roi_in_header.SurfaceMeshPrimitivesSequence.Item_1.TrianglePointIndexList = faces-1;
                    % 
                    %     % Change coordinate point data
                    %     num_points = length(verts_di(1, :));
                    %     verts_di = reshape(verts_di(1:3, :), 1, []);
                    %     roi_in_header.SurfacePointsSequence.Item_1.PointCoordinatesData = verts_di;
                    %     roi_in_header.SurfacePointsSequence.Item_1.NumberOfSurfacePoints = num_points;
                    % 
                    %     new_headr.SurfaceSequence = setfield(new_headr.SurfaceSequence, item_number, roi_in_header);
                    % end
                    continue;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DIRCOM header preparation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Change filename 
    filename = "RSconvertedToSEG_MR_Images.dcm";
    new_headr.Filename = filename;
    
    % Change file mod date
    new_headr.FileModDate = datetime('now');

    % Match tags from MR Images to SEG
    % These tages are constant across all slices=, take info from first MR
    mrMatchInfo = dicominfo(imageFilePath + "/" + filesArray(1));
        % Name
        new_headr.PatientName.FamilyName = mrMatchInfo.PatientName.FamilyName;
        
        % Birthdate
        new_headr.PatientBirthDate = mrMatchInfo.PatientBirthDate;

        % Patient ID
        new_headr.PatientID = mrMatchInfo.PatientID;

        % Accension Number
        new_headr.AccessionNumber = mrMatchInfo.AccessionNumber;

        % Study UID
        new_headr.StudyInstanceUID = mrMatchInfo.StudyInstanceUID;

        % Study Date
        new_headr.StudyDate = mrMatchInfo.StudyDate;


        % Referenced Image Sequence 
        for j = 1:numFiles
            % Get Item_i index strings
            segItemIdx = sprintf('Item_%d', j);
            
            % Get MR Files 
            tempInfo = dicominfo(imageFilePath + "/" + filesArray(j) + "");

            new_headr = setfield(new_headr, 'ReferencedImageSequence', segItemIdx, 'ReferencedSOPClassUID', tempInfo.SOPClassUID);
            new_headr = setfield(new_headr, 'ReferencedImageSequence', segItemIdx, 'ReferencedSOPInstanceUID', tempInfo.SOPInstanceUID);
        end


    % Write dicom file 
    dicomwrite([], filename, new_headr, 'CreateMode', 'copy', 'WritePrivate', true);
    fprintf('\t RT STruct written to: ' + filename + '\n');

    % Switch dictionary back to factory settings 
    dicomdict('factory');
    dicomdict('get');
    
end











