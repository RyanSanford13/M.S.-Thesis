%% Organ Mesh Generation Script %%

rtStructInfo = dicominfo('9_6_24_MR_Struct_ang/RS.qCRQMi6g4AYjASZTXUYnXMY3z.9424_axt2 3d.dcm');
imageFilePath = '9_6_24_MR_Struct_ang';

 if(strcmp(rtStructInfo.SOPClassUID,'1.2.840.10008.5.1.4.1.1.481.3')) 
        disp('Starting the conversion...');
        
        % Process depends on homogenous coordinate systems and matrix / vector algebra
        % For a specifc contour set, collect points from each contour ROI into a single '3D contour'
        % Find the number of ROIs that need to be converted 
        nROI = length(fieldnames(rtStructInfo.ROIContourSequence));
        disp(nROI);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GET VOLUME INFO FOR IMAGe AND STRUCTURE BEFORE LOOPING THROUGH CONTOURS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Load MR files
        filesArray = string({dir(imageFilePath + "/MR*").name});
        numFiles = length(filesArray);

        % Sort MR files by number once read inclo
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
        % Defines translation and is constant for any given volume + image combo 
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

        % Below info doesnt change between MR images
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
        save('in2di_rotMat.mat', "in2di_rotationMat");
        
        % Define blank 3d mask
        mask3d = zeros(image_width, image_height, numFiles);

        for i = 1:nROI
            roiIndex = sprintf('Item_%d', i);

            % Get ROI name
            roi_name = getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName;
        
            % Skip 'Body' ROI or other undesired ROIs if needed
            if strcmp(roi_name, 'BODY')
                continue;
            end
        
            % Check if processing the second-to-last ROI (Outer Vagina)
            if i == nROI - 3
                fprintf('Processing combined binary mask for Outer and Inner Vagina...\n');
        
                % Create binary masks for Outer Vagina and Inner Vagina
                outerMask = createBinaryMask(rtStructInfo, i, numFiles, image_height, image_width, filesArray, imageFilePath, in2di_rotationMat);
                innerMask = createBinaryMask(rtStructInfo, i + 1, numFiles, image_height, image_width, filesArray, imageFilePath, in2di_rotationMat);
        
                % % Initialize combined mask
                % combinedMask = false(size(innerMask));
                % 
                % % Find slices where both inner and outer masks exist
                % sharedSlices = find(any(innerMask, [1, 2]) & any(outerMask, [1, 2]));
                % 
                % % Handle shared slices
                % for sliceIdx = sharedSlices
                %     innerSlice = innerMask(:, :, sliceIdx);
                %     outerSlice = outerMask(:, :, sliceIdx);
                % 
                %     % Ensure innerSlice is concentric within outerSlice
                %     combinedMask(:, :, sliceIdx) = ~innerSlice & outerSlice;
                % end
                % 
                % % Append extra inner vagina slices
                % extraInnerSlices = find(~any(innerMask, [1, 2]) & any(outerMask, [1, 2]));
                % for sliceIdx = extraInnerSlices
                %     combinedMask(:, :, sliceIdx) = outerMask(:, :, sliceIdx);
                % end
                combinedMask = outerMask & ~innerMask;
               
        
                % Extract the isosurface from the combined binary mask
                [faces, verts] = extractIsosurface(combinedMask, 0.5);
        
                % Convert vertices to DICOM coordinates
                verts = verts.';
                verts(4, :) = 1; % Homogeneous dimension
                verts_di = in2di_rotationMat * verts;
        
                % Save the combined structure
                vagina = struct();
                vagina.verts = verts_di;
                vagina.faces = faces;

                vagina.mask = combinedMask;
        
                save('vagina.mat', 'vagina');
        
                % Break the loop early since the final item is already processed
                %break;
            end
        
            % For all other ROIs, process normally
            fprintf('Processing ROI: %s\n', roi_name);
            
            % if ~strcmp(roi_name, 'Cervix + HRCTV')
                % Create binary mask for the current ROI
                mask3d = createBinaryMask(rtStructInfo, i, numFiles, image_height, image_width, filesArray, imageFilePath, in2di_rotationMat);
                if strcmp(roi_name, 'Cervix + HRCTV')
                    for i = 36:57 % contour range
                        if ~any(mask3d(:, :, i))
                            % Shift the remaining slices forward
                            mask3d(:, :, i:end-1) = mask3d(:, :, i+1:end);
                            
                            % Pad the last slice with empty values
                            mask3d(:, :, end) = false;
                        end
                    end
                end
                % Extract the isosurface from the binary mask
                [faces, verts] = extractIsosurface(mask3d, 0.5);
            
                % Convert vertices to DICOM coordinates
                verts = verts.';
                verts(4, :) = 1; % Homogeneous dimension
                verts_di = in2di_rotationMat * verts;
            % else
            %     [mask3d_cerv, mask3d_hrctv] = createBinaryMask(rtStructInfo, i, numFiles, image_height, image_width, filesArray, imageFilePath, in2di_rotationMat);
            % 
            %     % cervix
            %     [faces_cerv, verts_cerv] = extractIsosurface(mask3d_cerv, 0.5);
            % 
            %     % Convert vertices to DICOM coordinates
            %     verts_cerv = verts_cerv.';
            %     verts_cerv(4, :) = 1; % Homogeneous dimension
            %     verts_cerv_di = in2di_rotationMat * verts_cerv;
            % 
            %     % hrctv
            %     [faces_hrctv, verts_hrctv] = extractIsosurface(mask3d_hrctv, 0.5);
            % 
            %     % Convert vertices to DICOM coordinates
            %     verts_hrctv = verts_hrctv.';
            %     verts_hrctv(4, :) = 1; % Homogeneous dimension
            %     verts_hrctv_di = in2di_rotationMat * verts_hrctv;
            % end
        
            switch roi_name
                %
                case 'Cervix + HRCTV'
                    
                    cerv_hrctv = struct();
                    cerv_hrctv.verts = verts_di;
                    cerv_hrctv.faces = faces;
                    cerv_hrctv.mask = mask3d;
                    
                    % % cervix
                    % cerv_hrctv.cerv.verts = verts_cerv_di;
                    % cerv_hrctv.cerv.faces = faces_cerv;
                    % cerv_hrctv.cerv.mask = mask3d_cerv;
                    % 
                    % % hrctv
                    % cerv_hrctv.hrctv.verts = verts_hrctv_di;
                    % cerv_hrctv.hrctv.faces = faces_hrctv;
                    % cerv_hrctv.hrctv.mask = mask3d_hrctv;

                    
        
                    save('cerv_hrctv.mat', 'cerv_hrctv');
        
                case 'Uterus'
        
                    uterus = struct();
                    uterus.verts = verts_di;
                    uterus.faces = faces;
                    uterus.mask = mask3d;
        
                    save('uterus.mat', 'uterus');
        
                case 'Rectum_ang'
        
                    rectum = struct();
                    rectum.verts = verts_di;
                    rectum.faces = faces;
                    rectum.mask = mask3d;
        
                    save('rectum.mat', 'rectum');
                case 'Cervix + HRCTV + margin'
                    hrctv_marg = struct();
                    hrctv_marg.verts = verts_di;
                    hrctv_marg.faces = faces;
                    hrctv_marg.mask = mask3d;

                    save('hrctv_marg.mat', 'hrctv_marg');
            end
        end 
 end 

