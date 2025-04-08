function varargout = createBinaryMask(rtStructInfo, roiIdx, numFiles, image_height, image_width, filesArray, imageFilePath, in2di_rotationMat)
    mask3d = false(image_height, image_width, numFiles); % Initialize 3D mask
    
    mask3d_cerv = false(image_height, image_width, numFiles);
    mask3d_hrctv = false(image_height, image_width, numFiles);

    roiIndex = sprintf('Item_%d', roiIdx);
    roi = getfield(rtStructInfo.ROIContourSequence, roiIndex);
    nConts = length(fieldnames(roi.ContourSequence));

    for j = 1:nConts
        itemIdx = sprintf('Item_%d', j);
        layer = getfield(roi.ContourSequence, itemIdx);
        di_points = reshape(layer.ContourData, 3, []);
        
        for k = 1:length(di_points)
            ang = atan2(di_points(2, k), di_points(1, k));
            di_points(4, k) = ang;
        end

        di_points = sort(di_points, 4, 'ascend');
        di_points(4, :) = 1;

        for l = 1:numFiles
            tempInfo = dicominfo(imageFilePath + "/" + filesArray(l) + "");

            if strcmp(string(tempInfo.SOPInstanceUID), ...
                      string(getfield(roi.ContourSequence, itemIdx).ContourImageSequence.Item_1.ReferencedSOPInstanceUID))
                imgIdx = in2di_rotationMat \ [tempInfo.ImagePositionPatient; 1];
                in_points = in2di_rotationMat \ di_points;

                % Special handling for the last four contours (concentric layers)
                if strcmp(getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName, 'Outer Vagina') && j >= 10
                    if j == 10 || j == 12
                        % Outer contour for the slice
                        bw_outer = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                    elseif j == 11 || j == 13
                        % Inner contour for the slice
                        bw_inner = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                        
                        % Subtract inner from outer to create concentric mask
                        bw_concentric = bw_outer & ~bw_inner;
                        mask3d(:, :, round(imgIdx(3) + 1)) = bw_concentric;
                    end
                elseif strcmp(getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName, 'Cervix + HRCTV') && (j > 5) && (j < 8)
                    if j == 6
                        bw_outer = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                    end
                    if j == 7
                        bw_inner = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                        bw_combo = bw_outer & ~bw_inner;
                        mask3d(:, :, round(imgIdx(3) + 1)) = bw_combo;
                    end
                    % if j < 8
                    %     bw = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                    %     mask3d_cerv(:, :, round(imgIdx(3) + 1)) = bw;
                    % else % j >= 8
                    %     bw = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                    %     mask3d_hrctv(:, :, round(imgIdx(3) + 1)) = bw;
                    % end
                else
                    % Normal contour processing
                    bw = poly2mask(in_points(1, :), in_points(2, :), image_height, image_width);
                    mask3d(:, :, round(imgIdx(3) + 1)) = bw;
                    
                    if strcmp(getfield(rtStructInfo.StructureSetROISequence, roiIndex).ROIName, 'Cervix + HRCTV')
                        figure;
                        hold on;
                            subplot(121)
                                hold on;
                                imshow(mask3d(:,:,round(imgIdx(3) + 1)));
                                plot(in_points(1,:), in_points(2,:), 'b', 'LineWidth',2)
                                hold off
                            subplot(122)
                                hold on;
                                imshow(dicomread(tempInfo), [])
                                plot(in_points(1,:), in_points(2,:), 'b', 'LineWidth',2)
                                hold off;
                    end 
                end
            end
        end
    end
    
    switch nargout
        case 1
            varargout{1} = mask3d;
        case 2
            varargout{1} = mask3d_cerv;
            varargout{2} = mask3d_hrctv;
    end
end