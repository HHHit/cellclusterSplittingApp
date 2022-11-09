function allMinDistRatio = calMinDistRatioGatherStatInsThres_v1(mask_ss, im_Smoothed_ss, lenx,leny, distTransform)
    % calculate min Ratio based on the thresholding of the intensity 
    im_ss = im_Smoothed_ss;
    maxIns = max(im_ss(:));
    minIns = min(im_ss(mask_ss(:)  == 1));
    thresAll = [maxIns:-(maxIns - minIns)/19:minIns];
    allMinDistRatio = [];
    skelx = bwskel(logical(mask_ss));
    skelxid = find(skelx(:) == 1);
    [skelxidx, skelxidy] = ind2sub([lenx, leny], skelxid);
    selectedRegionBound = zeros(lenx, leny);
%     selectedRegionBound(2:(lenx-1), 2:(leny - 1)) = 1;
    for i = 1:length(thresAll)
        % check if the shortest path goes through only the two seed
        binMask = double(im_ss >= thresAll(i));
        binMaskROI = bwlabel(binMask);
        binMaskROIidx = label2idx(binMaskROI);
        if(length(binMaskROIidx) > 1)
            for m =1:(length(binMaskROIidx)-1)
                for n = 2:length(binMaskROIidx)
                    idLabel1Now = binMaskROIidx{m};
                    idLabel2Now = binMaskROIidx{n};
                    bwdistTwo = bwdistgeodesic(logical(mask_ss), [idLabel1Now(:)],'cityblock') + bwdistgeodesic(logical(mask_ss), [idLabel2Now(:)],'cityblock');
                    bwdistTwo(isnan(bwdistTwo)) = inf;
                    paths = double(bwdistTwo <= (min(bwdistTwo(:)) + 0.001));
                    paths_thinned_many = bwmorph(paths, 'thin', inf);
                    thinLineid = find(paths_thinned_many(:) == 1);
                    thinLineid = setdiff(thinLineid, [idLabel1Now;idLabel2Now]);
                    thinLineLabel = unique(binMaskROI(thinLineid));
                    if(thinLineLabel == 0)
%                         [obj1IDX, obj1IDY] = ind2sub([lenx,leny],idLabel1Now);
%                         [obj2IDX, obj2IDY] = ind2sub([lenx,leny],idLabel2Now);
                        [thinLineidx, thinLineidy] = ind2sub([lenx, leny], thinLineid);
                        distanceAll = sqrt((thinLineidx(:) - skelxidx(:)').^2 + (thinLineidy(:) - skelxidy(:)').^2);
                        [~,skelSelectedID] = min(distanceAll,[],2);  
                        
                        radSeed = estimateRadius(idLabel1Now, idLabel2Now, lenx, leny);
                        minRatio = radSeed/min(distTransform(skelxid(skelSelectedID)));
                        allMinDistRatio = [allMinDistRatio;[minRatio,radSeed,min(distTransform(skelxid(skelSelectedID)))]];
                        
                        
                        
                    end
                end
            end
                
                
        
        end
    end
end
