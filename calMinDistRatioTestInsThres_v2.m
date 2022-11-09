function [pGeometryMap_ss,thresIns_ss, gapPvaluePair] = calMinDistRatioTestInsThres_v2(mask_ss, im_Smoothed_ss, lenx,leny, distTransform, pGeometryMap_ss,muRatio,varRatio,seedGraphMx,seedMaskROI)
    % calculate min Ratio based on the thresholding of the intensity 
    % Also output the matrix of the seeds relationship
    [sLabel, tLabel] = find(seedGraphMx == 1);
    seedGraphMxSym = seedGraphMx + seedGraphMx';
    gapPvaluePair = zeros(length(sLabel),3);
    gapPvaluePair(:,1) = sLabel(:);
    gapPvaluePair(:,2) = tLabel(:);
    gapPvaluePair(:,3) = nan;
    im_ss = im_Smoothed_ss;
    thresIns_ss = zeros(lenx, leny);
    maxIns = max(im_ss(:));
    minIns = min(im_ss(mask_ss(:)  == 1));
    thresAll = [minIns:(maxIns - minIns)/19:maxIns];
    allMinDistRatio = [];
    skelx = bwskel(logical(mask_ss));
    skelxid = find(skelx(:) == 1);
    [skelxidx, skelxidy] = ind2sub([lenx, leny], skelxid);
    selectedRegionBound = zeros(lenx, leny);
    selectedRegionBound(2:(lenx-1), 2:(leny - 1)) = 1;
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


                    bwdistTwo = bwdistgeodesic(logical(mask_ss), [idLabel1Now(:)],'quasi-euclidean') + bwdistgeodesic(logical(mask_ss), [idLabel2Now(:)],'quasi-euclidean');
                    bwdistTwo(isnan(bwdistTwo)) = inf;
                    paths = double(bwdistTwo <= (min(bwdistTwo(:)) + 0.001));
                    paths_thinned_many = bwmorph(paths, 'thin', inf);
                    thinLineid = find(paths_thinned_many(:) == 1);
                    thinLineid = setdiff(thinLineid, [idLabel1Now;idLabel2Now]);
                    thinLineLabel = unique(binMaskROI(thinLineid));
                    if(min(thinLineLabel == 0))
                        [obj1IDX, obj1IDY] = ind2sub([lenx,leny],idLabel1Now);
                        [obj2IDX, obj2IDY] = ind2sub([lenx,leny],idLabel2Now);
                        [thinLineidx, thinLineidy] = ind2sub([lenx, leny], thinLineid);
                        distanceAll = sqrt((thinLineidx(:) - skelxidx(:)').^2 + (thinLineidy(:) - skelxidy(:)').^2);
                        [~,skelSelectedID] = min(distanceAll,[],2);  
                        radSeed = estimateRadius(idLabel1Now, idLabel2Now, lenx, leny);
                        minRatio = radSeed/min(distTransform(skelxid(skelSelectedID)));
                        %
                        minDist = min(min(floor(sqrt((obj1IDX - obj2IDX').^2 + (obj1IDY - obj2IDY').^2))));
                        centerDist = floor(sqrt((mean(obj1IDX) - mean(obj2IDX)).^2 + (mean(obj1IDY) - mean(obj2IDY)).^2));
                        growIter = max(floor((centerDist + minDist)/2),3);
                        idLabel1Grow = regionGrowMatrix(idLabel1Now, growIter,lenx, leny);
                        idLabel2Grow = regionGrowMatrix(idLabel2Now, growIter,lenx, leny);
                        thinLineid = find(paths_thinned_many(:) == 1);


                        idSeed = [idLabel1Now(:);idLabel2Now(:);thinLineid];
                        overlappedid = intersect(idLabel1Grow, idLabel2Grow); 
                        overlappedid(binMask(overlappedid) > 0) = []; % largest possible gap region 
                        %             overlappedid(pr_ss(overlappedid) < maxThresLocal) = [];
                        overlappedid(paths_thinned_many(overlappedid) == 1) = [];
                        overlappedid(mask_ss(overlappedid) == 0) = [];
                        overlappedid(selectedRegionBound(overlappedid) == 0) = [];
                        gapfinid = gapRegionMinPerimeter(overlappedid, idSeed, lenx, leny);
                        thinLineid = setdiff(thinLineid, [idLabel1Now;idLabel2Now]);
                        gapfinid = [gapfinid;thinLineid];
                        maskWith2Seed = zeros(lenx, leny);
                        maskWith2Seed(idLabel1Now) = 1;
                        maskWith2Seed(idLabel2Now) = 2;

                        maskWithGap = zeros(lenx, leny);
                        maskWithGap(gapfinid) = 1;
                        maskWithGapROI = bwlabel(maskWithGap,8);
                        maskWithGapROIidx = label2idx(maskWithGapROI);
                        se = strel('disk',1);
                        maskWith2SeedROI = imdilate(maskWith2Seed,se);

                        if(length(maskWithGapROIidx) > 1)
                            
                            labelNum = cellfun(@(c) sum(unique(maskWith2SeedROI(c)) > 0), maskWithGapROIidx);
                            maskWithGapROIidx(labelNum<2) = [];
                            gapfinid = cell2mat(maskWithGapROIidx(:));
                        end
                        
                        
                        pvalueTmp = normcdf(muRatio - log(minRatio),0, varRatio);
                        if(pvalueTmp < min(pGeometryMap_ss(gapfinid)) )
                            
                            pGeometryMap_ss(gapfinid) = pvalueTmp;
                            thresIns_ss(gapfinid) = thresAll(i);
                            seedLabel1All = unique(seedMaskROI(idLabel1Now));
                            seedLabel1All(seedLabel1All == 0) = [];
                            seedLabel1All = seedLabel1All(:);
                            seedLabel2All = unique(seedMaskROI(idLabel2Now));
                            seedLabel2All(seedLabel2All == 0) = [];
                            seedLabel2All = seedLabel2All(:);
                            [label1PairId, label2PairId] = find(seedGraphMxSym(seedLabel1All,seedLabel2All) == 1);
                            possiblePairs = [seedLabel1All(label1PairId(:)), seedLabel2All(label2PairId(:))];
                            possiblePairs = sort(possiblePairs,2);
                            if(~isempty(possiblePairs))
                                gapPvaluePair(ismember(gapPvaluePair(:,1:2), possiblePairs,'rows'),3) = pvalueTmp;
                            end
                        end
%                         allMinDistRatio = [allMinDistRatio;minRatio];
                        
                        
                        
                    end
                end
            end
                
                
        
        end
    end
end
