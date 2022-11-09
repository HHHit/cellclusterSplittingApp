

function [growedROICell,minPvalue, maxPvalue] = regionSplit(mask_zz_ROI_idx,imCurzG,combinedPvalue,seedMaskAll,boundingBoxCell,lenx, leny )

growedROICell = cell(length(mask_zz_ROI_idx),1);
parfor ii = 1:length(mask_zz_ROI_idx)
        mask_ss = zeros(lenx, leny);
        mask_ss(mask_zz_ROI_idx{ii}) = 1;
        imCurzG_ss = imCurzG.*mask_ss;
        boundingBox = boundingBoxCell{ii};
        mask_ss = mask_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
    if(~isempty(combinedPvalue{ii}))
        imCurzG_ss = imCurzG_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
        gapPvaluePair = combinedPvalue{ii};
        linkGraph = gapPvaluePair(gapPvaluePair(:,3) > 0.05,1:2); %not significant case
        seedMaskROI = seedMaskAll{ii};
        seedMaskROIidx = label2idx(seedMaskROI);
        seedMaskROIidx = seedMaskROIidx(:);
        linkGraph = [linkGraph;[[1:length(seedMaskROIidx)]',[1:length(seedMaskROIidx)]']];
        Gx = graph(linkGraph(:,1), linkGraph(:,2));
        bins = conncomp(Gx);
        bins = bins(:);
        [lenxtmp, lenytmp] = size(mask_ss);
        newROI = zeros(lenxtmp, lenytmp);
        for n = 1:max(bins)
            newROI(cell2mat(seedMaskROIidx(bins == n))) = n;
        end
        [growedROI,~,~,~] = minCutkernel(newROI, mask_ss, 1./imCurzG_ss);
    else
        growedROI = mask_ss;
    end
    growedROICell{ii} = growedROI;
end


minPvalue = zeros(length(mask_zz_ROI_idx),1);
maxPvalue = zeros(length(mask_zz_ROI_idx),1);
for ii = 1:length(mask_zz_ROI_idx)
   if(~isempty(combinedPvalue{ii})) 
       cptmp = combinedPvalue{ii};
       pvalueTmp = cptmp(:,3);
       
       if(isempty( max(pvalueTmp(pvalueTmp <0.05))))
           minPvalue(ii) = nan;
       else
           minPvalue(ii) = max(pvalueTmp(pvalueTmp <0.05));
       end
       
       if(isempty( min(pvalueTmp(pvalueTmp > 0.05))))
           maxPvalue(ii) = nan;
       else
           maxPvalue(ii) = min(pvalueTmp(pvalueTmp > 0.05));
       end
   else
       minPvalue(ii) = nan;
       maxPvalue(ii) = nan;
   end    
end
end