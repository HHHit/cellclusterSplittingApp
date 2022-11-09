

function [growedROICell] = regionSplit_withThres(mask_zz_ROI_idx,imCurzG,combinedPvalue,seedMaskAll,boundingBoxCell,lenx, leny, curThres)

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
        linkGraph = gapPvaluePair(gapPvaluePair(:,3) > curThres,1:2); %not significant case
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
end