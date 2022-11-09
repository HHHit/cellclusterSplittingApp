function [muRatioThinNeck, sigmaRatioThinNeck,boundingBoxCell] = buildEmpDstrThinNeck(maskAll,imSmoothed, zz, lenx, leny)

mask_zz = maskAll(:,:,zz);
imCurzG = imSmoothed(:,:,zz);
imCurzG = imCurzG.*mask_zz;
mask_zz_ROI = bwlabel(mask_zz);
mask_zz_ROI_idx = label2idx(mask_zz_ROI);
mask_zz_ROI_idx(cellfun(@length, mask_zz_ROI_idx) < 30)= [];
mask_zz_ROI = zeros(lenx, leny);
for j = 1:length(mask_zz_ROI_idx)
    mask_zz_ROI(mask_zz_ROI_idx{j}) = j;
end

distZZ = bwdist(1 - mask_zz);
ratioAllzzCell = cell(length(mask_zz_ROI_idx),1);
boundingBoxCell = cell(length(mask_zz_ROI_idx),1);
parfor ii = 1:length(mask_zz_ROI_idx)
mask_ss = zeros(lenx, leny);
mask_ss(mask_zz_ROI_idx{ii}) = 1;
dist_ss = distZZ.*mask_ss;
imCurzG_ss = imCurzG.*mask_ss;
ssid = mask_zz_ROI_idx{ii};

[ssidx, ssidy] = ind2sub([lenx, leny], ssid);

boundingBox = [max(min(ssidx)-1,1), min(max(ssidx)+1, lenx); max(min(ssidy)-1,1),min(max(ssidy)+1, leny)];
boundingBoxCell{ii} = boundingBox;
mask_ss = mask_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
dist_ss = dist_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
imCurzG_ss = imCurzG_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
[lenxtmp,lenytmp] = size(mask_ss);
minDistRatio_ss = calMinDistRatioGatherStatInsThres_v1(mask_ss, imCurzG_ss,lenxtmp,lenytmp, dist_ss); 
minDistRatio_ss = unique(minDistRatio_ss,'rows');
ratioAllzzCell{ii} = double(minDistRatio_ss);
end
ratioAllzz = cell2mat(ratioAllzzCell);
ratioAllzzLog= log(ratioAllzz(:,1));
xx = ratioAllzzLog(ratioAllzz(:,3) >sqrt(2));
[muRatioThinNeck] = double(mean(xx));
sigmaRatioThinNeck = sqrt(double(fitTruncatedGaussian(xx, 0.05)));
end
