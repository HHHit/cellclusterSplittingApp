function [seedMaskAll,combinedPvalue,mask_zz_ROI_idx,boundingBoxCell] = testingTwoFeature_bulk(im,maskAll,prCurvature,imDiffAll,imDiffAll2,imSmoothed, zz, lenx, leny, muRatio, sigmaRatio)
addpath('Maxtree-Processing-Toolbox-master\')
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift(i+2,j+2) = i;
        yyshift(i+2,j+2) = j;
    end
end
imCurz = im(:,:,zz);
mask_zz = maskAll(:,:,zz);
pr_zz = prCurvature(:,:,zz);
imCurzG = imSmoothed(:,:,zz);
imDiff_zz = imDiffAll(:,:,zz);
pr_zz = pr_zz.*mask_zz;

imCurzG = imCurzG.*mask_zz;
mask_zz_ROI = bwlabel(mask_zz);
mask_zz_ROI_idx = label2idx(mask_zz_ROI);
mask_zz_ROI_idx(cellfun(@length, mask_zz_ROI_idx) < 30)= [];
mask_zz_ROI = zeros(lenx, leny);
for j = 1:length(mask_zz_ROI_idx)
    mask_zz_ROI(mask_zz_ROI_idx{j}) = j;
end
% gapZscoreOneZlice = zeros(lenx, leny);
% prGAPValuezslice = zeros(lenx, leny) + inf;
% seedMask = zeros(lenx, leny);
varDiff = cellfun(@(c) var(imDiff_zz(c)), mask_zz_ROI_idx);
imDiff_zz2 = imDiffAll2(:,:,zz);
varDiff2 = cellfun(@(c) var(imDiff_zz2(c)), mask_zz_ROI_idx);
corrFactor = median(sqrt((varDiff./varDiff2)));
distZZ = bwdist(1 - mask_zz);
combinedPvalue = cell(length(mask_zz_ROI_idx),1);
seedMaskAll = cell(length(mask_zz_ROI_idx),1);
boundingBoxCell = cell(length(mask_zz_ROI_idx),1);
parfor ii =1:length(mask_zz_ROI_idx)
% disp(ii)
mask_ss = zeros(lenx, leny);
mask_ss(mask_zz_ROI_idx{ii}) = 1;
imDiff_ss = imDiff_zz2.*mask_ss;
im_ss = imCurz.*mask_ss;
pr_ss = pr_zz.*mask_ss;
dist_ss = distZZ.*mask_ss;
imCurzG_ss = imCurzG.*mask_ss;
ssid = mask_zz_ROI_idx{ii};
[ssidx, ssidy] = ind2sub([lenx, leny], ssid);
boundingBox = [max(min(ssidx)-1,1), min(max(ssidx)+1, lenx); max(min(ssidy)-1,1),min(max(ssidy)+1, leny)];
boundingBoxCell{ii} = boundingBox;
mask_ss = mask_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
imDiff_ss = imDiff_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
im_ss = im_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
pr_ss = pr_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
dist_ss = dist_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
imCurzG_ss = imCurzG_ss(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2));
diffHalfTmp = imDiff_ss(mask_ss(:) > 0);
diffHalfTmp = diffHalfTmp(diffHalfTmp <= 0 );
diffTmp = [diffHalfTmp(:);-diffHalfTmp(:)];
svarx = sqrt(var(diffTmp))*corrFactor;
[lenxtmp,lenytmp] = size(im_ss);
gapZscoreOneZlice_ss = zeros(lenxtmp, lenytmp);
prGAPValuezslice_ss = zeros(lenxtmp, lenytmp) +inf;
[~, ~,seedGraphMx,seedMaskROI, gapZscorepair_darkGap] = calGapZscoreLocalThres_minPerimeter_componentTree(im_ss, mask_ss, pr_ss, lenxtmp,lenytmp, gapZscoreOneZlice_ss,prGAPValuezslice_ss, xxshift, yyshift,svarx);
seedMaskAll{ii} = seedMaskROI;
if(~isempty(gapZscorepair_darkGap))
pGeometryMap_ss = ones(lenxtmp, lenytmp);
[~, ~,gapPvaluePair_geometry] = calMinDistRatioTestInsThres_v2(mask_ss, imCurzG_ss,lenxtmp,lenytmp,dist_ss,pGeometryMap_ss,muRatio,sigmaRatio,seedGraphMx,seedMaskROI); 
gapZscorepair_darkGap = sortrows(gapZscorepair_darkGap);
gapPvaluepair_darkGap_x = normcdf(-gapZscorepair_darkGap(:,3));
gapPvaluepair_darkGap_x(isnan(gapPvaluepair_darkGap_x)) = 1;
gapPvaluePair_geometry = sortrows(gapPvaluePair_geometry);
gapPvaluePair_geometry_x = gapPvaluePair_geometry(:,3);

gapPvaluePair_geometry_x(isnan(gapPvaluePair_geometry_x)) = 1;
gapPvaluePair_combined = -2*(log(gapPvaluepair_darkGap_x) + log(gapPvaluePair_geometry_x));
gapZscorePair_combined = 1 - chi2cdf(gapPvaluePair_combined,4);
combinedPvalue{ii} = [gapPvaluePair_geometry(:,1:2), gapZscorePair_combined(:)];

end
end


end