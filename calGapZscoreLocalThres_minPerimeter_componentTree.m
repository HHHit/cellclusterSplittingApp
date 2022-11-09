function [gapZscoreOneZlice, prGAPValuezslice, seedGraphMx,seedMaskROI,gapZscorePair] = ...
calGapZscoreLocalThres_minPerimeter_componentTree(imCurz, mask_ss, pr_ss,lenx,leny,gapZscoreOneZlice,prGAPValuezslice,xxshift, yyshift,svarx)

pr_ss2 = pr_ss;
pr_ss2(pr_ss2 <0) = 0;
pr_ss2 = (max(pr_ss2(:)) - pr_ss2).*mask_ss;
pr_ss2 = pr_ss2./max(pr_ss2(:)).*32000;
seedMask = double(pr_ss <= 0).*mask_ss;
seedMaskROI = bwlabel(seedMask);
seedMaskROIidx = label2idx(seedMaskROI);
seedMaskROIidx(cellfun(@length, seedMaskROIidx) < 2) = [];



seedMaskROI = zeros(lenx, leny);
for j = 1:length(seedMaskROIidx)
    seedMaskROI(seedMaskROIidx{j}) = j;

end
seedMask = double(seedMaskROI > 0);
gapZscorePair = [];
seedGraphMx = zeros(length(seedMaskROIidx));
if(length(seedMaskROIidx) > 1)
    [seedNeighborInfo] = searchThres4PartitionComponentTree(seedMask, seedMaskROIidx,pr_ss2, mask_ss);

    if(~isempty(seedNeighborInfo))
    sLabel = seedNeighborInfo(:,1);
    tLabel = seedNeighborInfo(:,2);
    maxThres = seedNeighborInfo(:,3);
    [~, sortedID] = sort(maxThres ,'descend');
    sLabel  =sLabel(sortedID);
    tLabel = tLabel(sortedID);
    sLabel = sLabel(:);
    tLabel = tLabel(:);
    gapZscorePair = zeros(length(sLabel),3);
    gapZscorePair(:,1) = sLabel(:);
    gapZscorePair(:,2) = tLabel(:);
    gapZscorePair(:,3) = nan;

    maxThres = maxThres(sortedID) + 1;
    thresRecord = cell(length(maxThres),1);
    for i = 1:length(maxThres)
        if(maxThres(i) == 32000)
            thresRecord{i} = 31999;
        else
            thresRecord{i} = [31999:-(31999 - maxThres(i))/9:(maxThres(i))];
        end
    end
    iniGraph = [sLabel(:), tLabel(:), maxThres(:)];
    if(~isempty(sLabel))
        seedGraphMx = zeros(length(seedMaskROIidx));
    for i = 1:length(sLabel)
        seedGraphMx(sLabel(i), tLabel(i)) = 1; 
    end
    if(length(seedMaskROIidx) > 1)
        edgeWeight = sort(unique(maxThres), 'descend');
        for i = 1:length(edgeWeight)
            if(i > 1)
                idSelected = find(iniGraph(:,3) >= edgeWeight(i-1));
                g_tmp = graph([sLabel(idSelected);[1:length(seedMaskROIidx)]'], [tLabel(idSelected);[1:length(seedMaskROIidx)]']);
                bins = conncomp(g_tmp);
                bins = bins(:);
                extGraph0 = [bins(sLabel), bins(tLabel)]; 
                %when the new pair of id are the same need to make sure the
                %original pairs with the same pair of new ID are given the same
                %z-score update
                iniGraph2 = [iniGraph, extGraph0];
                extGraph0 = sort(extGraph0,2);
                ibrm = find(extGraph0(:,1) == extGraph0(:,2));
                iniGraph2(ibrm,:) = [];
                thresRecord2 = thresRecord;
                thresRecord2(ibrm) = [];
                extGraph1 = extGraph0;
                extGraph1(ibrm, :) = [];
                [extGraph2, ia, ic] = unique(extGraph1,'rows');
                iniGraph3 = [iniGraph2(ia,1:3), extGraph2];           
                thresRecord2 = thresRecord2(ia);   
            else
                iniGraph3 = iniGraph;
                thresRecord2 = thresRecord;
                bins = [];
            end
            for j = 1:size(iniGraph3,1)
                label1ws = iniGraph3(j,1);
                label2ws = iniGraph3(j,2);
                idLabel1seed = seedMaskROIidx{label1ws};
                idLabel2seed = seedMaskROIidx{label2ws};
                thresAll = thresRecord2{j};
                thresAll(thresAll < edgeWeight(i)) = [];
                if(i>1)
                    thresAll(thresAll >= edgeWeight(i-1)) = [];
                end
                for k = 1:length(thresAll)
                [gapfinid,fgid] = searchGapfinid_min_perimeter_componentTree(mask_ss, pr_ss2,seedMask,thresAll(k),lenx,leny, xxshift,yyshift,idLabel1seed,idLabel2seed);
                gapfinid(mask_ss(gapfinid) == 0)=[];
                fgid(mask_ss(fgid) == 0) = [];
                    if(length(gapfinid) >= 5 && length(fgid)>=5)
                        signal = imCurz(fgid);
                        signalNei = imCurz(gapfinid);
                        [mutmp, sigmatmp] = ksegments_orderstatistics_v2(signal, signalNei);
                        meanDiff = mean(signal) - mean(signalNei);
                        zscoretmp = (meanDiff - mutmp.*svarx)/(svarx.*sigmatmp);
                        if(isempty(bins))
                            curPairZscoretmp = gapZscorePair((gapZscorePair(:,1) == label1ws & gapZscorePair(:,2) == label2ws), 3);
                        else
                            refPair = [bins(label1ws), bins(label2ws)];
                            refPair = sort(refPair);
                            ccPairZscoretmp = gapZscorePair(ismember(extGraph0, refPair,'rows'), 3);
                            ccPairZscoretmp(isnan(ccPairZscoretmp)) = 0;
                            curPairZscoretmp =  min(ccPairZscoretmp);
                        end

                        if((isnan(curPairZscoretmp) && (zscoretmp > 0))||((zscoretmp > curPairZscoretmp) && (zscoretmp > 0)))
                            gapZscoreOneZlice(gapfinid) = max(gapZscoreOneZlice(gapfinid),zscoretmp);
                            prGAPValuezslice(gapfinid(gapZscoreOneZlice(gapfinid) == zscoretmp)) = thresAll(k);
                            if(isempty(bins))
                                gapZscorePair((gapZscorePair(:,1) == label1ws & gapZscorePair(:,2) == label2ws), 3) = zscoretmp;
                            else
                                refPair = [bins(label1ws), bins(label2ws)];
                                refPair = sort(refPair);
                                gapZscorePair(ismember(extGraph0, refPair,'rows'), 3) = max(gapZscorePair(ismember(extGraph0, refPair,'rows'), 3), zscoretmp);
                            end
                        end
                    end 
                end
            end
        end
    end
    end
    end    
end
    

end

