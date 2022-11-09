function [gapfinid,fgid] = searchGapfinid_min_perimeter_componentTree(mask_ss, pr_ss2, seed_ss, maxThresLocal,lenx,leny, xxshift,yyshift,idLabel1seed,idLabel2seed)
% 8 connectivity
        se = strel('disk',1);
        mask_ss_erosion = imerode(mask_ss, se); % grow one round to ensure the shortest path can be connected to the overlapped region
        tmpMask = mask_ss.*(pr_ss2 > maxThresLocal);
        tmpMaskROi = bwlabel(tmpMask,8);
        tmpMaskROiidx = label2idx(tmpMaskROi);
        isWithSeed = cellfun(@(c) max(seed_ss(c)), tmpMaskROiidx);
        tmpMaskROiidx = tmpMaskROiidx(isWithSeed > 0);
        tmpMaskROi = zeros(lenx,leny);
        for i = 1:length(tmpMaskROiidx)
            tmpMaskROi(tmpMaskROiidx{i}) = i;
        end
        Label1Now = unique(tmpMaskROi(idLabel1seed));
        Label1Now(Label1Now == 0) = [];
        Label2Now = unique(tmpMaskROi(idLabel2seed));
        Label2Now(Label2Now == 0) = [];
        if(Label1Now == Label2Now)
           warning('not well separated!'); 
           gapfinid = [];
           fgid = [];
           disp(maxThresLocal);
        else
            tmpMaskROiidx = label2idx(tmpMaskROi);
            idLabel1Now = tmpMaskROiidx{Label1Now};
            idLabel2Now = tmpMaskROiidx{Label2Now};
            [obj1IDX, obj1IDY] = ind2sub([lenx,leny],idLabel1Now);
            [obj2IDX, obj2IDY] = ind2sub([lenx,leny],idLabel2Now);
            bwdistTwo = bwdistgeodesic(logical(mask_ss), [idLabel1Now(:)],'quasi-euclidean') + bwdistgeodesic(logical(mask_ss), [idLabel2Now(:)],'quasi-euclidean');
            bwdistTwo(isnan(bwdistTwo)) = inf;
            paths = imregionalmin(bwdistTwo);
            paths_thinned_many = bwmorph(paths, 'thin', inf);
            minDist = min(min(floor(sqrt((obj1IDX - obj2IDX').^2 + (obj1IDY - obj2IDY').^2))));
            centerDist = floor(sqrt((mean(obj1IDX) - mean(obj2IDX)).^2 + (mean(obj1IDY) - mean(obj2IDY)).^2));
            growIter = max(floor((centerDist + minDist)/2),3);
            idLabel1Grow = regionGrowMatrix(idLabel1Now, growIter,lenx, leny);
            idLabel2Grow = regionGrowMatrix(idLabel2Now, growIter,lenx, leny);
            thinLineid = find(paths_thinned_many(:) == 1);
            idSeed = [idLabel1Now(:);idLabel2Now(:);thinLineid];
            overlappedid = intersect(idLabel1Grow, idLabel2Grow); 
            overlappedid(tmpMaskROi(overlappedid) > 0) = []; % largest possible gap region 
            overlappedid(paths_thinned_many(overlappedid) == 1) = [];
            overlappedid(mask_ss(overlappedid) == 0) = [];
            gapfinid = gapRegionMinPerimeter(overlappedid, idSeed, lenx, leny);
            thinLineid = setdiff(thinLineid, [idLabel1Now;idLabel2Now]);
            gapfinid = [gapfinid;thinLineid];
            fgid = [];
            iters = 0;
            %filter any gap regions that are not in between two seeds
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
            gapfinid(mask_ss(gapfinid) == 0) = [];
            gapfinidGrow = gapfinid;
            if(~isempty(gapfinid))    
                while((length(fgid) < 5) && (iters <10))  
                    gapfinidGrow = regionGrowxx(gapfinidGrow,1, lenx,leny,xxshift, yyshift);
                    fgid = intersect(gapfinidGrow, unique([idLabel1Now;idLabel2Now]));
                    iters = iters + 1;
                end  
            else
                gapfinid = [];
                fgid = [];
            end
            
            
        end
 
       
end
function filteredGap = gapRegionMinPerimeter(gapID, objectID, lx, ly)
    SS = [];
    TT = [];
    weight = [];
    %consider 8-connectivity
    xxshift = zeros(3,3);
    yyshift = zeros(3,3);
    for i = -1:1
        for j = -1:1
            xxshift(i+2,j+2) = i;
            yyshift(i+2,j+2) = j;
        end
    end
    nodeMap = zeros(lx, ly);
    nodeMap(gapID) = 1:length(gapID);
    maxlabelGap = length(gapID);
    [idx, idy] = ind2sub([lx,ly],gapID); 
    idxN = idx(:) + xxshift(:)';
    idyN = idy(:) + yyshift(:)';
    idxN(idxN(:)< 1) = 1;
    idxN(idxN(:) > lx) = lx;
    idyN(idyN(:) < 1) =1;
    idyN(idyN(:) > ly) = ly;
    idxN = idxN';
    idyN = idyN';
    idNall = sub2ind([lx,ly], idxN(:), idyN(:));
    idall = repelem(gapID, 9);
    objectID2 = intersect(objectID, idNall);
    
    nodeMap(objectID2) = (maxlabelGap+1):(maxlabelGap+length(objectID2));
    maxlabelObject = maxlabelGap + length(objectID2);
    bgID = setdiff(idNall, [gapID;objectID]);
    nodeMap(bgID) = (maxlabelObject+1):(maxlabelObject+length(bgID));
    maxlabelbg = maxlabelObject + length(bgID);
    nodelist = [nodeMap(idall(:)), nodeMap(idNall(:))];
    nodelist = sort(nodelist,2);
    nodelist = unique(nodelist, 'rows');
    nodelist = [nodelist, ones(size(nodelist,1),1)];
    nodelist(nodelist(:,1) == nodelist(:,2),:) = [];
    source = maxlabelbg + 1;
    nodelistsource = [repmat(source, length(objectID2),1), [(maxlabelGap+1):(maxlabelGap+length(objectID2))]',inf( length(objectID2),1)];
    sink = maxlabelbg + 2;
    nodelistsink = [repmat(sink, length(bgID),1), [(maxlabelObject+1):(maxlabelObject+length(bgID))]',inf(length(bgID),1)];
    nodelist = [nodelist;nodelistsource;nodelistsink];
    if(isempty(bgID))
        filteredGap = gapID;
        
    else
        G = graph(nodelist(:,1), nodelist(:,2), nodelist(:,3));
        [~, ~,cs1,~] = maxflow(G,source,sink);
        cs1(cs1>maxlabelGap)=[];
        filteredGap = gapID(cs1);
    end

end
