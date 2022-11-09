function filteredGap = gapRegionMinPerimeter(gapID, objectID, lx, ly)
    SS = [];
    TT = [];
    weight = [];
    %now we only consider the 4-connectivity
    %8 connectivity will be biased towards the backgrounds and each seed
%     xxshift = [1,0,-1,0];
%     yyshift = [0,1,0,-1];
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
%     idall = repelem(gapID, 4);
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
    %     G2 = addedge(G,1,2,inf);
        [~, ~,cs1,~] = maxflow(G,source,sink);
        cs1(cs1>maxlabelGap)=[];
        filteredGap = gapID(cs1);
    end

end
