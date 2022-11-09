function [xxximg, G, maxID, allids3] = minCutkernel(mask2roi,foreground,score2)
%allids = label2idx(mask2roi);
mask2roix  = mask2roi;
mask2roix(foreground(:) > 0 & mask2roi(:) == 0) = -1;
mask2roix = mask2roix + 1;
allids2 = label2idx(mask2roix);
maxID = length(allids2);%%% This refers to all the nodes that could be linked to source or sink 
%%%(one background and all the other seeds)
allids3 = cell(maxID + sum(mask2roix(:) == 0),1);
allids3(1:maxID) = allids2;
idCur0 = find(mask2roix(:) == 0);
allids3(maxID+1:end) = num2cell(idCur0);
mask2roix(idCur0) = (maxID+1):length(allids3);

%If considering 8 neighbor
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift(i+2,j+2) = i;
        yyshift(i+2,j+2) = j;
    end
end
%if considering 4 neighbors
% xxshift = [1,0,-1,0];
% yyshift = [0,1,0,-1];
[headNodes, tailNodes, weight] = adj8WtListInvScoreFast(mask2roix,maxID+1,allids3,score2, xxshift, yyshift);
G = graph(headNodes, tailNodes, weight);
maxGnode = size(G.Nodes,1);

xxx = cell(max(mask2roi(:)),1);

for i = 1:max(mask2roi(:))
%     disp(i);
    source = i+1;
    sink = [1:maxID];
    sink(i+1) = [];
    G2 = addnode(G, 1); %add two nodes which are the source and sink nodes
    G2 = addedge(G2,sink,repmat(maxGnode+1,length(sink),1),inf(length(sink),1));
    [~, ~,cs1,~] = maxflow(G2,source,maxGnode+1);
    cs1(cs1 == source) = [];
    xxx{i} = cs1;
end


xxximg = zeros(size(mask2roi));
for kk = 1: length(xxx)
    labelx = xxx{kk};
    idx = cell2mat(allids3(labelx));
    idx = [idx;allids3{kk+1}];
    xxximg(idx) = kk;
end
end


function [SS, TT, weight] = adj8WtListInvScoreFast(foreground2, startingLabel,labelIdx, score2, xxshift, yyshift)
%The edge only exists between 
    SS = [];
    TT = [];
    weight = [];
    [lx, ly] = size(foreground2);
    maxNorm = 1/min(score2(score2(:) > 0));
    
    for i =startingLabel:max(foreground2(:))
        curIdx = labelIdx{i};
        curScore = score2(curIdx);
        [idx, idy] = ind2sub([lx,ly],curIdx);
        idxN = idx + xxshift;
        idyN = idy + yyshift;
        idxN(idxN(:)< 1) = 1;
        idxN(idxN(:) > lx) = lx;
        idyN(idyN(:) < 1) =1;
        idyN(idyN(:) > ly) = ly;
        Nid = sub2ind([lx,ly], idxN, idyN);
        Nid = unique(Nid);
        Nid(Nid(:) == curIdx) = [];
        fOrb = foreground2(Nid);
        unifOrb = fOrb;
        %If one of the node comes from the background or from the other
        %seeds, we assign the weight to both direction. Otherwise, we only
        %assign the weight in one direction, because the other direction
        %will be given by another pixel (node)
        for j = 1:length(unifOrb)
            scoreN = score2(Nid(fOrb == unifOrb(j)));
            if(unifOrb(j) < startingLabel)
                scoreN = mean(scoreN);
%                 mlpx = (scoreN + curScore)/2;
                mlpx = scoreN * curScore;
%                 mlpx = 1./mlpx;
%                 mlpx = 1./mlpx;
                mlpx = 1./(mlpx).^4;
                SS = [SS,i,unifOrb(j)];
                TT = [TT,unifOrb(j),i];
                weight = [weight,mlpx,mlpx];
            else
%                 mlpx = (scoreN + curScore)/2;
                mlpx = scoreN * curScore;
%                 mlpx = 1./mlpx;
                mlpx = 1./(mlpx).^4;
                SS = [SS i];
                TT = [TT unifOrb(j)];
                weight = [weight,mlpx];
                
            end
        end
    end
    weight = weight./maxNorm;
    
end