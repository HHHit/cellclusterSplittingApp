function [seedNeighborInfo] = searchThres4PartitionComponentTree(seedMask, seedMaskROIidx,scoreMap, mask_ss)
% the first column of the output is sLabel the second is the tLabel the
% third is the max threshold that the two are met
[lenx, leny] = size(seedMask);
score_combined_int = int16(scoreMap);
%all the threshold pair between seeds
% score_combined_int = int16(scoreMap./max(scoreMap(:)).*32000);
tree = maxtree_of_image(score_combined_int,8);
nodeIDmap = zeros(lenx, leny);
for i = 1:length(tree)
   nodeIDmap(double(tree(i).Pixels) + 1) = i; 
end
startingNode = cellfun(@(c) max(nodeIDmap(c)), seedMaskROIidx);
parents_Nodes = cell(length(startingNode),1);
for i = 1:length(startingNode)
    parent_cur = [];
    startNode_tmp = startingNode(i);
    parent_tmp = startNode_tmp;
   while(~isempty(parent_tmp))
       parent_tmp = tree(startNode_tmp).Parent;
       parent_cur = [parent_cur;parent_tmp];
       startNode_tmp = parent_tmp;
   end
    parents_Nodes{i} = double(parent_cur);
end

thresMerged = zeros(length(startingNode), length(startingNode)) + nan;

for m = 1:length(startingNode)-1
    for n = (m+1):length(startingNode)
        nodeID = max(intersect(parents_Nodes{m}, parents_Nodes{n}));
        thresMerged(m,n) = tree(nodeID).GrayLevel;
    end
end
thresPair = thresMerged;


% watershed to know the neighbor info
scoreMapInv = (1./scoreMap);
seedMaskROI = bwlabel(seedMask);

scoreMapInv = imimposemin(scoreMapInv, seedMask);
scoreMapInvWD = watershed(scoreMapInv);
scoreMapInvWD = double(scoreMapInvWD).*mask_ss;
scoreMapInvWDidx = label2idx(scoreMapInvWD);
scoreMapInvWDx = zeros(lenx, leny);
scoreMapInvWDGrowedidx = cell(length(scoreMapInvWDidx),1);
for kk = 1:length(scoreMapInvWDidx)
    newLabel = max(seedMaskROI(scoreMapInvWDidx{kk}));
    scoreMapInvWDx(scoreMapInvWDidx{kk}) = newLabel;
    scoreMapInvWDGrowedidx{newLabel} = regionGrowMatrix(scoreMapInvWDidx{kk}, 1, lenx, leny);
end

seedNeighborInfo = [];
for k1 = 1:(length(scoreMapInvWDGrowedidx) - 1)
   for k2= (k1+1):length(scoreMapInvWDGrowedidx) 
       intersecID = intersect(scoreMapInvWDGrowedidx{k1}, scoreMapInvWDGrowedidx{k2});
       intersecID(scoreMapInvWDx(intersecID) > 0) = [];
       if(~isempty(intersecID))
           seedNeighborInfo = [seedNeighborInfo;[k1, k2,thresPair(k1, k2)]];
       end
   end
end


end