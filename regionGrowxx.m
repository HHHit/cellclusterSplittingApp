function outputID = regionGrowxx(inputID,iters, lx,ly,xxshift, yyshift)
   %%% input the numbers of the circles to grow
    [idX, idY] = ind2sub([lx,ly],inputID);
    for i = 1:iters
        idX = min(max(idX + xxshift(:)',1), lx);
        idY = min(max(idY + yyshift(:)',1),ly);
        outputID = sub2ind([lx,ly],idX(:),idY(:));
        outputID = unique(outputID);
        [idX, idY] = ind2sub([lx, ly], outputID);
    end
    
end