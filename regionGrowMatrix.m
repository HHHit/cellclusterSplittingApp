function outputID = regionGrowMatrix(inputID,iters, lx,ly)
   %%% input the numbers of the circles to grow
    [idX, idY] = ind2sub([lx,ly],inputID);
    xxshift = zeros(2*iters+1,2*iters+1);
    yyshift = zeros(2*iters+1,2*iters+1);
    for i = -iters:iters
        for j = -iters:iters
            xxshift(i+iters + 1,j+iters + 1) = i;
            yyshift(i+iters + 1,j+iters + 1) = j;
        end
    end
    
    idX = min(max(idX + xxshift(:)',1), lx);
    idY = min(max(idY + yyshift(:)',1),ly);
    outputID = sub2ind([lx,ly],idX(:),idY(:));
%     outputID = unique(outputID(:));

    
end