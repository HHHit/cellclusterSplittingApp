function radOut = estimateRadius(idLabel1, idLabel2, lenx, leny)
[idLabel1x, idLabel1y] = ind2sub([lenx, leny], idLabel1);
idLabel1x = lenx - idLabel1x;
[idLabel2x, idLabel2y] = ind2sub([lenx, leny], idLabel2);
idLabel2x = lenx - idLabel2x;
center1x = mean(idLabel1x);
center1y = mean(idLabel1y);
center2x = mean(idLabel2x);
center2y = mean(idLabel2y);
if(center1y == center2y)
    p1 = idLabel1y;
    p1 = p1 - mean(p1);
    distance1 = idLabel1x;
    [id1,~] = discretize(distance1,max(floor(sqrt(length(p1))/2),1));
    hist1 = zeros(max(id1),1);
    for k = 1:length(hist1)
        hist1(k) = max(p1(id1 ==k)) - min(p1(id1 == k));
    end
    p2 = idLabel2y;
    p2 = p2 - mean(p2);
    distance2 = idLabel2x;
    [id2,~] = discretize(distance2,max(floor(sqrt(length(p2))/2),1));
    hist2 = zeros(max(id2),1);
    for k = 1:length(hist2)
        hist2(k) = max(p2(id2 ==k)) - min(p2(id2 == k));
    end
    radOut = max((max(hist2) + max(hist1))/4,1);
elseif(center1x == center2x)
    p1 = idLabel1x;
    p1 = p1 - mean(p1);
    distance1 = idLabel1y;
    [id1,~] = discretize(distance1,max(floor(sqrt(length(p1))/2),1));
    hist1 = zeros(max(id1),1);
    for k = 1:length(hist1)
        hist1(k) = max(p1(id1 ==k)) - min(p1(id1 == k));
    end
    p2 = idLabel2x;
    p2 = p2 - mean(p2);
    distance2 = idLabel2y;
    [id2,~] = discretize(distance2,max(floor(sqrt(length(p2))/2),1));
    hist2 = zeros(max(id2),1);
    for k = 1:length(hist2)
        hist2(k) = max(p2(id2 ==k)) - min(p2(id2 == k));
    end
    radOut = max((max(hist2) + max(hist1))/4,1);


else
    a = (center2x - center1x)/(center2y - center1y);
     b = center2x - a*center2y;
    % y = -1/a*x+c
    ccx = (center1x + center2x)/2;
    ccy = (center1y + center2y)/2;
    c = ccx + (1/a)*ccy;
    % projection
    p1y = (a*idLabel1y + c - idLabel1x)/(1/a + a);
    p1x = -1/a.*p1y + c;
    p1abs = sqrt((p1x - ccx).^2 + (p1y - ccy).^2);
    p1abs(p1x < ccx) = -(p1abs(p1x < ccx));
    p1 = p1abs;
    distance1 = sqrt((idLabel1x - p1x).^2 + (idLabel1y - p1y).^2);
    [id1,~] = discretize(distance1,max(floor(sqrt(length(p1))/2),1));
    hist1 = zeros(max(id1),1);
    for k = 1:length(hist1)
        hist1(k) = max(p1(id1 ==k)) - min(p1(id1 == k));
    end
    p2y = (a*idLabel2y + c - idLabel2x)/(1/a + a);
    p2x = -1/a.*p2y + c;
    p2abs = sqrt((p2x - ccx).^2 + (p2y - ccy).^2);
    p2abs(p2x < ccx) = -(p2abs(p2x < ccx));
    p2 = p2abs;
    distance2 = sqrt((idLabel2x - p2x).^2 + (idLabel2y - p2y).^2);
    [id2,~] = discretize(distance2,max(floor(sqrt(length(p2))/2),1));
    hist2 = zeros(max(id2),1);
    for k = 1:length(hist2)
        hist2(k) = max(p2(id2 ==k)) - min(p2(id2 == k));
    end
    radOut = max((max(hist2) + max(hist1))/4,1);
end

end