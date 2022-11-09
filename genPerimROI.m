function perimROIzz = genPerimROI(growedROICell, boundingBoxCell, lenx, leny)
perimROIzz = zeros(lenx, leny);
for i =  1:length(growedROICell)
    boundingBox = boundingBoxCell{i};
    growedROI = growedROICell{i};
    if(max(growedROI(:)) == 1)
        perimROItmp = double(bwperim(growedROI));
        perimROIzz(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2)) = perimROIzz(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2)) + perimROItmp;
    else
        perimROItmp = zeros(size(growedROI));
        for j = 1:max(growedROI(:))
            perimROItmp = perimROItmp + bwperim(double(growedROI == j));
        end
        perimROIzz(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2)) = perimROIzz(boundingBox(1,1):boundingBox(1,2), boundingBox(2,1):boundingBox(2,2)) + perimROItmp;

    end

end
end