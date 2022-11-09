function [diffCurvatureScore, u1, v1] = principalCurvatureSeg_2D(im, sigma2D)
% use 2D principal curvature to remove 2D gaps only
[lenx, leny, lenz] = size(im);
imSmoothed = zeros(size(im));

for i = 1:size(im,3)
    
   imSmoothed(:,:,i) = imgaussfilt(im(:,:,i),sigma2D); 
    
end
% use 2D principal curvature to target at the gap at each z stacks
% parfor i = 1:size(diffCurvatureScore,3)
% [lx, ly] = gradient(imSmoothed(:,:,i));
% [lxx,lyx] = gradient(lx);
% [lxy, lyy] = gradient(ly);
% eig2 = zeros(lenx, leny);
% u2 = zeros(lenx, leny);
% v2 = zeros(lenx, leny);
% for ii  =1: lenx*leny
%     MM = [lxx(ii) lxy(ii);lyx(ii) lyy(ii)];
%     [Evec,Eval] = eig(MM);
%     dEval = diag(Eval);
%     [c , id] = sort(dEval);
%     eig2(ii) = dEval(id(2));
%     u2(ii) = Evec(1,id(2));
%     v2(ii) = Evec(2,id(2));
% end
% 
% diffCurvatureScore(:,:,i) = eig2;
% end
diffCurvatureScore = zeros(size(imSmoothed));
u1 = zeros(size(imSmoothed));
v1 = zeros(size(imSmoothed));
for i = 1:size(diffCurvatureScore,3)
[lx, ly] = gradient(imSmoothed(:,:,i));
[lxx,lyx] = gradient(lx);
[lxy, lyy] = gradient(ly);
b = (lxx + lyy);
c = lxx.*lyy - lxy.*lyx;
sumx = b.^2 - 4.*c;
sumx(sumx<0) = 0;
diffCurvatureScore(:,:,i) = (b+sqrt(sumx))/2;
u1tmp= -lxy./(lxx - diffCurvatureScore(:,:,i));
u1tmp(isnan(u1tmp)) = 0;
v1tmp = ones(lenx, leny);
normuv = sqrt(u1tmp.^2 + v1tmp.^2);
u1(:,:,i) = u1tmp./normuv;
v1(:,:,i) = v1tmp./normuv;
% diffCurvatureScore(:,:,i) = (b+sqrt(b.^2 - 4.*c))/2;
end




% mask2 = double(diffCurvatureScore >threshold);
% initialMask3 = maskAll.*(1 - mask2);
% gap = mask2;
% maskAllNew = initialMask3;




end