function diffout = caldiffMap_v2(imregion1)
diffout = zeros(size(imregion1));
localmean5 = zeros(size(imregion1));
localmean3 = zeros(size(imregion1));
for kk = 1:size(imregion1,3)
localmean3(:,:,kk) = imboxfilt(imregion1(:,:,kk),[3,3]);    
end

imregion2 = imregion1(3:end-2, 3:end-2,:);
localmean2 = localmean3(3:end-2, 3:end-2,:).*9 - imregion2;
diffx = (imregion2 - localmean2./8);
diffout(3:end-2, 3:end-2,:) = diffx;
end