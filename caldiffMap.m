function diffout = caldiffMap(imregion1)
diffout = zeros(size(imregion1));
localmean5 = zeros(size(imregion1));
localmean3 = zeros(size(imregion1));
for kk = 1:size(imregion1,3)
localmean5(:,:,kk) = imboxfilt(imregion1(:,:,kk),[5,5]);
localmean3(:,:,kk) = imboxfilt(imregion1(:,:,kk),[3,3]);    
end

imregion2 = imregion1(3:end-2, 3:end-2,:);
localmean2 = 25.*localmean5(3:end-2, 3:end-2,:) - 9.*localmean3(3:end-2, 3:end-2,:);
localmean2 = localmean2./(25-9);
diffx = (imregion2 - localmean2);
diffout(3:end-2, 3:end-2,:) = diffx;
end