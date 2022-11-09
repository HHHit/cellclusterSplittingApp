function im_tiff = tiffRead(file_path,digits)
info = imfinfo(file_path);
num_im = numel(info);
width = info.Width;
height = info.Height;
if digits == 16
    im_tiff = zeros(height,width,num_im,'uint16');
else
    im_tiff = zeros(height,width,num_im,'uint8');
end
for i = 1:num_im
   im_tiff(:,:,i) = imread(file_path,i);
end
end