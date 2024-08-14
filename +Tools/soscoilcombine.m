function Image = soscoilcombine(Image_Stack)

%Function to take a stack of images 
%(Size =ImageSize(1),ImageSize(2),ImageSize(3),NCoils)
%Return sum of squares combination of the images following GPI's
%implementation of this

if ndims(Image_Stack) == 4
    nCoil = size(Image_Stack,4);

    tempSOS = zeros(size(Image_Stack));
    for i = 1:nCoil
        tempSOS(:,:,:,i) = abs(Image_Stack(:,:,:,i).^2);
    end

    Image = squeeze(sqrt(sum(tempSOS,4)));
elseif ndims(Image_Stack) == 3
    nCoil = size(Image_Stack,3);
    
    tempSOS = zeros(size(Image_Stack));
    for i = 1:nCoil
        tempSOS(:,:,i) = abs(Image_Stack(:,:,i).^2);
    end

    Image = squeeze(sqrt(sum(tempSOS,3)));
end