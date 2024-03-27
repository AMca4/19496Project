%% Image Slicing of Umbra GEC Data

clear all

%% Open Tiff Object
t = Tiff('2024-01-20-12-26-56_UMBRA-04_GEC.tif'); % Create Tiff Object
imageData = read(t); % Read the image data contained in the Tiff Object
editImage = imageData; % Create a copy of image to perform slicing on
% Display the image
%imshow(imageData);
close(t);

%% Slice Image

focusedSection = editImage(3001:4200,6701:8000); % Target area of image containing desired targets
imshow(focusedSection); % Show section of dataset

images=[]; % Array to hold frames
count = 1;

% Since Image size is equivalent to 1300x1200 = 13x12 of 100x100 pixel
% images
nColumns = 13;
nRows = 12;

for j = 1:1:nColumns % Extract 100x100 Sliced Frame from focused area column j
    for i = 1:1:nRows % Extract 100x100 Sliced Frame from focused area row i
        idxXEnd=i*100;
        idxYEnd=j*100;
        idxXStart=idxXEnd-99;
        idxYStart=idxYEnd-99;
        fracImage=focusedSection(idxXStart:idxXEnd, idxYStart:idxYEnd); % Generate 100x100 frame from indexes
        images(:,:,count) = fracImage; % Add image to frame data array
        count = count+1;
    end
end

%% Write Images to files

nFrames = 156; % NFrames to be produced

for z = 1:1:nFrames
    filename = "Images/ImageNo" + int2str(z)+".gif";
    imwrite(images(:,:,z),filename);
end