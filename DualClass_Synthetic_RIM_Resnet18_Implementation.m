
%%RESNET TEST LIVESCRIPT


%% Sort Folders 

% Clear Input Data Folders and recreate
clear all
% Specify the folder where the files live.
myFolderBlanks = "InputData/BlanksFrames";
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if isfolder(myFolderBlanks)
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolderBlanks, '*.png'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(myFolderBlanks, baseFileName);
      delete(fullFileName);
    end
end

myFolderTargets = "InputData/TargetsFrames";
% Check to make sure that folder actually exists
if isfolder(myFolderTargets)
    % Get a list of all files in the folder with the desired file name
    filePattern = fullfile(myFolderTargets, '*.png'); 
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(myFolderTargets, baseFileName);
      delete(fullFileName);
    end
end



% Clean folders, all included incase of old data present from other
% implementation

if isfolder(myFolderBlanks)
    rmdir("InputData\BlanksFrames")
end
% if isfolder(myFolderDots)
%     rmdir("InputData\DotsFrames");
% end
% if isfolder(myFolderTanks)
%     rmdir("InputData\TanksFrames");
% end
% if isfolder(myFolderPlanes)
%     rmdir("InputData\PlanesFrames");
% end

if isfolder(myFolderTargets)
    rmdir("InputData\TargetsFrames");
end

if isfolder("InputData\");
    rmdir("InputData\");
end



imageDataset = {};
dataImage = [];
newFilename = "InputImage";
filecount = 0;

% Make and add folders to path


blanksFolder = "InputData/BlanksFrames/";
targetsFolder = "InputData/TargetsFrames/";
addpath("InputData");
mkdir(blanksFolder);
mkdir(targetsFolder);
addpath(blanksFolder);
addpath(targetsFolder);


%% Read in synthetic data and format as a 3d image (RGB)

load("SyntheticDatasets\SyntheticData02Mar2024013607\CompleteDataset.mat"); % Select dataset to be used

dataset = noisyData;

[xSize,ySize,zSize] = size(dataset);

dataReal = real(dataset);
dataImag = imag(dataset);
dataMag = abs(dataset);

for i = 1:zSize
    dataImage(:,:,1) = dataReal(:,:,i);
    dataImage(:,:,2) = dataImag(:,:,i);
    dataImage(:,:,3) = dataMag(:,:,i);
    %imageDataset(i) = {dataImage};
    filecount = string(i);
    indFilename = "InputData/" + newFilename + filecount + ".png";
    imwrite(dataImage, indFilename);

    if targetsInFrame(i) == 1
        movefile(indFilename, blanksFolder);
    else 
        movefile(indFilename, targetsFolder);
    end
end


save("ImageDataset.mat", "imageDataset");

%imds = imageDatastore("ImageDataset.mat");
%[imdsTrain, imdsValidation] = splitEachLabel(imds,0.7);


%% Create Image Data Store
imds = imageDatastore('InputData', ...
    'IncludeSubfolders',true, ...
    'LabelSource','foldernames'); 
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.7);


%% CNN Implementation
DualClassRes18Net = resnet18;
%net = googlenet;
inputSize = DualClassRes18Net.Layers(1).InputSize;

lgraph = layerGraph(DualClassRes18Net);

[learnableLayer,classLayer] = findLayersToReplace(lgraph);



%% Replace Final Layers


numClasses = numel(categories(imdsTrain.Labels));

if isa(learnableLayer,'nnet.cnn.layer.FullyConnectedLayer')
    newLearnableLayer = fullyConnectedLayer(numClasses, ...
        'Name','new_fc', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
    
elseif isa(learnableLayer,'nnet.cnn.layer.Convolution2DLayer')
    newLearnableLayer = convolution2dLayer(1,numClasses, ...
        'Name','new_conv', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
end

lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);

newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);

figure('Units','normalized','Position',[0.3 0.3 0.4 0.4]);
plot(lgraph)
ylim([0,10])



%% Freeze Initial Layers


layers = lgraph.Layers;
connections = lgraph.Connections;

layers(1:10) = freezeWeights(layers(1:10));
lgraph = createLgraphUsingConnections(layers,connections);


%% Train Network
pixelRange = [-30 30];
scaleRange = [0.9 1.1];
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandXTranslation',pixelRange, ...
    'RandYTranslation',pixelRange, ...
    'RandXScale',scaleRange, ...
    'RandYScale',scaleRange);
augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain, ...
    'DataAugmentation',imageAugmenter);

augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation);


miniBatchSize = 10;
valFrequency = floor(numel(augimdsTrain.Files)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',6, ...
    'InitialLearnRate',3e-4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',augimdsValidation, ...
    'ValidationFrequency',valFrequency, ...
    'Verbose',false, ...
    'Plots','training-progress');


DualClassRes18Net = trainNetwork(augimdsTrain,lgraph,options);


%% Classify Validation Images
[YPred,probs] = classify(DualClassRes18Net,augimdsValidation);
accuracy = mean(YPred == imdsValidation.Labels);


idx = randperm(numel(imdsValidation.Files),4);
figure
for i = 1:4
    subplot(2,2,i)
    I = readimage(imdsValidation,idx(i));
    imshow(I)
    label = YPred(idx(i));
    title(string(label) + ", " + num2str(100*max(probs(idx(i),:)),3) + "%");
end
