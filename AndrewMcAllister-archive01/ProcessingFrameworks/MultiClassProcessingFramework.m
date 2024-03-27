%% Multi-Class Processing Framework

clear all 

%% Input Data and Processing Paramaters
inputDataFile = "SyntheticDatasets/0dB_LandNew/CompleteDataset.mat"; % Synthetic Data Source
load(inputDataFile);
dataset = noisyData; % Uses data with added noise

% Parameter Selection
formatFlag = 1; % 0 = RIM Implementation, 1 = RIP Implementation
pathToInputDataBlanks = "C:\Users\mcall\Documents\19496Project\InputData\BlanksFrames\InputImage"; % Specific path to the folder and file where blank frames are stored
pathToInputDataPlanes = "C:\Users\mcall\Documents\19496Project\InputData\PlanesFrames\InputImage"; % Specific path to the folder and file where the plane frames are stored
pathToInputDataTanks = "C:\Users\mcall\Documents\19496Project\InputData\TanksFrames\InputImage"; % Specific path to the folder and file where the tank frames are stored
processProducts = 1; % 0 to disable Range-Doppler Algorithm processing, 1 to enable

%% Data Folder Formatting
% Delete Images in Blank Folder
myFolderBlanks = "InputData/BlanksFrames"; % Generic Path to Blank Folder
if isfolder(myFolderBlanks) % Check if folder exists
    fileString = fullfile(myFolderBlanks, '*.png'); % Search folder for blank frame .png files
    filePath = dir(fileString); % Find number of png files in Blanks folder
    for idx = 1:1:length(filePath) 
      fileName = filePath(idx).name; % Find filename
      pathFileName = fullfile(myFolderBlanks, fileName); % Generate full path to filename in blank folder
      delete(pathFileName); % delete file
    end
end
% Delete Images in Planes Folder
myFolderPlanes = "InputData/PlanesFrames"; % Generic Path to Planes Folder
if isfolder(myFolderPlanes) % Check if folder exists
    fileString = fullfile(myFolderPlanes, '*.png'); % Search folder for plane frame .png files
    filePath = dir(fileString); % Find number of png files in plane folder
    for idx = 1:1:length(filePath) 
      fileName = filePath(idx).name;% Find filename
      pathFileName = fullfile(myFolderPlanes, fileName); % Generate full path to filename in plane folder
      delete(pathFileName); % delete file
    end
end
% Delete Images in Tanks Folder
myFolderTanks = "InputData/TanksFrames"; % Generic Path to tanks Folder
if isfolder(myFolderTanks) % Check if folder exists
    fileString = fullfile(myFolderTanks, '*.png'); % Search folder for tanks frame .png files
    filePath = dir(fileString); % Find number of png files in tanks folder
    for idx = 1:1:length(filePath) 
      fileName = filePath(idx).name;% Find filename
      pathFileName = fullfile(myFolderTanks, fileName); % Generate full path to filename in tanks folder
      delete(pathFileName); % delete file
    end
end

% Remove Folders, Included to throw warning incase any data has been left
% over that needs to be removed. If anything left in warning appears as dont have permission to
% remove

if isfolder(myFolderBlanks)
    rmdir("InputData\BlanksFrames")
end

if isfolder(myFolderPlanes)
    rmdir("InputData\PlanesFrames");
end
if isfolder(myFolderTanks)
    rmdir("InputData\TanksFrames");
end

if isfolder("InputData\");
    rmdir("InputData\");
end

dataImage = [];
newFilename = "InputImage";
filecount = 0;

% Make and add data folders to the workspace

blanksFolder = "InputData/BlanksFrames/";
planesFolder = "InputData/PlanesFrames/";
tanksFolder = "InputData/TanksFrames/";
addpath("InputData");
mkdir(blanksFolder);
mkdir(planesFolder);
mkdir(tanksFolder);
addpath(blanksFolder);
addpath(planesFolder);
addpath(tanksFolder);

%% Image Formatting

[xSize,ySize,zSize] = size(dataset); % Find size of the dataset

dataReal = real(dataset); % Extract real parts of data
dataImag = imag(dataset); % Extract imag parts of data

if (formatFlag == 0) % If RIM formatting selected
    dataOpt = abs(dataset); % 3rd channel = Magnitude
elseif (formatFlag == 1) % If RIP formatting selected
    dataOpt = angle(dataset); % 3rd Channel = Phase
end

for i = 1:zSize % For numebr of frames in data set to be generated
    dataImage(:,:,1) = dataReal(:,:,i); % Real part in channel 1
    dataImage(:,:,2) = dataImag(:,:,i); % Imag part in channel 2
    dataImage(:,:,3) = dataOpt(:,:,i); % Selected fotmat data in channel 3
    filecount = string(i); % Number Image
    indFilename = "InputData/" + newFilename + filecount + ".png"; % Generate filename
    imwrite(dataImage, indFilename); % Write image to file

    if targetsInFrame(i) == 1 % If true frame is blank from ground truth
        movefile(indFilename, blanksFolder);
    elseif  targetsInFrame(i) == 2 % If true frame is plane from ground truth
        movefile(indFilename, planesFolder);
    else % If true frame is tank from ground truth
        movefile(indFilename, tanksFolder); 
    end
end

%% CNN Implementation - Based on Mathworks Example "Train Deep Learning Network to Classify New Images"

%% Create Image Data Store
imds = imageDatastore('InputData','IncludeSubfolders',true,'LabelSource','foldernames'); % Generate Image Datastore
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.7); % Split data into training and validation

%% CNN Implementation
MultiClassRes18Net = resnet18; % Load Resnet-18 Network
inputSize = MultiClassRes18Net.Layers(1).InputSize; % Set input to sze to ResNet-18 requirement 224x224x3
lgraph = layerGraph(MultiClassRes18Net); % Convert the network to a layer graph
[learnableLayer,classLayer] = findLayersToReplace(lgraph); % Replace layers that produce output with new layers for new learning and output


%% Replace Final Layers
numClasses = numel(categories(imdsTrain.Labels)); % Number of classes in dataset (2 for Dual Class)
newLearnableLayer = fullyConnectedLayer(numClasses, 'Name','New_FC_Layer','WeightLearnRateFactor',10,'BiasLearnRateFactor',10); % Create learnable layer
lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer); % Replace old Fully connected layer with new fully conneted layer
newClassLayer = classificationLayer('Name','Output_Layer'); % Define new output layer
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer); % Replace old Output Layer with new output Layer

%% Freeze Initial Layers
layers = lgraph.Layers;
connections = lgraph.Connections;
layers(1:10) = freezeWeights(layers(1:10)); % Freeze the weights of the first ten layers to prevent overfitting to new dataset
lgraph = createLgraphUsingConnections(layers,connections); %


%% Data Augmentation
% Another step to prevent overfitting since small network is being used
pixelRange = [-30 30]; % Pixel range for translation
scaleRange = [0.9 1.1];% Scale range of pixels
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ... % Randomly flip around verical axis
    'RandXTranslation',pixelRange, ... 
    'RandYTranslation',pixelRange, ...
    'RandXScale',scaleRange, ...
    'RandYScale',scaleRange);
augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain,'DataAugmentation',imageAugmenter); % Augment training data
augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation); % Augment validation data

%% Train Network

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


MultiClassRes18Net = trainNetwork(augimdsTrain,lgraph,options); % Train network on training data with specified training options


%% Classify Validation Images
[YPred,confidenceOfFrame] = classify(MultiClassRes18Net,augimdsValidation); % Perfrom Data Classification
accuracyOverall = mean(YPred == imdsValidation.Labels); % Calculate Overall Accuracy


%% Filter Frames 

planeFrames = [];
tankFrames = [];
infoFrameCount = 0;
infolessFrames = [];
planeFrameCount = 0;
tankFrameCount = 0;
infolessFrameCount = 0;
[nValFrames,nValLabels] = size(confidenceOfFrame); % Find Length of validation dataset

for i = 1:nValFrames
    if (confidenceOfFrame(i,2) > confidenceOfFrame (i,1)) && (confidenceOfFrame(i,2) > confidenceOfFrame (i,3))  % If true Frame has target
        frameNo = strrep(imdsValidation.Files(i),pathToInputDataPlanes,"");
        frameNo = strrep(frameNo,".png","");
        frameNo = str2double(frameNo); % Find frame no from filename
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataBlanks,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
         if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataTanks,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        planeFrameCount = planeFrameCount + 1;
        planeFrames(:,:,planeFrameCount) = noisyData(:,:,frameNo); % Add frame to predicted target frame dataset
    elseif (confidenceOfFrame(i,1) > confidenceOfFrame (i,2)) && (confidenceOfFrame(i,1) > confidenceOfFrame (i,3)) % If true Frame is blank
        frameNo = strrep(imdsValidation.Files(i),pathToInputDataBlanks,"");
        frameNo = strrep(frameNo,".png","");
        frameNo = str2double(frameNo); % Find frame no from filename
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataPlanes,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataTanks,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        infolessFrameCount = infolessFrameCount + 1;
        infolessFrames(:,:,infolessFrameCount) = noisyData(:,:,frameNo); % Add frame to predicted blank frame dataset
    else % If true Frame is tank
        frameNo = strrep(imdsValidation.Files(i),pathToInputDataTanks,"");
        frameNo = strrep(frameNo,".png","");
        frameNo = str2double(frameNo); % Find frame no from filename
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataPlanes,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataBlanks,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        tankFrameCount = tankFrameCount + 1;
        tankFrames(:,:,tankFrameCount) = noisyData(:,:,frameNo); % Add frame to predicted blank frame dataset
    end

end

save("TestNewFileMC.mat",'confidenceOfFrame', 'accuracyOverall', 'YPred', 'imdsValidation', 'MultiClassRes18Net','tankFrames','planeFrames', 'infolessFrames'); % Save Workspace Data for results