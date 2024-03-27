%% Dual Class Processing Framework
clear all 

%% Input Data and Processing Paramaters
inputDataFile = "SyntheticDatasets/0dB_LandNew/CompleteDataset.mat"; % Synthetic Data Source
load(inputDataFile);
dataset = noisyData; % Uses data with added noise

% Parameter Selection
formatFlag = 1; % 0 = RIM Implementation, 1 = RIP Implementation
pathToInputDataBlanks = "C:\Users\mcall\Documents\19496Project\InputData\BlanksFrames\InputImage"; % Specific path to the folder and file where blank frames are stored
pathToInputDataTargets = "C:\Users\mcall\Documents\19496Project\InputData\TargetsFrames\InputImage"; % Specific path to the folder and file where the target frames are stored
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
% Delete Images in Targets Folder
myFolderTargets = "InputData/TargetsFrames"; % Generic Path to Target Folder
if isfolder(myFolderTargets) % Check if folder exists
    fileString = fullfile(myFolderTargets, '*.png'); % Search folder for target frame .png files
    filePath = dir(fileString); % Find number of png files in Target folder
    for idx = 1:1:length(filePath) 
      fileName = filePath(idx).name;% Find filename
      pathFileName = fullfile(myFolderTargets, fileName); % Generate full path to filename in target folder
      delete(pathFileName); % delete file
    end
end

% Remove Folders, Included to throw warning incase any data has been left
% over that needs to be removed. If anything left in warning appears as dont have permission to
% remove

if isfolder(myFolderBlanks)
    rmdir("InputData\BlanksFrames")
end

if isfolder(myFolderTargets)
    rmdir("InputData\TargetsFrames");
end

if isfolder("InputData\");
    rmdir("InputData\");
end

dataImage = [];
newFilename = "InputImage";
filecount = 0;

% Make and add data folders to the workspace

blanksFolder = "InputData/BlanksFrames/";
targetsFolder = "InputData/TargetsFrames/";
addpath("InputData");
mkdir(blanksFolder);
mkdir(targetsFolder);
addpath(blanksFolder);
addpath(targetsFolder);

%% Image Formatting

[xSize,ySize,zSize] = size(dataset); % Find size of the dataset

dataReal = real(dataset); % Extract real parts of data
dataImag = imag(dataset); % Extract imag parts of data

if (formatFlag == 0) % If RIM formatting selected
    dataOpt = abs(dataset); % 3rd channel = Magnitude
elseif (formatFlag == 1) % If RIP formatting selected
    dataOpt = angle(dataset); % 3rd Channel = Phase
end

for i = 1:zSize % For numebr of frames in datasetto be generated
    dataImage(:,:,1) = dataReal(:,:,i); % Real part in channel 1
    dataImage(:,:,2) = dataImag(:,:,i); % Imag part in channel 2
    dataImage(:,:,3) = dataOpt(:,:,i); % Selected fotmat data in channel 3
    filecount = string(i); % Number Image
    indFilename = "InputData/" + newFilename + filecount + ".png"; % Generate filename
    imwrite(dataImage, indFilename); % Write image to file

    if targetsInFrame(i) == 1 % If frame is blank from ground truth
        movefile(indFilename, blanksFolder);
    else % Else frame contains a target
        movefile(indFilename, targetsFolder);
    end
end

%% CNN Implementation - Based on Mathworks Example "Train Deep Learning Network to Classify New Images"

%% Create Image Data Store
imds = imageDatastore('InputData','IncludeSubfolders',true,'LabelSource','foldernames'); % Generate Image Datastore
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.7); % Split data into training and validation

%% CNN Implementation
DualClassRes18Net = resnet18; % Load Resnet-18 Network
inputSize = DualClassRes18Net.Layers(1).InputSize; % Set input to sze to ResNet-18 requirement 224x224x3
lgraph = layerGraph(DualClassRes18Net); % Convert the network to a layer graph
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


DualClassRes18Net = trainNetwork(augimdsTrain,lgraph,options); % Train network on training data with specified training options


%% Classify Validation Images
[YPred,confidenceOfFrame] = classify(DualClassRes18Net,augimdsValidation); % Perfrom Data Classification
accuracyOverall = mean(YPred == imdsValidation.Labels); % Calculate Overall Accuracy


%% Filter Frames 

infoFrames = [];
infolessFrames = [];
infoFrameCount = 0;
infolessFrameCount = 0;
[nValFrames,nValLabels] = size(confidenceOfFrame); % Find Length of validation dataset

for i = 1:nValFrames
    if confidenceOfFrame(i,2) > confidenceOfFrame (i,1) % If true Frame has target
        frameNo = strrep(imdsValidation.Files(i),pathToInputDataTargets,"");
        frameNo = strrep(frameNo,".png","");
        frameNo = str2double(frameNo); % Find frame no from filename
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataBlanks,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        infoFrameCount = infoFrameCount + 1;
        infoFrames(:,:,infoFrameCount) = noisyData(:,:,frameNo); % Add frame to predicted target frame dataset
    elseif confidenceOfFrame(i,1) > confidenceOfFrame (i,2) % If true Frame is blank
        frameNo = strrep(imdsValidation.Files(i),pathToInputDataBlanks,"");
        frameNo = strrep(frameNo,".png","");
        frameNo = str2double(frameNo); % Find frame no from filename
        if isnan(frameNo) % If true frame has been classified incorrectly
            frameNo = strrep(imdsValidation.Files(i),pathToInputDataTargets,"");
            frameNo = strrep(frameNo,".png","");
            frameNo = str2double(frameNo); % Find frame no from filename
        end
        infolessFrameCount = infolessFrameCount + 1;
        infolessFrames(:,:,infolessFrameCount) = noisyData(:,:,frameNo); % Add frame to predicted blank frame dataset
    end

end

%% Move Infoless Frames to storage

%save("InformationlessFrames.mat", "infolessFrames"); Commented out to save
%filling storage

%% Process Info Raw Data Frames to Images using Range-Doppler Algorithm
% Based on provided Simulation2D.m provided by Supervisor

if processProducts == 1 % If RDA processign has been selected

    productImages = [];
    for i = 1:1:infoFrameCount
        %% Range Reference Signal
        td0=t-2*(Xc/c);
        pha20=pi*Kr*((td0.^2)-td0*Tp);
        s0=exp(cj*pha20).*(td0 >= 0 & td0 <= Tp);
        
        fs0=fty(s0); % Reference Signal in frequency domain
        
        %% Power equalization
        amp_max=1/sqrt(2); % Maximum amplitude for equalization
        afsb0=abs(fs0);
        P_max=max(afsb0);
        
        I=find(afsb0 >= amp_max*P_max);
        fs0(I)=((amp_max*(P_max^2)*ones(1,length(I)))./afsb0(I))...
        .*exp(cj*angle(fs0(I)));
        deltaR=(2*lambda^2*(Xc).*(Ka*(dur*0.5-eta)).^2)/(8*vp^2); % RCM
        cells=round(deltaR/.56); % .56 meters/cell in range direction
        rcm_max=9; %maximum range cell migration
        fs=zeros(PRF*dur,rbins); fsm=fs; fsmb=fs; smb=fs; fsac=fs; sac=fs;
        
        %% Range Compression
        for k=1:(PRF*dur);
        fs(k,:)=fty(infoFrames(k,:,i)); % Range FFT
        fsm(k,:)=fs(k,:).*conj(fs0);% Range Matched Filtering
        smb(k,:)=ifty(fsm(k,:)); % Range IFFT
        end;
        % Range Doppler Domain
        smb0=exp(cj*pi*Ka.*eta.*(2*eta(PRF*dur/2+1)-eta));
        fsmb0=ftx(smb0); % Azimuth Matched Filter Spectrum
        %% Azimuth FFT
        for l=1:rbins;
        fsmb(:,l)=ftx(smb(:,l));  
        end;
        
        % Range Cell Migration Correction (RCMC)
        fsmb2=fsmb;
        for k=1:dur*PRF
        for m=1:rbins-rcm_max
        fsmb2(k,m)=fsmb(k,m+cells(k));
        end
        end
        for l=1:rbins
        fsac(:,l)=iftx(fsmb2(:,l)); % Azimuth IFFT
        end
        f1 = fsac; 
        % Azimuth Compression
        for l=1:rbins
        fsac(:,l)=fsmb2(:,l).*conj(fsmb0); % Azimuth Matched Filtering
        sac(:,l)=iftx(fsac(:,l)); % Azimuth IFFT / Final Target Image
        end
        productImages(:,:,i) = sac;
    end
end


%% Store Results and workspace

save("TestNewFile.mat",'confidenceOfFrame', 'accuracyOverall', 'YPred', 'imdsValidation', 'DualClassRes18Net','infoFrames', 'infolessFrames','productImages'); % Save Workspace Data for results