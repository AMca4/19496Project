%% RES-NET 18 NEURAL NETWORK

% Initially store data into a 3D matrix of (Real, Imag, Phase)

dataset = [1+1i,2-2i,3; 4,5+9i,6;7,8,9-8i];
dataset(:,:,2) = [10,11+9i,13;14-2i,15,16+1i;17-4i,18,19+6i];

dataReal = real(dataset);
dataImag = imag(dataset);
dataPhase = angle(dataset);

% Perform any manipulation of data that may be required, resising etc


% Restructure Dataset for training/validation

imds = imageDatastore('');
[imdsTrain, imdsValidation] = splitEachLabel(imds,0.7);


% Implement ResNet

net = deepNetworkDesigner(resnet18);


inputSize = net.Layers(1).InputSize;
