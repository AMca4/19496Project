%% Synthetic Data Generation
clear all
clc

%% IMAGE FORMATION PROCESS

%% Read in blank frame
I = imread("DataFolder/SyntheticDataImageFiles/BlankFrame.png");  % Read in blank (black frame for targets to be overlayed onto)
I = rgb2gray(I); % Convert frame to 3 channel image

%% Read in Template Target Images
planeImage = imread('DataFolder/SyntheticDataImageFiles/airplaneSimple.png'); % Read in plane target
tankImage = imread('DataFolder/SyntheticDataImageFiles/tankSimple.png'); % Read in tank target

%% Target Image Options
nFrames = 30; % Number of frames to be generated
shapeSelection = "filled-circle"; % Target Type in Frame, for dot dataset

[xsize,ysize] = size(I); % Size of frame
Frames = zeros(xsize,ysize); % Create Equivalent Blank Frame
for i = 1:nFrames
    Frames(:,:,i) = 0.05 + 0.05i; % Add base complex noise, required for adding complex noise later
end
%% Target Generation In Image

% Currently relies on two methods, A. Dot generation or B generic vehicle
% generation. Will produce a random integer corresponding to a
% method/target or blank and generate frame as rewuired.

% Target Type in Frame

targetsInFrame = []; % Array of target type present in frame if any

% Ground Truth Generation

for i = 1:nFrames % Iterate through nFrames (lenght of image dataset)
         value = randi([1 3]); % Generate random integer from 1 to 3 where 1 = blank, 2 = plane and 3 = tank
        targetsInFrame(i) = value; % Popuate data with target type in frame
end

finalFrames = []; % 3d matrix to contain frame data


for j = 1:nFrames

    if targetsInFrame(j) == 1
         finalFrames(:,:,j) = Frames(:,:,j); % Add blank frame to 

    elseif targetsInFrame(j) == 2 || targetsInFrame(j) == 3 % Frame contains a target

        if targetsInFrame(j) == 2 % Plane Target in frame
            targetImage = planeImage;
        elseif targetsInFrame(j) == 3 % Tank Target in frame
            targetImage = tankImage;
        end
        pos = randi(9); % Place in random position from nine options
        targetImage = rgb2gray(targetImage); % Convert target to 1 channel
        sza = size(targetImage);

        switch(pos)

            case 1    % Top Left
                
                Frames(1:sza(1),1:sza(2),j) = targetImage; 

            case 2    % Top Right

                Frames(1:sza(1),end-sza(2)+1:end,j) = targetImage; 
            case 3    % Bottom Left
               

                Frames(end-sza(1)+1:end,1:sza(2),j) = targetImage; 

            case 4    % Bottom Right

                Frames(end-sza(1)+1:end,end-sza(2)+1:end,j) = targetImage; 

            case 5

                Frames((end/2)-sza(1)+1:(end/2),end-sza(2)+1:end,j) = targetImage; 
                        
                      % Middle Left
            case 6

                 Frames((end/2)-sza(1)+1:(end/2),1:sza(2),j) = targetImage; 

            case 7   % Top Middle
                
                Frames(1:sza(1),(end/2)-sza(2)+1:(end/2),j) = targetImage; 

            case 8    % Bottom Right

                Frames(end-sza(1)+1:end,(end/2)-sza(2)+1:(end/2),j) = targetImage; 

            case 9    % Middle

                Frames((end/2)-sza(1)+1:(end/2),(end/2)-sza(2)+1:(end/2),j) = targetImage; 
        end
            
        finalFrames(:,:,j) = Frames(:,:,j);
    end

end


%% Folder Genration for Target Images
date = string(datetime("now")); % Generate Time of data generation for foldername
date = replace(date,' ',''); % String operations to make name directory friendly
date = replace(date,':',''); % String operations to make name directory friendly
foldername = strcat("SyntheticData", date); % String operations to add Synthetic data title to dataset
foldername = replace(foldername,'-','');   % String operations to make name directory friendly
mkdir(foldername);   % Generate folder directory in project
addpath(foldername); % add path to project so its in scope for future use of synthetic data generation

%% Move Dataset to Project Folder Storage

for i = 1:nFrames    %
    % imageMatrixTruecolor = cell2mat(finalFrames(i));
    imageMatrixTruecolor = finalFrames(:,:,i); % Read in frame from dataset
    imageMatrixGS = imageMatrixTruecolor*1000; % Convert from RGB image to a greyscale image (** Multiply by 1000 to put into valid range for imwrite??)
    filename = "FrameNo" + num2str(i) + ".gif"; % generate image and save as a .gif file
    imwrite(imageMatrixGS, filename);  % Write greyscale image in matrix form to gif filel generated
    movefile(filename, foldername);  % Move file to data storage folder
end


%% Synthetic SAR Data Generation

% Main SAR Parameters
PRF=300; % Pulse Repetition Frequency (Hz)
dur=3; % Time of Flight (sec), PRF*dur = received echoes
vp=200; % Velocity of platform
fo=4.5e9; % Carrier frequency (4.5GHz)
La=2; % Antenna length actual
Xc=20000; % Range distance to center of target area
X0=100; % Half Target Area Width (Target is located within [Xc-X0,Xc+X0])
Tp=.25e-5; % Chirp Pulse Duration
B0=100e6; % Baseband bandwidth is plus/minus B0\

dataFolder = foldername;

% Find Number of Images in folder
folderSearch = strcat(dataFolder,"/*.gif");
a=dir(folderSearch);
out=size(a,1);
RAWDATA = [];
snrValues = [];

for z = 1:out
    
    target_name = dataFolder + "/FrameNo" + int2str(z); % Target Image filename
    
    %% General Variables
    cj=sqrt(-1);
    c=3e8; % Propagation speed
    ic=1/c; % Propagation frequency
    lambda=c/fo; % Wavelength (60cm for fo = 4.5e9)
    eta=linspace(0,dur,PRF*dur)'; % Slow Time Array
    
    %% Range Parameters
    Kr=B0/Tp; % Range Chirp Rate
    dt=1/(2*B0); % Time Domain Sampling Interval
    Ts=(2*(Xc-X0))/c; % Start time of sampling
    Tf=(2*(Xc+X0))/c+Tp; % End time of sampling
    
    %% Azimuth Parameters
    Ka=(2*vp^2)./(lambda*(Xc)); % Linear Azimuth FM rate
    %% Measurement Parameters
    rbins=2*ceil((.5*(Tf-Ts))/dt); % Number of time (Range) samples
    t=Ts+(0:rbins-1)*dt; % Time array for data acquisition
    s=zeros(PRF*dur,rbins); % Echoed signal array
    
    %% Target Initialization
    target=imread(target_name,'gif'); %Select Input Target Profile
    [M, N] = size(target); 
    ntarget=M*N;
    
    tnum=1; xn=zeros(ntarget,1); yn=xn; Fn=xn; % Target IntializationVariables
    for m=1:M
        for n=1:N
            xn(tnum)=(n-N/2);
            yn(tnum)=(M/2-m+1);
            Fn(tnum)=double(target(m,n))/255;
            tnum=tnum+1;
        end
    end
    stretch=1;
    xn=xn*stretch; yn=yn*stretch; %Stretch out Target Profile
    
    %% GENERATE ECHOES
    for j=1:(PRF*dur)
        
        for i=1:ntarget
            
            wa=sinc(La*(atan(vp*(eta(j)-dur/2+yn(i)/vp)/Xc))/lambda).^2; % Calculate Antenna Pattern
            R=sqrt((Xc+xn(i))^2+vp^2*((eta(j)-dur/2+yn(i)/vp)^2)); % Calculate Range
            td=t-2*R/c; % Calculate Time Doppler
            s(j,:)=s(j,:)+Fn(i)*wa*exp(-cj*(4*pi*fo*ic*R)+cj*pi*Kr*(td.^2-td*Tp)).*(td >= 0 & td <= Tp); % Calculate SAR Demodulated Baseband Signal
            
        end
        
    end
    
    RAWDATA(:,:,z) = s;

    % Noise Implementation

    [xData, yData]  = size(RAWDATA(:,:,z));
    snrValues(z) = randi([-30 0]); % Calculate Random value for noise in frame
    powerArray = xData*yData;
    dataResized = reshape(RAWDATA(:,:,z),[1,powerArray]);
    signalPower = 10*log10(rms(dataResized)^2);
    noisyData(:,:,z) = awgn(RAWDATA(:,:,z),-30, signalPower); % Add noise to raw data, optionally replace SNR Value with snrValues(z) for random SNR
    

    dataName = target_name + ".png"; % Generate Filename for equivalent image
    outFolder = "SyntheticDatasets\" + dataFolder; % 
    copyfile(dataFolder, outFolder);
    final_target_name = target_name + ".gif";
    delete(final_target_name);
    delete(foldername);
end

save("CompleteDataset");
movefile("CompleteDataset.mat", outFolder);