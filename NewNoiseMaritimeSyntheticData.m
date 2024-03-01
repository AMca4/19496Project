clear all
clc

%% IMAGE FORMATION PROCESS

%% Read in blank frame
I = imread("DataFolder/SyntheticDataImageFiles/BlankFrame.png");  % Read in blank (black frame for targets to be overlayed onto)
I = rgb2gray(I);

%% Read in Template Target Images

planeImage = imread('DataFolder/SyntheticDataImageFiles/airplaneSimple.png');
tankImage = imread('DataFolder/SyntheticDataImageFiles/tankSimple.png');

%% Target Image Options
nFrames = 5; % Number of frames to be generated
nTargets = 1; % Max Number of Targets present in each frame
shapeSelection = "filled-circle"; % Target Type in Frame

[xsize,ysize] = size(I);
Frames = zeros(xsize,ysize);
for i = 1:nFrames
    Frames(:,:,i) = I;
end
size(j)
%% Target Generation In Image

% Currently relies on two methods, A. Dot generation or B generic vehicle
% generation. Will produce a random integer corresponding to a
% method/target or blank and generate frame as rewuired.





% Target Type in Frame

targetsInFrame = []; % Array of target type present in frame if any

for i = 1:nFrames % Iterate through nFrames (lenght of image dataset)
        %value = randi([0 3]); % Generate random integer from 0 to 3 where 0 = blank, 1 = dot, 2 = plane and 3 = tank
        value = randi([1 3]); % Generate random integer from 0 to 3 where 0 = blank, 1 = dot, 2 = plane and 3 = tank
        targetsInFrame(i) = value; % Popuate data with target type in frame
end
     
    % finalFrames = {}; % 3d cell array to contain frame data
    finalFrames = []; % 3d cell array to contain frame data
    baseFrames = [];
    snrValues = [];
    %lengthA = length(targetsInFrame); % Number of frames in dataset to be generated

    for j = 1:nFrames

        if targetsInFrame(j) == 0
            %noisyFrame = imnoise(Frames(:,:,j),"gaussian",0, 0.01);
            %noisyFrame = imnoise(Frames(:,:,j),"gaussian",0, 0.1);
            %baseFrames(:,:,j) = Frames(:,:,j);
            %snrValues(j) = snr(baseFrames(:,:,j),noisyFrame);
            %finalFrames(:,:,j) = noisyFrame;
            finalFrames(:,:,j) = Frames(:,:,j);
            
        elseif targetsInFrame(j) == 1
            frame = Frames(:,:,j); % Variable for individual frame to be operated on
            RGB =  Frames(:,:,j);
            plainFrame = 0;
            x = randi(xsize); % Genrated random x postion on frame
            y = randi(ysize); % generate random y postion on frame
            sizeShape = 1; % Size of Shape Input (INVESTIGATE MEANING OF THIS INPUT, RATIO TO FRAME ETC????)
            RGB = insertShape(frame,"filled-circle",[x y sizeShape],ShapeColor="white",Opacity=1); % I to be modified to be currently modified frame, needs to be modified n times then saves as new PNG
            baseFrames(:,:,j) = rgb2gray(RGB);
            %noisyFrame = imnoise(RGB,"gaussian",0, 0.01);
            %noisyFrame = imnoise(RGB,"gaussian",0, 0.1);
            %snrValues(j) = snr(baseFrames(:,:,j),rgb2gray(noisyFrame));
            %finalFrames(:,:,j) = rgb2gray(noisyFrame); % Insert Frame with targets into dataset 
            finalFrames(:,:,j) = baseFrames(:,:,j);

        elseif targetsInFrame(j) == 2 || targetsInFrame(j) == 3

            if targetsInFrame(j) == 2
                targetImage = planeImage;
            elseif targetsInFrame(j) == 3
                targetImage = tankImage;
            end
            pos = randi(9);
                % if shapeSelection == "plane"
                %     targetImage = planeImage;
                % elseif shapeSelection == "tank"
                %     targetImage = tankImage;
                % end
                % simply insert the smaller image into the SE corner by direct
                % indexing Random position algorithm to be completed
                targetImage = rgb2gray(targetImage);
                sza = size(targetImage);
    
                switch(pos)
    
                    case 1    % Top Left
                        
                        Frames(1:sza(1),1:sza(2),j) = targetImage; % insert into ROI
    
                    case 2    % Top Right
    
                        Frames(1:sza(1),end-sza(2)+1:end,j) = targetImage; % insert into ROI
    
                    case 3    % Bottom Left
                       
    
                        Frames(end-sza(1)+1:end,1:sza(2),j) = targetImage; % insert into ROI
    
                    case 4    % Bottom Right
    
                        Frames(end-sza(1)+1:end,end-sza(2)+1:end,j) = targetImage; % insert into ROI

                    case 5

                        Frames((end/2)-sza(1)+1:(end/2),end-sza(2)+1:end,j) = targetImage; % insert into ROI) = targetImage;
                                
                              % Middle Left
                    case 6

                         Frames((end/2)-sza(1)+1:(end/2),1:sza(2),j) = targetImage; % insert into ROI

                    case 7   % Top Middle
                        
                        Frames(1:sza(1),(end/2)-sza(2)+1:(end/2),j) = targetImage; % insert into ROI

                    case 8    % Bottom Right
    
                        Frames(end-sza(1)+1:end,(end/2)-sza(2)+1:(end/2),j) = targetImage; % insert into ROI
    
                    case 9    % Middle
    
                        Frames((end/2)-sza(1)+1:(end/2),(end/2)-sza(2)+1:(end/2),j) = targetImage; % insert into ROI
                end
                
                %noisyFrame = imnoise(Frames(:,:,j),"gaussian",0, 0.01);
                %noisyFrame = imnoise(Frames(:,:,j),"gaussian",0, 0.1);
                %baseFrames(:,:,j) = Frames(:,:,j);
                %snrValues(j) = snr(baseFrames(:,:,j),noisyFrame);
                %finalFrames(:,:,j) = noisyFrame;
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


%% Ground Truth Data Generation

% Count Target Occurences For random Dot Case
shapeType =1;
if shapeType == 1
    count1 = 0;
    count2 = 0;
    count3 = 0;

    for z = 1:nFrames
        nDots = targetsInFrame(z);
        switch nDots
            case 1
                count1 = count1 +1;
            case 2
                count2 = count2 +1;
            case 3
                count3 = count3 +1;
        end
    end
    dotCount = [count1, count2, count3];
    

    % N Targets in Frame for Plane Data
elseif shapetype == 2
        for z = 1:nFrames
            targetPresent = targetsInFrame(z);
        end

end





% Synthetic SAR Data Generation

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

% Target

%target_name='pole1'; %Name of Target Profile Image (GIF GrayscaleImage)
%target_name='singlePS';
%target_name= 'airplaneSimple';

% Query Folder containing Images to be processed

prompt = "Name of Folder containing Images: ";
dataFolder = foldername;

% Find Number of Images in folder
folderSearch = strcat(dataFolder,"/*.gif");
a=dir(folderSearch);
out=size(a,1);
rawData = [];
noisyData = [];
SNRTestValues = [];

for z = 1:out
    
    target_name = dataFolder + "/FrameNo" + int2str(z);

    %%Noise
    noise=0; % Set this flag to add noise to signal
    std_dev=.2; % Standard Deviation of Noise
    
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
    % Measurement Parameters
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
    for j=1:(PRF*dur);
        
        for i=1:ntarget;
            
            wa=sinc(La*(atan(vp*(eta(j)-dur/2+yn(i)/vp)/Xc))/lambda).^2;
            
            R=sqrt((Xc+xn(i))^2+vp^2*((eta(j)-dur/2+yn(i)/vp)^2));
            
            td=t-2*R/c;
            
            if noise==1
                s(j,:)=s(j,:)+std_dev*randn(size(s(j,:)))...
                    +Fn(i)*wa*exp(-cj*(4*pi*fo*ic*R)+cj*pi*Kr*...
                    (td.^2-td*Tp)).*(td >= 0 & td <= Tp);
                
            else
                s(j,:)=s(j,:)+Fn(i)*wa*exp(-cj*(4*pi*fo*ic*R)+cj*pi*Kr*...
                    (td.^2-td*Tp)).*(td >= 0 & td <= Tp);
            
            end
            
        end
        
    end
    
    %rawDataReal = real(s);
    %rawDataImag = imag(s);

    rawData(:,:,z) = s;


    %% Noise Implementation

    % noisyData = awgn(rawData, 5);
    % noise = noisyData - rawData;
    % snrCheck = snr(rawData, noise);

    [xData, yData]  = size(rawData(:,:,z));

    powerArray = xData*yData;
    dataResized = reshape(rawData(:,:,z),[1,powerArray]);
    signalPower = 10*log10(rms(dataResized)^2);
    noisyData = awgn(rawData(:,:,z), -15, signalPower);
    
    
    noise = noisyData - rawData(:,:,z);
    snrValues(z) = snr(rawData(:,:,z), noise);

    
    %figure(z), imagesc(real(s))
    %xlabel('Range, samples'), ylabel('Azimuth, samples')
    %title('Raw Data'), colormap('gray');

    %save("CompleteDataset");

    dataName = target_name + ".png";
   % exportgraphics(figure, dataName);
    outFolder = "SyntheticDatasets\" + dataFolder;
    copyfile(dataFolder, outFolder);
    %movefile("CompleteDataset.mat", outFolder);
    %(dataName, outFolder);
    final_target_name = target_name + ".gif";
    delete(final_target_name);
    delete(foldername);
    %pause();

end

save("CompleteDataset");
movefile("CompleteDataset.mat", outFolder);