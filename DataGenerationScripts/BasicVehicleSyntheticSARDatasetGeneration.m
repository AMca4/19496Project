clear all
clc

%% IMAGE FORMATION PROCESS

%% Read in blank frame
I = imread("BlankFrame.png");  % Read in blank (black frame for targets to be overlayed onto)
I = rgb2gray(I);

%% Target Image Options

nFrames = 5; % Number of frames to be generated
nTargets = 1; % Max Number of Targets present in each frame
shapeSelection = "tank"; % Target Type in Frame, Options are "plane" or "tank"

[xsize,ysize] = size(I);
Frames = zeros(xsize,ysize);
for i = 1:nFrames
    Frames(:,:,i) = I;
end
size(j)
%% Target Generation In Image

 % Generate random nTargets for nFrames
    targetsInFrame = []; % Array containing N targets present in each individual frame
    
    
    for i = 1:nFrames % Iterate through nFrames (lenght of image dataset)
             value = randi([0 1]); % Generate random integer from 1 to 3
            targetsInFrame(i) = value; % Popuate data with number of targets present in each frame
    end


    % finalFrames = {}; % 3d cell array to contain frame data
    finalFrames = []; % 3d cell array to contain frame data
    lengthA = length(targetsInFrame); % Number of frames in dataset to be generated

    for j = 1:nFrames 
      %if (shapeSelection == "plane")
        %for a=1:lengthA
            if targetsInFrame(j) ==1
                pos = randi(4);
                if shapeSelection == "plane"
                    targetImage = imread('airplaneSimple.png');
                elseif shapeSelection == "tank"
                    targetImage = imread('tankSimple.png');
                end
                % simply insert the smaller image into the SE corner by direct
                % indexing Random position algorithm to be completed
                targetImage = rgb2gray(targetImage);
                sza = size(targetImage);
    
                switch(pos)
    
                    case 1    % Top Left
                        
                        Frames(1:sza(1),end-sza(2)+1:end,j) = targetImage; % insert into ROI
    
                    case 2    % Top Right
    
                         Frames(1:sza(1),end-sza(2)+1:end,j) = targetImage; % insert into ROI
    
                    case 3    % Bottom Left
                       
    
                        Frames(end-sza(1)+1:end,1:sza(2),j) = targetImage; % insert into ROI
    
                    case 4    % Bottom Right
    
                        Frames(end-sza(1)+1:end,end-sza(2)+1:end,j) = targetImage; % insert into ROI
                end
                
                finalFrames(:,:,j) = Frames(:,:,j);
            end
        %end
    end
    %end

%% Folder Genration for Target Images

date = string(datetime("now")); % Generate Time of data generation for foldername
date = replace(date,' ',''); % String operations to make name directory friendly
date = replace(date,':',''); % String operations to make name directory friendly
foldername = strcat("SyntheticData", date); % String operations to add Synthetic data title to dataset
foldername = replace(foldername,'-','');   % String operations to make name directory friendly
mkdir(foldername);   % Generate folder directory in project
addpath(foldername); % add path to project so its in scope for future use of synthetic data generation

%% Move Dataset to Project Folder Storage

for i = 1:lengthA    %
    % imageMatrixTruecolor = cell2mat(finalFrames(i));
    imageMatrixTruecolor = finalFrames(:,:,i); % Read in frame from dataset
    imageMatrixGS = imageMatrixTruecolor*1000; % Convert from RGB image to a greyscale image (** Multiply by 1000 to put into valid range for imwrite??)
    filename = "FrameNo" + num2str(i) + ".gif"; % generate image and save as a .gif file
    imwrite(imageMatrixGS, filename);  % Write greyscale image in matrix form to gif filel generated
    movefile(filename, foldername);  % Move file to data storage folder
end


%% Ground Truth Data Generation

% 
% % N Targets in Frame for Plane Data
% for z = 1:lengthA
%     targetPresent = targetsInFrame(z);
% end

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
RAWDATA = [];

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
    
    rawDataReal = real(s);
    rawDataImag = imag(s);
    
    RAWDATA(:,:,z) = s;
  
    figure(z), imagesc(real(s))
    xlabel('Range, samples'), ylabel('Azimuth, samples')
    title('Raw Data'), colormap('gray');

    %save("CompleteDataset");

    dataName = target_name + ".png";
   % exportgraphics(figure, dataName);
    outFolder = "SyntheticDatatsets\" + dataFolder;
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