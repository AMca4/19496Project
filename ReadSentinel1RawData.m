clc; close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%  Index Data Component 
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wd = cd;
% Go to the location where you can get file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
location1=wd;
% Enter the number of files you need to process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(location1);  
    % Get the file
    %%%%%%%%%%%%%%%%%%%%%
cd(wd)  
      ss='Pick the .dat file';   
     [filename1, pathname1] = uigetfile('*_index.mat', ss);
     if isequal(filename1,0) || isequal(pathname1,0);    disp('User pressed cancel');      end
fullname=fullfile(pathname1, filename1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullname,'dataUnitsOffset','NIm')
    % 
    fullname=fullfile(pathname1, [filename1(1:40),'_Header.mat']);
    load(fullname,'NSAMP')
    load(fullname,'TableHeader')
    Nr = max(NSAMP);
    NIms=NIm;   %% NIm
    Na = dataUnitsOffset(NIms);
    echo = single(zeros(Nr,Na)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use clock tick to measure CPU time
inittic = tic;
updr = inittic;
updateTime=100;  % time between backprojection processing status updates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
IXAZ =find(TableHeader.DFT==4);
ab=IXAZ(1) ;
%
for jj = 1:NIms-1
    %
   sst = dataUnitsOffset(jj);
    if jj~=NIm-1
   een = dataUnitsOffset(jj+1)-1;
    else
   een = dataUnitsOffset(jj+1);    
    end
    % 
    fullname=fullfile(pathname1, [filename1(1:40),'_',num2str(sst),'_',num2str(een),'.mat']);
    load(fullname)
    echo(1:size(complex_Svalue,1),sst:een) = complex_Svalue;
    clear complex_Svalue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% updated time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            tdr = toc(updr);        % plot an updated image
            if tdr > updateTime
                fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',een,Na,toc(inittic)/60,(toc(inittic)*(Na-een)/een)/60)
                updr = tic; 
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
  fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',een,Na,toc(inittic)/60,(toc(inittic)*(Na-een)/een)/60)
echo = echo(:,IXAZ); %%% start from 409 : end
fullname=fullfile(pathname1, [filename1(1:40),'_echo.mat']);
disp('Save data')
save(fullname,'echo','ab','IXAZ')
%
disp('display data in time domain')
imagesc(abs(echo));colormap jet
% pause
% 
% % plot(abs(fftshift(fft(echo(:,fix(Na/3))))))
% %imagesc(abs(fftshift( fft(echo,[],1),1 )));colormap jet
% frecho= zeros(size(echo));
% disp('display data in freq domain in range')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % use clock tick to measure CPU time
% inittic = tic;
% updr = inittic;
% updateTime=100;  % time between backprojection processing status updates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for a=1:Na
%  frecho(:,a) =  abs(fftshift(fft(echo(:,a))));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% updated time
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%             tdr = toc(updr);        % plot an updated image
%             if tdr > updateTime
%                 fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',a,Na,toc(inittic)/60,(toc(inittic)*(Na-a)/a)/60)
%                 updr = tic; 
%             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',a,Na,toc(inittic)/60,(toc(inittic)*(Na-a)/a)/60)
% imagesc(frecho);colormap jet
% clear frecho
% pause
% 
% % plot(abs(fftshift(fft(echo(fix(Nr/2),:)))))
% % disp('display data in freq domain in azimuth')
% % imagesc(abs(fftshift( fft(echo,[],2),2 )));colormap jet
% %
% faecho= zeros(size(echo));
% disp('display data in freq domain in azimuth')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % use clock tick to measure CPU time
% inittic = tic;
% updr = inittic;
% updateTime=100;  % time between backprojection processing status updates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for a=1:Nr
%  faecho(a,:) =  abs(fftshift(fft(echo(a,:))));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% updated time
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%             tdr = toc(updr);        % plot an updated image
%             if tdr > updateTime
%                 fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',a,Nr,toc(inittic)/60,(toc(inittic)*(Nr-a)/a)/60)
%                 updr = tic; 
%             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% fprintf(1,'%d of %d complete, %.3G min elapsed, %8.3G min estimated\n',a,Nr,toc(inittic)/60,(toc(inittic)*(Nr-a)/a)/60)
% imagesc(faecho);colormap jet
% clear frecho
% 
% pause
