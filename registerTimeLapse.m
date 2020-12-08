% This script registers a time lapse image series using discrete Fourier
% transforms.  It requires the subfunction dftregistration.
%
% Images are iteratively registered to the previous registered image. This
% will remove any net motion of the features in the images and is more
% robust to small changes in features over time than other schemes.
%
% Images are assumed to be multi-image tif files with channels separated
% into separate files.  One channel can optionally be used as the reference
% channel for other channels.
%
% Note that the script was specficially designed for actin time lapse
% images of cells that are moving slowy, so check assumptions before using
% it for other data.

clc, clear, close all % Start with a clean workspace for reproducibility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder full of tif files to analyze
imDir = 'E:\2020-12-07_JitterCorrection_forMatt\100x Brightfield and 488\';
% Which channel should the jitter correction be peformed on?
imTag = 'Fluor'; % Set to '' to analyze all images in the folder.
% Where should the jitter corrected files be saved?
saveDir = 'E:\2020-12-07_JitterCorrection_forMatt\Registered_100x Brightfield and 488\';
% Should the script re-register images or skip over already completed work?
% overwrite = 1 will *replace* existing registered files.
% overwrite = 0 will skip files that already have a registered counterpart
overwrite = 0;
% Images are registered to within 1/usfac of a pixel.
usfac = 10; 
% If you would like other channels to registered based on the imTag chosen
% above, list their tags here.  Otherwise set otherChannels = {};
otherChannels = {'Brightf'};

%%%%%% OUTPUTS
% -- Multiple-image fit files are saved to the savedir with 'Registered_'
% added to the original image name.
% -- A .mat file name 'RegistrationInformation_' is saved for each image
% with the recorded offsets as well as original images and registered
% images for any additional comparisons.
% -- If other channels are optionally specified, multi-image tif files will
% also be saved to the savedir with ['Registered_basedOn' imTag '_'] added
% to the original tif file name.

%%%%%% DEPENDENT FUCNTIONS:
% dftregistration can be downloaded from the MATLAB file exchange:
% Efficient subpixel image registration by cross-correlation by Manuel Guizar 
% https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

%%%%%% CHANGE LOG
% 2020/12/08 RML adding cropping, commented code

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESS IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Clean up user inputs
% If no saveDir specified, save in the original directory
if isempty(saveDir) 
    saveDir = imDir;
end
% Make sure the directories end in a \ or / as appropriate
if ~strcmp(imDir(end), filesep)
    imDir = [imDir filesep];
end
if ~strcmp(saveDir(end), filesep)
    saveDir = [saveDir filesep];
end
% Make the save location if necessary
if ~exist(saveDir,'file')
    mkdir(saveDir)
end

%%%%%% Loop through each images to be registered
% Find the images to jitter correct
list = dir([imDir '*' imTag '*.tif']);

for jj = 1:length(list)    

    imname = list(jj).name(1:end-4); % Name minus .tif 
    
    % Check for already existing analysis
    if exist([saveDir 'RegistrationInformation_' imname '.mat'],'file') && ~overwrite
        disp(['Registration already complete for ' imname ' (' num2str(jj) '/' num2str(length(list)) ')'])
        continue
    end    
    
    % Number of frames to iteratively register
    N = length(imfinfo([imDir imname '.tif']));
    
    disp(['Registering ' imname ' (' num2str(N) ' frames, ' num2str(jj) '/' num2str(length(list)) ')'])

    % Set up first frame
    imOG = imread([imDir imname '.tif'],'tif',1);
    imDFT = imOG;
    
    %%%%%% Run registration
    output = NaN*ones(N,4); % Preallocate storage
    for kk = 2:N

        if ~mod(kk,50)
            disp(['   Frame ' num2str(kk) ' completed'])
        end

        % Load in the next image
        imOG(:,:,kk) = imread([imDir imname '.tif'],'tif',kk);

        f = fft2(imDFT(:,:,kk-1)); % Reference image
        g = fft2(imOG(:,:,kk)); % Current image

        % Output = information registration
        % Greg = fourier transform of the registered image
        [output(kk,:), Greg] = dftregistration(f,g,usfac);

        imDFT(:,:,kk) = abs(ifft2(Greg)); % Convert back to image space        

    end
    
    %%%%%% Cropping
    imCrop = imDFT;
    % Ceiling & floor lead to conversative cropping
    maxOff = ceil(max(output));
    minOff = floor(min(output));
    % Row limits
    r1 = max([1,maxOff(3)]); % If maxOff(3) is negative, no cropping is needed on the top, set r1 = 1
    r2 = min([size(imCrop,1),size(imCrop,1)+minOff(3)]); % If minOff(3) is positive, no cropping is needed on the bottom, set r2 = image size
    % Column limits
    c1 = max([1,maxOff(4)]);
    c2 = min([size(imCrop,1),size(imCrop,1)+minOff(4)]);
    % Crop the images    
    imCrop = imCrop(r1:r2,c1:c2,:);
    % Save the images
    imwrite(imCrop(:,:,1),[saveDir 'Registered_' imname '.tif'],'tif')
    for kk = 2:N
        % WriteMode Append makes a multi-image tif
        imwrite(imCrop(:,:,kk),[saveDir 'Registered_' imname '.tif'],'tif','WriteMode','Append')
    end
    
    %%%%%% Save registration information
    save([saveDir 'RegistrationInformation_' imname '.mat'],...
        'output','imOG','imDFT','imCrop')
    
    %%%%%% Optionally register other channels
    if ~isempty(otherChannels)
        nameParts = strsplit(imname,imTag);
        if length(nameParts) > 2
            error('Unexpected naming convention.  Expected imTag to only occur once in the image name.')
        end
        % Loop through specified other channels
        for ii = 1:length(otherChannels)
            nowTag = otherChannels{ii};
            nowName = [nameParts{1} nowTag nameParts{2}];
            disp(['   Registering ' nowName ' based on ' imTag])
            
            % Load in the first image
            nowDFT = imread([imDir nowName '.tif'],'tif',1);
            
            % Loop through each frame
            for kk = 2:N
                nowOG = imread([imDir nowName '.tif'],'tif',kk);
                buf2ft = fft2(nowOG);
                
                [nr,nc]=size(buf2ft);
                Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);                
                [Nc,Nr] = meshgrid(Nc,Nr);
                
                % Get shifts from main channel of interest
                diffphase = output(kk,2);
                row_shift = output(kk,3);
                col_shift = output(kk,4);
                
                Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
                Greg = Greg*exp(1i*diffphase);
                
                nowDFT(:,:,kk) = abs(ifft2(Greg)); % Convert back to image space
                
            end
            
            % Crop
            nowDFT = nowDFT(r1:r2,c1:c2,:);
            % Save the images
            imwrite(nowDFT(:,:,1),[saveDir 'Registered_basedOn' imTag '_' nowName '.tif'],'tif')
            for kk = 2:N
                % WriteMode Append makes a multi-image tif
                imwrite(nowDFT(:,:,kk),[saveDir 'Registered_basedOn' imTag '_' nowName '.tif'],'tif','WriteMode','Append')
            end
            
        end
        
    end

end

% Clean up and tell the user you are done
disp(' ')
disp('registerTimeLapse complete')
disp(' ')