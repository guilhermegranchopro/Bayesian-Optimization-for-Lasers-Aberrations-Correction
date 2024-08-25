clc
clear all
close all

%%
H=1280; V=1024; %%Number of Horizontal and Vertical pixels
PixelSize = 12.5e-3; % [mm] pixel
x = PixelSize.*(-H/2:(H/2-1)); %[mm] x width of window
y = PixelSize.*(-V/2:(V/2-1)); %[mm] y width of window
[X,Y] = meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);

Holo = zeros(V,H); Holo_Blank = zeros(V,H);
Nx = 220; Ny = 0;%%Number of horizontal and vertical grooves

%CAMERA STUFF
% Load TLCamera DotNet assembly. The assembly .dll is assumed to be in the 
% same folder as the scripts.
NET.addAssembly([pwd, '\Thorlabs.TSI.TLCamera.dll']);
disp('Dot NET assembly loaded.');

tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;

% Get serial numbers of connected TLCameras.
serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
disp([num2str(serialNumbers.Count), ' camera was discovered.']);

if (serialNumbers.Count > 0)

    % Open the first TLCamera using the serial number.
    disp('Opening the first camera')
    cam = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);

     % Check if the camera supports setting "Gain"
    gainRange = cam.GainRange;
    if (gainRange.Maximum > 0)
        cam.Gain = 0;
    end
    
    % Set the FIFO frame buffer size. Default size is 1.
    cam.MaximumNumberOfFramesToQueue = 1;

    fig=figure(1);

    % Start software triggered image acquisition
    disp('Starting software triggered image acquisition.');

    % Set the number of frames per software trigger and start trigger
    % acquisition
    cam.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
    cam.FramesPerTrigger_zeroForUnlimited = 1;
    cam.Arm;

    %GET CORRECT EXPOSURE TIME
    
    cam.ExposureTime_us = 200000;  %define a large exposure to begin
    maxPixelValue = double(2^cam.BitDepth - 3);  

    %get frame
    cam.IssueSoftwareTrigger;
    % Wait for image buffer to be filled to prevent sending too many
    % software triggers.
    while (cam.NumberOfQueuedFrames == 0)
        pause(0.01);
    end
    imageFrame = cam.GetPendingFrameOrNull; 
    if ~isempty(imageFrame)
        imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);
        disp("ok")
    end

    %find correct value for exposure
    exposure = cam.ExposureTime_us;
    while(max(imageData) >= maxPixelValue)
        cam.ExposureTime_us = exposure/2;
        exposure = exposure/2;
        disp(cam.ExposureTime_us);
        %take new frame
        cam.IssueSoftwareTrigger;
        while (cam.NumberOfQueuedFrames == 0)
            pause(0.01);
        end
        imageFrame = cam.GetPendingFrameOrNull;
        if ~isempty(imageFrame)
            imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);
        end     
    end
    disp('Exposure time in  microseconds: ');
    disp(cam.ExposureTime_us);
    

    %START ACQUISITION

    for ai=0:255   %for the aquisition

        a = ai/255; %between 0 and 1
        Holo = zeros(V,H);
        for i=1:H
            for j=1:V
                if i > H/2 
                    Holo(j,i)=a;
                end
            end
        end 
        % Grascale normalization from [0, 2Pi]to [0 255]
        Holo = Holo - min(Holo(:)); %Make all pixel values positive
        HoloSLM = Holo;
        device_number=2;
        fullscreen(HoloSLM, device_number)

        %get frame
        pause(1); %pausa de 1 segundo para dar tempo ao holograma de chegar ao slm e ao feixe de reagir ao novo holograma
        cam.IssueSoftwareTrigger; %carregar no botao de tirar foto
        while (cam.NumberOfQueuedFrames == 0)
            pause(0.01);
        end

        % Get the pending image frame. 
        imageFrame = cam.GetPendingFrameOrNull;  %objeto do tipo Thorlabs.TSI.TLCameraInterfaces.Frame
        if ~isempty(imageFrame)

            % For color images, the image data is in BGR format.
            imageData = imageFrame.ImageData.ImageData_monoOrBGR;

            disp(['Image frame number: ' num2str(imageFrame.FrameNumber)]);

            % TODO: custom image processing code goes here
            imageHeight = imageFrame.ImageData.Height_pixels;
            imageWidth = imageFrame.ImageData.Width_pixels;

            imageData2D = reshape(uint16(imageData), [imageWidth, imageHeight]);
            fig; imagesc(imageData2D'), colormap(gray), colorbar 
            str = "calib"+ ai + ".png";
            str2 = "calib"+ ai + ".tif";
            saveas(fig,str);
            imwrite(imageData2D,str2);

        end

    end %for

    % Stop software triggered image acquisition
    disp('Stopping software triggered image acquisition.');
    cam.Disarm;

    % Release the camera
    disp('Releasing the camera');
    cam.Dispose;
    delete(cam);

end %if

% Release the serial numbers
delete(serialNumbers);

% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);

