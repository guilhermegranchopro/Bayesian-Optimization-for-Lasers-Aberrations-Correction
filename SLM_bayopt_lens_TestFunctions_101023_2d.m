clear all
close all

% addpath()

%Bayesian optimisation of a focal spot using SLM and Zernika polynomial
%corrections

%% Set-up the camera
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
end
cam.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered; 
cam.IssueSoftwareTrigger;
cam.ExposureTime_us = 1000;  %define a large exposure to begin
cam.BitDepth
maxPixelValue = double(2^cam.BitDepth - 3);
%% Main loop
Zern_amplitudes_1 = optimizableVariable("Zernikamplitudes_foc",[-0.7,-0.5],"Type","real");
Zern_amplitudes_2 = optimizableVariable("Zernikamplitudes_astig1",[0.25,0.35],"Type","real");

max_count_h = @(x)peak_intensity(x,cam);
results = bayesopt(max_count_h,[Zern_amplitudes_1,Zern_amplitudes_2]);
 
%% Close camera
disp('Releasing the camera');
cam.Dispose;
delete(cam);
% Release the serial numbers
delete(serialNumbers);
% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);
%% Objective function
function max_count = peak_intensity(xx,cam)%TAKE PICTURE
    %% SLM and laser parameters
    H=1280; V=1024; %%Number of Horizontal and Vertical pixels
    PixelSize = 12.5e-3; % [mm] pixel
    x = PixelSize.*(-H/2:(H/2-1)); %[mm] x width of window
    y = PixelSize.*(-V/2:(V/2-1)); %[mm] y width of window
    [X,Y] = meshgrid(x,y);
    [PHI,RHO] = cart2pol(X,Y);
    
    Holo = zeros(V,H); Holo_Blank = zeros(V,H);
    lambda=660e-6;
    %% Hologram part
    ff =- 5000;  %mm 
    k = 2*pi/lambda;
    T = k/ff*(X.^2+Y.^2);
    lens = T;
    %zernike corrections
    pist = 0; x_tilt = 0; y_tilt = 0; h_coma = 0.0; v_coma = 0.0; spher = 0;
    
    HoloL = lens + zernike(pist, x_tilt , y_tilt , xx.Zernikamplitudes_foc , xx.Zernikamplitudes_astig1,  0.08 , h_coma , v_coma , spher, H,V,RHO,PHI);
    
    %Calibration
    poly_file = fopen('polyfit_par_diode.txt','r');
    formatSpec = '%f';
    header = fgetl(poly_file);
    polyvalues_cali = str2num(fgetl(poly_file));
    fclose(poly_file);
    GREYVAL = polyval(polyvalues_cali, mod(HoloL,2*pi)); 
    
    % Grascale normalization from [0, 2Pi]to [0 255]
    GREYVAL = GREYVAL./255; %Normalise holograms
    HoloSLM = [GREYVAL];
    % display onto SLM
    device_number=2;
    fullscreen(HoloSLM, device_number)

    %% CAMERA part
        
        % Set the FIFO frame buffer size. Default size is 1.
        cam.MaximumNumberOfFramesToQueue = 1;
    
        % Start software triggered image acquisition
        % disp('Starting software triggered image acquisition.');
    
        % Set the number of frames per software trigger and start trigger
        % acquisition
        cam.FramesPerTrigger_zeroForUnlimited = 1;
        cam.Arm;  
    
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
            imageHeight = imageFrame.ImageData.Height_pixels;
            imageWidth = imageFrame.ImageData.Width_pixels;
            imageData2D = reshape(uint16(imageData), [imageWidth, imageHeight]);
            imageData2D(858:876, 398:429) = zeros(19,32);
            fig3 = figure(3);
            fig3; imagesc(imageData2D'), colormap(gray), colorbar 
            % disp("ok")
        end
    
        %find correct value for exposure
        exposure = cam.ExposureTime_us;
        max_imageData2D = max(max(imageData2D,[],"all"));
        maxPixelValue = double(2^cam.BitDepth - 3);
        if (max_imageData2D >= maxPixelValue)
            % Release the camera
            disp('Releasing the camera');
            cam.Dispose;
            delete(cam);
            % Release the serial numbers
            delete(serialNumbers);
            % Release the TLCameraSDK.
            tlCameraSDK.Dispose;
            delete(tlCameraSDK);
            error("Camera saturated")
        end
        max_count = maxPixelValue-max_imageData2D;
        % Stop software triggered image acquisition
        % disp('Stopping software triggered image acquisition.');
        cam.Disarm;
    
        
    end %if

%%

%zernike corrections function definition
function  Holo_correc  = zernike( pist, x_tilt , y_tilt ,  dfoc ,ver_ast, ob_ast, h_coma , v_coma , spher, H,V, RHO,PHI)

    Holo_0 = (zeros(V,H)+1)*pist;
    Holo_1 = x_tilt * (RHO.*cos(mod(PHI,2*pi))) ;
    Holo_2 = y_tilt* (RHO.*sin(mod(PHI,2*pi)) );
    Holo_3 = dfoc * ((2* ( RHO.* RHO) -1)) ;
    Holo_4 = ver_ast * (RHO.* RHO.*cos(2 * mod(PHI,2*pi)));
    Holo_5 = ob_ast * (( RHO.* RHO).*sin(2 * mod(PHI,2*pi)));
    Holo_6 = h_coma * ((3*( RHO.* RHO)-2) .* RHO .* cos(mod(PHI,2*pi)));
    Holo_7 = v_coma * ((3*( RHO.* RHO)-2) .* RHO .* sin(mod(PHI,2*pi)));
    Holo_8 = spher * (6*(RHO.*RHO.*RHO.*RHO) - 6*(RHO.*RHO) + 1);
    
    
    Holo_correc =  Holo_0 + Holo_1 + Holo_2+ Holo_3+ Holo_4+ Holo_5+ Holo_6+ Holo_7+ Holo_8;   
end 


    