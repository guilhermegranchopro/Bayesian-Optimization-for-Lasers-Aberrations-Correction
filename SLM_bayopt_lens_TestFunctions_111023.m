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
Zern_amplitudes_1 = optimizableVariable("Zernikamplitudes_foc",[-0.8,-0.4],"Type","real");
Zern_amplitudes_2 = optimizableVariable("Zernikamplitudes_astig1",[0.1,0.3],"Type","real");
Zern_amplitudes_3 = optimizableVariable("Zernikamplitudes_astig2",[-0.05,0.2],"Type","real");
Zern_amplitudes_4 = optimizableVariable("Zernikamplitudes_comma1",[-0.02,0.0],"Type","real");
Zern_amplitudes_5 = optimizableVariable("Zernikamplitudes_comma2",[-0.02,0.0],"Type","real");


foc_max_count_h = @(x)foc_peak_intensity(x,cam);
results = bayesopt(foc_max_count_h,Zern_amplitudes_1,'MaxObjectiveEvaluations',30);
opt_foc = results.XAtMinObjective.Zernikamplitudes_foc;

astig_max_count_h = @(x)astig_peak_intensity(x,cam,opt_foc);
results1 = bayesopt(astig_max_count_h,[Zern_amplitudes_2,Zern_amplitudes_3],'MaxObjectiveEvaluations',30);
opt_astig1 = results1.XAtMinObjective.Zernikamplitudes_astig1;
opt_astig2 = results1.XAtMinObjective.Zernikamplitudes_astig2;

comma_max_count_h = @(x)comma_peak_intensity(x,cam,opt_foc,opt_astig1,opt_astig2);
results2 = bayesopt(comma_max_count_h,[Zern_amplitudes_4,Zern_amplitudes_5],'MaxObjectiveEvaluations',30);
opt_comma1 = results2.XAtMinObjective.Zernikamplitudes_comma1;
opt_comma2 = results2.XAtMinObjective.Zernikamplitudes_comma2;

%% 5D Opt
Zern_amp_opt_1 = optimizableVariable("OPTamplitudes_foc",[opt_foc-(abs(opt_foc)*0.1),opt_foc+(abs(opt_foc)*0.1)],"Type","real");
Zern_amp_opt_2 = optimizableVariable("OPTamplitudes_astig1",[opt_astig1-(abs(opt_astig1)*0.1),opt_astig1+(abs(opt_astig1)*0.1)],"Type","real");
Zern_amp_opt_3 = optimizableVariable("OPTamplitudes_astig2",[opt_astig2-(abs(opt_astig2)*0.1),opt_astig2+(abs(opt_astig2)*0.1)],"Type","real");
Zern_amp_opt_4 = optimizableVariable("OPTamplitudes_comma1",[opt_comma1-(abs(opt_comma1)*0.1),opt_comma1+(abs(opt_comma1)*0.1)],"Type","real");
Zern_amp_opt_5 = optimizableVariable("OPTamplitudes_comma2",[opt_comma2-(abs(opt_comma2)*0.1),opt_comma2+(abs(opt_comma2)*0.1)],"Type","real");


final_max_count_h = @(x)final_peak_intensity(x,cam);
final_results = bayesopt(final_max_count_h,[Zern_amp_opt_1,Zern_amp_opt_2,Zern_amp_opt_3,Zern_amp_opt_4,Zern_amp_opt_5],'MaxObjectiveEvaluations',30);



%% Close camera
disp('Releasing the camera');
cam.Dispose;
delete(cam);
% Release the serial numbers
delete(serialNumbers);
% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);


%% Objective function for focus 
function foc_max_count = foc_peak_intensity(xx,cam)%TAKE PICTURE
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
    pist = 0; x_tilt = 0; y_tilt = 0; v_astig=0.25 ; o_astig=0.05  ; h_comma = 0.0; v_comma = 0.0; spher = 0;
    
    HoloL = lens + zernike(pist, x_tilt , y_tilt , xx.Zernikamplitudes_foc , v_astig,  o_astig , h_comma ,v_comma , spher, H,V,RHO,PHI);
    
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
            if 1
                fig3 = figure(3);
                fig3; imagesc(imageData2D'), colormap(gray), colorbar 
            end
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
        foc_max_count = maxPixelValue-max_imageData2D;
        % Stop software triggered image acquisition
        % disp('Stopping software triggered image acquisition.');
        cam.Disarm;
    
        
    end %if

%% Objective function for astigmatism 
function astig_max_count = astig_peak_intensity(xx,cam,opt_foc)%TAKE PICTURE
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
    pist = 0; x_tilt = 0; y_tilt = 0; foc = opt_foc; h_comma = 0.0; v_comma = 0.0; spher = 0;
    
    HoloL = lens + zernike(pist, x_tilt , y_tilt , foc , xx.Zernikamplitudes_astig1,  xx.Zernikamplitudes_astig2 ,h_comma , v_comma , spher, H,V,RHO,PHI);
    
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
            if 1
                fig3 = figure(3);
                fig3; imagesc(imageData2D'), colormap(gray), colorbar 
            end
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
        astig_max_count = maxPixelValue-max_imageData2D;
        % Stop software triggered image acquisition
        % disp('Stopping software triggered image acquisition.');
        cam.Disarm;
    
        
    end %if

%% Objective function for comma 
function comma_max_count = comma_peak_intensity(xx,cam,opt_foc,opt_astig1,opt_astig2)%TAKE PICTURE
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
    pist = 0; x_tilt = 0; y_tilt = 0; foc = opt_foc; v_astig = opt_astig1; o_astig = opt_astig2; h_comma = 0.0; v_comma = 0.0; spher = 0;
    
    HoloL = lens + zernike(pist, x_tilt , y_tilt , foc , v_astig,  o_astig , xx.Zernikamplitudes_comma1 , xx.Zernikamplitudes_comma2 , spher, H,V,RHO,PHI);
    
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
            if 1
                fig3 = figure(3);
                fig3; imagesc(imageData2D'), colormap(gray), colorbar 
            end
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
        comma_max_count = maxPixelValue-max_imageData2D;
        % Stop software triggered image acquisition
        % disp('Stopping software triggered image acquisition.');
        cam.Disarm;
    
        
    end %if


%% Objective function for 5d 
function final_max_count = final_peak_intensity(xx,cam)%TAKE PICTURE
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
    pist = 0; x_tilt = 0; y_tilt = 0; spher = 0;
    
    HoloL = lens + zernike(pist, x_tilt , y_tilt , xx.OPTamplitudes_foc , xx.OPTamplitudes_astig1,  xx.OPTamplitudes_astig2 , xx.OPTamplitudes_comma1 , xx.OPTamplitudes_comma2 , spher, H,V,RHO,PHI);
    
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
            if 1
                fig3 = figure(3);
                fig3; imagesc(imageData2D'), colormap(gray), colorbar 
            end
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
        final_max_count = maxPixelValue-max_imageData2D;
        % Stop software triggered image acquisition
        % disp('Stopping software triggered image acquisition.');
        cam.Disarm;
    
        
    end %if    
%%

%% Zernike corrections function definition
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


    