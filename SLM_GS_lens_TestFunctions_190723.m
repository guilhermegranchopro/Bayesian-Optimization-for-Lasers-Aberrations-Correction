clc
clear all
close all

% addpath()

%Every time this code is run, if we want to change the hologram on display in the
%SLM we set the variable runalg to 1 - to run the GS algorithm - and we need to change the varable GS_name 
%to not overwrite on the previous one; if we want to take a new photo we
%also need to change the name of the image - variables str and str2

%%
H=1280; V=1024; %%Number of Horizontal and Vertical pixels
PixelSize = 12.5e-3; % [mm] pixel
x = PixelSize.*(-H/2:(H/2-1)); %[mm] x width of window
y = PixelSize.*(-V/2:(V/2-1)); %[mm] y width of window
[X,Y] = meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);

Holo = zeros(V,H); Holo_Blank = zeros(V,H);
%Nx = 220; Ny = 0;%%Number of horizontal and vertical grooves
%Gx = Nx/(H*PixelSize); Gy = Ny/(V*PixelSize); %spatial frequencies in units of lines per screen

lambda=660e-6;

runalg = 1;

%% Define the target intesity distribution
phasetype = 3; % 0 = two dots, 1 = complex object

switch phasetype
    case 0 % Two dots

        p_hole = 0.25;
        x_hole = [p_hole,p_hole-0.1]; %[mm]
        y_hole = [-p_hole,p_hole]; %[mm]
        r_hole = 0.25; %[mm]
        w_hole = r_hole/10; %[mm] Beam waist of Gaussian
        A = zeros(V,H);
        for jhole = 1:length(x_hole)
            x_ap = x_hole(jhole);
            y_ap = y_hole(jhole);
            RHO = ((X-x_ap).^2 + (Y-y_ap).^2).^0.5;
            PHI = atan2((Y-y_ap),(X-x_ap));
            A(RHO <= r_hole) = exp(-(RHO(RHO <= r_hole)).^2/w_hole^2);
        end
        Azo = A;
        GS_name = 'holo_complex1_GS_seq1';

    case 1 % complex object
        % Define Fourier parameters
        xf = linspace(1,H,H); % extent of Fourier plane
        yf = linspace(1,V,V); % extent of Fourier plane
        sigma_q = 5; % sigma of gaussian peaks
        c = 3; % complexity parameter (number of gaussian peaks = 2c)
        r_q = 100; % radial position of the peaks
        
        % Maths
        xq = xf - length(x)/2;
        yq = yf - length(y)/2;
        [Xf,Yf] = meshgrid(xq,yq);
        M = c*exp(-(Xf.^2 + Yf.^2)/2/sigma_q^2);
        for i = 1:1:2*c
            xx = cos(2*i*pi./(2*c))*r_q;
            yy = sin(2*i*pi/(2*c))*r_q;  
            M = M + exp(-((Xf-xx).^2 + (Yf-yy).^2)/2/sigma_q^2);
        end
        Azo = abs(fftshift(fft2(M)));
        GS_name = ['holo_complex_object_c',num2str(c)];

        case 2 %letter M
            Azo = zeros(V,H);
            M = imread('small.jpg');
            M = ~ (M/max(max(M)));
            Azo(V/2:V/2+size(M,1)-1,H/2:H/2+size(M,2)-1) = double(M);
            GS_name = 'holo_M_GS_seq1';
            imagesc(Azo)
            %pause
        case 3 % lambda
            Azo = zeros(V,H);
            lambda = imread('lambda_small.png');
            lambda = ~(lambda/max(max(lambda)));
            Azo(V/2-size(lambda,1)/2:V/2+size(lambda,1)/2-1,H/2-size(lambda,2)/2+1:H/2+size(lambda,2)/2) = double(lambda);
            GS_name = 'holo_lambda_GS_seq1';
            imagesc(Azo)
            %pause

        case 4 % Two lines
            linewidth_x=2;
            linewidth_y=100;
            separation=40 ;
            shift =100; %horizontal
            Azo = zeros(V,H);
            Azo(V/2-linewidth_y/2:V/2+linewidth_y/2,H/2-linewidth_x/2-separation/2+shift:H/2+linewidth_x/2-separation/2+shift) = ones();
            Azo(V/2-linewidth_y/2:V/2+linewidth_y/2,H/2-linewidth_x/2+separation/2+shift:H/2+linewidth_x/2+separation/2+shift) = ones();
            imagesc(Azo)
            GS_name = 'holo_lines_GS_seq1';

        case 5
            linewidth = 10*PixelSize;
            alpha = 45/180*pi;
            Azo = zeros(V,H);
            RHO = ((X).^2 + (Y).^2).^0.5;
            R = 50*PixelSize;
            Azo((Y < (tan(alpha)*X + linewidth)) & (Y > (tan(alpha)*X - linewidth)) & RHO < R) = ones();
            Azo((Y < (-tan(alpha)*X + linewidth)) & (Y > (-tan(alpha)*X - linewidth)) & RHO < R ) = ones();

            imagesc(Azo)
            GS_name = 'holo_x_GS_seq1';
            

end
    
        

if(runalg==1)

    %% Run GS algorithm
    maxiter = 50;
    
    Holo1 = ones(V,H).*exp(1i*rand(V,H)*2*pi); %Uniform amplitude with random phase
    Hol = Holo1;
    % qval = 0.1;
    % Holo = ones(V,H).*exp(1i*qval*(X.^2 + Y.^2)); %Uniform amplitude with quadratic phase
    
    E_FFT = fftshift(fft2(Hol));
    E_FFT = abs(Azo).*exp(1i*angle(E_FFT)); %Reimpose target amplitude and keep propagated phase
    Eout = E_FFT;
         
    jiter = 0;
    % pause
    while jiter < maxiter %Iterative loop
        m = 0;
        Hol = ifft2(ifftshift(E_FFT));    
        Hol = ones(V,H).*exp(1i*angle(Hol));
        
        E_FFT = fft2(Hol);
        E_FFT = fftshift(E_FFT);
        Eout = E_FFT;
        E_FFT = abs(Azo).*exp(1i*angle(E_FFT));
        
        jiter = jiter+1;
        disp(num2str(jiter))
        Ioutnorm = abs(Eout).^2/sum(sum(abs(Eout).^2));
        Iazonorm = abs(Azo).^2/sum(sum(abs(Azo).^2));
    
        discrepancy = norm(Ioutnorm - Iazonorm)/norm(Ioutnorm);
        disp(['Target-actual intensity discrepancy = ' num2str(discrepancy)])
    end
    
    Ioutnorm = abs(Eout).^2/sum(sum(abs(Eout).^2));
    Iazonorm = abs(Azo).^2/sum(sum(abs(Azo).^2));
    
    discrepancy = norm(Ioutnorm - Iazonorm)/norm(Ioutnorm);
    disp(['Target-actual intensity discrepancy = ' num2str(discrepancy)])
    
    Holo = angle(Hol); % this is the output hologram of the GS algorithm for the SLM
    save(GS_name, "Holo")

else 
    Holo = load(GS_name);
end

%%

ff =- 5000;  %mm 
lamdba = 660e-6;
k = 2*pi/lamdba;
T = k/ff*(X.^2+Y.^2);
lens = T;

%%
if runalg
    HoloL = Holo+lens;
else
    HoloL = Holo.Holo+lens;
end
%zernike corrections
pist = 0; x_tilt = 0; y_tilt = 0; dfoc = -0.13;  ver_ast=0.15; ob_ast = 0.1; h_coma = 0.0002; v_coma = 0.01; spher = 0;
HoloL = HoloL + zernike(pist, x_tilt , y_tilt , dfoc ,ver_ast, ob_ast, h_coma , v_coma , spher, H,V, RHO,PHI);

%Calibration
poly_file = fopen('polyfit_par_diode.txt','r');
formatSpec = '%f';
type polyfit_par_diode.txt
header = fgetl(poly_file);
polyvalues_cali = str2num(fgetl(poly_file));


fclose(poly_file);
disp(polyvalues_cali);
GREYVAL = polyval(polyvalues_cali, mod(HoloL,2*pi)); 

% Grascale normalization from [0, 2Pi]to [0 255]
GREYVAL = GREYVAL./255; %Normalise holograms
HoloSLM = [GREYVAL];




fig=figure(1);
imagesc(HoloSLM)
colormap gray
clim([0,1]);
axis equal off;

device_number=2;
fullscreen(HoloSLM, device_number)

%%
  
%TAKE PICTURE

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

    % Start software triggered image acquisition
    disp('Starting software triggered image acquisition.');

    % Set the number of frames per software trigger and start trigger
    % acquisition
    cam.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
    cam.FramesPerTrigger_zeroForUnlimited = 1;
    cam.Arm;

    %GET CORRECT EXPOSURE TIME
    
    cam.ExposureTime_us = 10000;  %define a large exposure to begin
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

        imageHeight = imageFrame.ImageData.Height_pixels;
        imageWidth = imageFrame.ImageData.Width_pixels;
        imageData2D = reshape(uint16(imageData), [imageWidth, imageHeight]);
        imageData2D(858:876, 398:429) = zeros(19,32);
        fig2 = figure(2);
        fig2; imagesc(imageData2D'), colormap(gray), colorbar 
        disp("ok")
    end

    %find correct value for exposure
    exposure = cam.ExposureTime_us;
    max_imageData2D = max(max(imageData2D,[],"all"));
    while(max_imageData2D >= maxPixelValue)
        cam.ExposureTime_us = exposure/2;
        exposure = exposure/2;
        max_imageData2D = max(max(imageData2D,[],"all"));
        disp(max_imageData2D);
        disp(cam.ExposureTime_us);
        %take new frame
        cam.IssueSoftwareTrigger;
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
        end 
    end
    if (max(max(imageData2D,[],"all"))< 800)
        cam.ExposureTime_us = 1000;
        
    end
    disp('Exposure time in  microseconds: ');
    disp(cam.ExposureTime_us);
    
    %GET FRAME

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
        fig3 = figure(3);
        fig3; imagesc(imageData2D'), colormap(gray), colorbar 
        str = "28_09_23_seq1.png";
        str2 = "28_09_23_seq1.tif";
        saveas(fig3,str);
        imwrite(imageData2D,str2);

    end

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

%%

%zernike corrections function definition
function  Holo_correc  = zernike( pist, x_tilt , y_tilt , dfoc ,ver_ast, ob_ast, h_coma , v_coma , spher, H,V, RHO,PHI)

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

    