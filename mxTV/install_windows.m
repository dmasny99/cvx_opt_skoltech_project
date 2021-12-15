%INSTALL_WINDOWS Script used to install mxTV on the Windows platform.
%
% Compiles and links the mex files for the mxTV package.
%
% Tested with the Lcc compiler bundled with Matlab.
%
% See readme.txt for further instructions.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

%
% If you want to be able to break the execution of the programs, try to
% set CTRLCBREAK = 1, which uses a non-documented MATLAB API.
% If you do not have libut.lib in your Matlab dir, try CTRLCBREAK = 2.
% Default.
%
CTRLCBREAK=0;

if CTRLCBREAK==0
    sbreak = '';
elseif CTRLCBREAK==1
    sbreak = ['-DLIBUT -L' matlabroot '\extern\lib\win32\lcc -llibut'];
elseif CTRLCBREAK==2
    sbreak = ['-DLIBUT -Lexternlib -llibut'];
else
    error('Not a valid option for CTRLCBREAK')
end

any_error = false;

% Denosing.
try
    cs = sprintf('mex %s c/tv_denoise.c c/tools.c c/tv_core.c',sbreak);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdenoise because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

% Inpainting.
try
    cs = sprintf('mex %s c/tv_inpaint.c c/tools.c c/tv_core.c',sbreak);
    eval(cs)
catch
    any_error= true;    
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVinpaint because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end


% Deblurring.
deblur_error= false;
try
    cs = sprintf('mex -Iexternlib c/mxtrp.c c/tools.c externlib/libfftw3-3.lib');
    eval(cs)	

    cs = sprintf('mex -DFFTW3 -Iexternlib %s c/tv_deblur.c c/tv_core.c c/tools.c externlib/libfftw3-3.lib',sbreak);
    eval(cs)
  
    cs = sprintf('mex -DFFTW3 -Iexternlib %s c/tv_deblur_rr.c c/tv_core.c c/tools.c externlib/libfftw3-3.lib',sbreak);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdeblur because the compilation failed,')
    disp('probably because the fftw3.h header file is missing. To locate the') 
    disp('problem, follow the above instructions or see the readme.txt file.')
    disp('Ignore this error if you do not need TVdeblur.')
    deblur_error = true;
end

%We will just end by trying out if fftw has the correct setup
if deblur_error == false
	try	
        xtemp=TVdeblur(1,1,0.1);
        clear xtemp
    catch
        any_error = true;
        disp('----------------------------------------------------------------')
        disp('You will not be able to use TVdeblur because the setup of fftw3')
        disp('is wrong or fftw3 is missing. Follow the instructions in readme.txt')
        disp('to solve the problem if you would like to use TVdeblur.')
        disp('Ignore this error if you do not need TVdeblur.')
	end
end

if any_error == false && deblur_error == false
    disp('Install completed successfully.')
else
    disp('Installation did NOT complete successfully.')
end