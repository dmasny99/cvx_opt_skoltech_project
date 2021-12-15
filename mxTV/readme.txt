  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHORS:

       JOACHIM DAHL
       ANYBODY TECHNOLOGY A/S, DENMARK
       E-MAIL: dahl.joachim@gmail.com
       
       PER CHRISTIAN HANSEN
       TECHNICAL UNIVERSITY OF DENMARK
       E-MAIL: pch@imm.dti.dk

       SOREN HOLDT JENSEN, TOBIAS LINDSTROM JENSEN
       AALBORG UNIVERSITY, DENMARK
       E-MAIL: shj@es.aau.dk, tlj@es.aau.dk

   REFERENCE:

    -  ALGORITHMS AND SOFTWARE FOR TOTAL VARIATION IMAGE RECONSTRUCTION VIA
       FIRST-ORDER METHODS
       NUMERICAL ALGORITHMS, 53 (2010), PP. 67-92

   SOFTWARE REVISION:

       V 1.0 - APRIL 2009

   SOFTWARE LANGUAGE:

       C and MATLAB 7.5 (or later)


****************************************************************
*                                                              *
*                  mxTV - a Matlab package for                 *
*             Total Variation Image Reconstruction             *
*                                                              *
*             Requires Matlab version 7.5 or later             *
*                                                              *
****************************************************************

This package includes Matlab and C codes for Total Variation (TV)
image reconstruction: denoising, inpainting, and deblurring.

If you use this package, please give reference to:

   J. Dahl, P.C. Hansen, S.H. Jensen and T.L. Jensen, Algorithms
   and Software for Total Variation Image Reconstruction via
   First-Order Methods, Numerical Algorithms, 53 (2010), pp 67-92.

The code is part of the project CSI: Computational Science in
Imaging, supported by the Danish Research Council for Technology
and Production Sciences.  The work was carried out at Aalborg
University and Technical University of Denmark.

Note that the FFTW software is needed only for TV deblurring.
If the FFTW code is not installed (as described below) the codes
for denoising and in-painting will still work.


Installation guide for Windows
------------------------------

1. Go to http://www.fftw.org/install/windows.html and get the pre-
   compiled FFTW for Windows.  To determine if your computer uses
   32 or 64 bits, execute the Matlab command mexext, which returns
   a string.  If the string includes the characters '64' then your
   computer uses 64 bits, otherwise it uses 32 bits.

2. Unzip the downloaded file and put it in a new folder, for example:
      c:\Program Files\fftw-zzz-dll
   where zzz denotes the version of the files.

3. Go to Start -> Control Panel -> System, select Advanced and then
   Environmental Variables.  In the bottom half of the window, under
   System variables, select Path, and append the full path to the
   directory where you unzipped the FFTW.  You do this by adding the
   path after the last ";" and finish with a ";".

4. Go to the directory where you keep your other Matlab toolboxes,
   and unzip the na28.zip files to a new folder mxTV.

5. Start Matlab and go to the above mxTV folder.  If Matlab is already
   running at this point, then restart Matlab to make the change in
   the path environment variable (step 3) effective in Matlab

6. Run "install_windows". If asked, select the Lcc-win32 compiler
   bundled with Matlab as compiler. If you have used another compiler
   previous, run "mbuild -setup" and select the Lcc-win32 compiler.

7. Add mxTV to Matlab's path: go to File -> Set Path -> Add Folder
   and choose the folder where mxTV is located.  Then save and close.
   Alternatively, you can use the addpath command in Matlab.

8. To learn more, try the three demos:
      TVdenoise_demo, TVinpaint_demo, and TVdeblur_demo.

This installation is tested with the Lcc compiler bundled with Matlab.
More options are available in install_windows.m.


Installation guide for Linux and Unix
-------------------------------------

1. Get the fftw3 software from www.fftw.org/download.html or via your
   package manager, and install it.  Two packages are needed: libfftw3-3
   (base software) and libfftw3-dev (development software, needed for
   header files).

2. Go to the directory where you keep your other Matlab toolboxes, and
   unzip the na28.zip files to a new directory mxTV.

3. Start Matlab and go to the above mxTV directory.

4. Run "install_linux".  Our experience is that when installing mxTV,
   you can ignore any warnings deriving from an officially unsupported
   version of gcc.

5. Add mxTV to Matlab's path: go to File -> Set Path -> Add Folder and
   choose the folder where mxTV is located. Then save and close.
   Alternatively, you can use the addpath command in Matlab.  You may
   have to restart Matlab at this point.

6. To learn more, try the three demos:
      TVdenoise_demo, TVinpaint_demo, and TVdeblur_demo.

If you have any problems on Windows, Linus or Unix, please check the
files install_windows.m and install_linux.m.
