
PatchTable
==========

Accelerates patch queries by taking up extra precomputation time (precomputation not optimized yet).

Build:

 - Install dependencies: libpng, png++, OpenMP.
 - Use make to build on Mac/Linux.
 - On Windows you will have to make a project that does something similar to the Makefile.
 - Run ./patchtable and you should see the command-line usage printed out.

Test:

 - Open MATLAB and run compare_patchtable
 - To get the PatchMatch comparison you will have to mexify PatchMatch by running one of
   build_mac.m/build_windows.bat/build_unix.sh in subdirectory patchmatch-2.1.



