# PatchTable

This is the source code for the project [PatchTable](http://www.connellybarnes.com/work/project_pages/patchtable/), presented at ACM SIGGRAPH 2015. PatchTable speeds up patch queries by making a precomputed index data structure. The precomputation is assumed to be an offline process that can take seconds to minutes.

Currently this code runs on Mac and Windows. 

License
-------

The code can be used for non-commercial research purposes. For a commercial license, contact the authors.

Citations
---------

If you use this in an academic setting, please cite our paper ([BibTeX](http://www.connellybarnes.com/work/bib/2015_patchtable.bib)):

 * *PatchTable: Efficient Patch Queries for Large Datasets and Applications*. Connelly Barnes, Fang-Lue Zhang, Liming Lou, Xian Wu, Shi-Min Hu. ACM Transactions on Graphics (Proc. SIGGRAPH) 2015.

Install Dependencies
--------------------

The dependency libraries should be installed:

 * libpng
 * png++
 * OpenCV
 * FLANN
 * ANN
 * Boost

Build Command-line
------------------

PatchTable can be built on Mac OS 10.10:

    cd patchtable/proj
    make -j                           % Parallel build

PatchTable can be built on Windows 7 64-bit Operating System with Visual Studio 2013:

 * The user needs to change project property->Configuration Properties->VC++ Directories to point to their own directory and add the dependency libraries.
 * The `patchtable.exe` binary will be created in a directory such as `patchtable\patchtable\VSproj\ConsoleApplication1\x64\Release\`.


Command-line Usage
------------------

For the documentation of the command-line tool, run:

    ./patchtable

Example of matching with accuracy roughly similar to PatchMatch:
    
    ./patchtable match vidpair0/a.png vidpair0/b.png out.pfm

Example of matching with increased coherence (mimics smooth NNFs of PatchMatch):

    ./patchtable match vidpair0/a.png vidpair0/b.png out.pfm -speed 3 -coherence_spatial 4

Simple C++ Code Example
-----------------------

Here is a simple code fragment of how PatchTable could be used in your own C++ program:

    PatchTableParams p(argc, argv);  /* Initialize with command-line arugments */
    
    Array<float> a(load_color_image<float>(a_filename));
    Array<float> b(load_color_image<float>(b_filename));

    double T_build_table_start = wall_time();
    PatchTable<> table(&p, b);
    printf("Build table time: %f seconds\n", wall_time() - T_build_table_start);

    Array<double> ann;

    double T_query = table.lookup(a, ann);
    printf("Query time: %f seconds\n", T_query);

    /* Explanation of ann (nearest neighbor field matching from a -> b):
       int b_patch_x       = ann(a_patch_y, a_patch_x, NNF_X);
       int b_patch_y       = ann(a_patch_y, a_patch_x, NNF_Y);
       double b_patch_dist = ann(a_patch_y, a_patch_x, NNF_DIST); */

    save_color_image(ann, out_filename);

For more detailed usage of PatchTable, see the `test_match` function in `patchtable_main.cpp`.

GUI (Windows)
-------------

The Windows GUI code for image editing and stitching is in `patchtable_apps/ImgCollection/Superes`.

The Light field super-resolution code is in `patchtable_apps/HybridSuperres/Superes`. The input file lists are `HybridSuperres/Superes/Highresolution.txt` and `lowresolutionFile.txt`. These contain the locations for the high and low resolution PNG files for each of the camera views.
When the program runs sucessfully, it will create a folder named result in `HybridSuperres/Superes/`
All the parameters are set in the code which are used to produce the same reference images in our paper. 
