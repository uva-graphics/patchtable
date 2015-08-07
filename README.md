# PatchTable: Efficient Patch Queries for Large Datasets and Applications

This is the source code for the project [PatchTable|http://www.connellybarnes.com/work/project_pages/patchtable/], presented at ACM SIGGRAPH 2015. PatchTable speeds up patch queries by making a precomputed index data structure. The precomputation is assumed to be an offline process that can take seconds to minutes.

Currently this code runs on Mac. It also can build on Windows. We hope to release Windows project files in a few days.

License
-------

The code can be used for non-commercial research purposes. For a commerical license, contact the authors.

Install Dependencies
--------------------

The dependency libraries should be installed:

 * libpng
 * png++
 * OpenCV
 * FLANN
 * ANN

Build Command-line
------------------

PatchTable can be built on Mac OS 10.10:

    cd patchtable/proj
    make -j                           % Parallel build

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

For more detailed usage of PatchTable, see the test_match function in patchtable_main.cpp.

GUI (Windows)
-------------

The Windows GUI for image editing and stitching is in patchtable_apps/ImgCollection/Patch.

