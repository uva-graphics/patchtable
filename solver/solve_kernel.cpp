
/* Main program for the kernel solver, which approximates a target kernel with smaller/simpler kernels.
 
   TODO:
    - Priority list:
      - Mandatory:
         - Run applications (Gaussian, Gradient integration, boundary interpolation)
           through our G.A.
         - Compare ours using C backend against Convolution Pyramids
         
    ---------------------------------------------------------------------------------------
 
         - Validate running times against compiler
            - Could analyze them each independently
            - Try using different features for the 4 different directions of IIR
            - Try using 5 x 5 grid of how many filters call each other filter
         - Application: Gradient integration, Gaussian blur
         - Have Jason help write paper
         - Final prefiltering:
            - Search sizes 5, 7 for prefilter, evaluate error using 8 samples of antialiasing
         - Comparison against straw man simulated annealing / GA
         - Compiler -- Add truncate option - Puneet?
         
    ---------------------------------------------------------------------------------------
     Done:
         - Final G.A. algorithm:
            - Add -resume to solve_kernel solve mode, make run_ga.py script that runs the main G.A. and then
              reruns the top solutions using -resume.

     Optional:

      - Genetic algorithm / combinatoric search:
          - Straw man simulated annealing comparison.
          - Compiler: Multiply final program by appropriate 'scale' constant.
          - Antialiasing -- and validate that programs work
      - OpenTuner
          - Get working
          - Validate running time model, or improve it
      - Applications:
          - Gaussian blur
          - Poisson solver
      - Paper
      
    ---------------------------------------------------------------------------------------

    - Finished:
      - Test summed area table (IIR+x, IIR+y), Gaussian as multiple box blurs ((IIR+x, IIR+y)*3) (add to test suite).
      - Finish solver:
          - Termination rule, epsilon or hit max OF evaluations
          - Properly compute time if downsample/upsample involved

    ---------------------------------------------------------------------------------------
 
      - Test suite of classic reductions:
          - Gaussian => box filter, Gaussian => SVD separable, Gaussian => IIR, Gaussian => quadratic approximation with IIRs+FIRs.
          - Gabor filter => SVD separable
          - Box filters => summed area tables
          - Plus any new topologies we discover for the above
      - New results (Poisson integration, etc)
      - Solver:
          - Contract: Replace entire subpart with identity
          - Adaptive
          - For downsample, reduce any estimated time by 4.
          - Topology search (genetic mutations) -- add, remove, replace, SVD, integral via IIR
             - Allow parallel (place state for optimization in a class) for topology search / genetic algorithm
      - More kernels
      - Working with compiler
      - Applications
          - Title? "SparseConv: Accelerating Convolutions Using Sparsity" (or Approximating. But accelerating sounds more exciting).
 
    (1) Carry forward/backward a "current" truncated/expanded filter and also merge this with the best one that has been
        seen so far. (?)
    (2) Avoid re-solving for the same sparsity pattern in the same optimization run.
    (3) "Rank all the pixels", run contract/expand on the top k ranked ones, and update rank only every m iterations.
    (4) Option to start from random guess or use previous guess when doing coefficient solve
 
      - For topology search, memoize. Do not search over two topologies that are equivalent. Use iterative deepening
        from all current best solutions (including always the 0-node graph) to search all topologies of size <= n.
        List all topologies on paper. Perhaps could determine if two topologies are equivalent by repeatedly merging
        data structures such as set for associative and commutative operators such as +, *.
        
 
   Could do coarse-to-fine when solving for coefficients and sparsity pattern.
 */

#include <Halide.h>
using namespace Halide;
#include "../../Halide/apps/support/image_io.h"

#include "solver.h"
#include "simplify.h"
#include "mutate.h"
#include "test_solver.h"
#include <ctime>
#include <dirent.h>

/* -------------------------------------------------------------------------------------------
   Unit tests
   ------------------------------------------------------------------------------------------- */

template<class real>
void printCImg(CImg<real> I){
    printf("\n[");
    for (int y = 0; y < I.height() ; y++){
	    printf("[");
    	for (int x = 0; x < I.width() ; x++){
    		printf("%s",num_to_str(I(x,y,0,0)).c_str());
    		if (x != I.width()-1) { printf(", "); }
    	}
	    printf("]");
	    if (y != I.height()-1) { printf("\n "); }
    }
    printf("]\n\n");
}

template<class real>
void printHalideImage(Image<real> I){
    printf("\n[");
    for (int y = 0; y < I.height() ; y++){
	    printf("[");
    	for (int x = 0; x < I.width() ; x++){
    		printf("%s",num_to_str(I(x,y)).c_str()); // This just prints the red channel since the third param is assumed to be 0
    		if (x != I.width()-1) { printf(", "); }
    	}
	    printf("]");
	    if (y != I.height()-1) { printf("\n "); }
    }
    printf("]\n\n");
}

template<class real>
void printCImgWindow(CImg<real> image, int previewWindowUpperLeftXCoordinate, int previewWindowUpperLeftYCoordinate, int previewWindowHeight, int previewWindowWidth){
    printf("\n[");
    for (int y = previewWindowUpperLeftYCoordinate; y < previewWindowUpperLeftYCoordinate+previewWindowHeight ; y++){
	    printf("[");
    	for (int x = previewWindowUpperLeftXCoordinate; x < previewWindowUpperLeftXCoordinate+previewWindowWidth ; x++){
    		printf("%s",num_to_str(image(x,y,0,0)).c_str());
    		if (x != previewWindowUpperLeftXCoordinate+previewWindowWidth-1) { printf(", "); }
    	}
	    printf("]");
	    if (y != previewWindowUpperLeftYCoordinate+previewWindowHeight-1) { printf("\n "); }
    }
    printf("]\n\n");
}

template<class real>
void printArrayWindow(Array<real> a, int previewWindowUpperLeftXCoordinate, int previewWindowUpperLeftYCoordinate, int previewWindowHeight, int previewWindowWidth){
	printf("\n[");
	for (int y = previewWindowUpperLeftYCoordinate; y < previewWindowUpperLeftYCoordinate+previewWindowHeight ; y++){
	    printf("[");
    	for (int x = previewWindowUpperLeftXCoordinate; x < previewWindowUpperLeftXCoordinate+previewWindowWidth ; x++){
    		printf("%s",num_to_str(a(y,x)).c_str());
    		if (x != previewWindowUpperLeftXCoordinate+previewWindowWidth-1) { printf(", "); }
    	}
	    printf("]");
	    if (y != previewWindowUpperLeftYCoordinate+previewWindowHeight-1) { printf("\n "); }
    }
    printf("]\n\n");
}

template<class real>
void printHalideImageWindow(Image<real> image, int previewWindowUpperLeftXCoordinate, int previewWindowUpperLeftYCoordinate, int previewWindowHeight, int previewWindowWidth){
    printf("\n[");
    for (int y = previewWindowUpperLeftYCoordinate; y < previewWindowUpperLeftYCoordinate+previewWindowHeight ; y++){
	    printf("[");
    	for (int x = previewWindowUpperLeftXCoordinate; x < previewWindowUpperLeftXCoordinate+previewWindowWidth ; x++){
    		printf("%s",num_to_str(image(x,y)).c_str());
    		if (x != previewWindowUpperLeftXCoordinate+previewWindowWidth-1) { printf(", "); }
    	}
	    printf("]");
	    if (y != previewWindowUpperLeftYCoordinate+previewWindowHeight-1) { printf("\n "); }
    }
    printf("]\n\n");
}

void test_paul() {
	
	printf("\n\n***************Test Paul***************\n\n");
	
	printf("\n\n***************Testing CImg functionality of Array Class***************\n\n");
	
    auto image1Name = string("gradient.png");
    auto image1 = CImg<double>(image1Name.c_str());
    
    int previewWindowUpperLeftXCoordinate = 50;
    int previewWindowUpperLeftYCoordinate = 50;
    int previewWindowHeight = 10;
    int previewWindowWidth = 10;
	
    printf("\n\nTesting Array initialization with PNG CImg input\n\n");
    auto image1Array = Array<double>(image1);
    
    printf("Image1 Name: %s\n\n", image1Name.c_str());
    
    printf("image1.height(): %d\n", image1.height());
    printf("image1.width(): %d\n\n", image1.width());
    
    printf("image1(%d:%d,%d:%d):", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printCImgWindow(image1, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
    printf("image1Array.height(): %d\n", image1Array.height());
    printf("image1Array.width(): %d\n\n", image1Array.width());
    
    printf("image1Array(%d:%d,%d:%d):", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printArrayWindow(image1Array, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
    printf("\n\nTesting Array to CImg Functionality\n\n");
    
    auto image2 = image1Array.getCImg();
    printf("image2 will be obtained from image1Array.\n\n");
    printf("image2(%d:%d,%d:%d): \n[", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printCImgWindow(image2, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
	printf("\n\n***************Testing Halide::Image functionality of Array Class***************\n\n");
	
    printf("\n\nTesting Array initialization with PNG CImg input\n\n");
    
    Halide::Image<double> image3 = load<double>(image1Name);
    save(image3, "out.png");
    
    printf("image3 is the same as image1.\n\n");
    printf("image3.height(): %d\n", image3.height());
    printf("image3.width(): %d\n\n", image3.width());
    
    printf("image3(%d:%d,%d:%d):", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printHalideImageWindow(image3, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
    auto image3Array = Array<double>(image3);
    
    printf("image3Array.height(): %d\n", image3Array.height());
    printf("image3Array.width(): %d\n\n", image3Array.width());
    
    printf("image3Array(%d:%d,%d:%d):", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printArrayWindow(image3Array, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
    printf("\n\nTesting Array to CImg Functionality\n\n");
    
    auto image4 = image3Array.getHalideImage();
    printf("image4 will be obtained from image3Array.\n\n");
    printf("image4(%d:%d,%d:%d): \n[", previewWindowUpperLeftXCoordinate, previewWindowUpperLeftXCoordinate+previewWindowHeight-1, previewWindowUpperLeftYCoordinate, previewWindowUpperLeftYCoordinate+previewWindowWidth-1);
    printHalideImageWindow(image4, previewWindowUpperLeftXCoordinate, previewWindowUpperLeftYCoordinate, previewWindowHeight, previewWindowWidth);
    
	printf("\n\n***************Testing Filters***************\n\n");
	
    static vector<Array<double> > temp_images;
    
    auto I1 = Array<double>::random({4, 6});
    auto I2 = Array<double>::random({4, 6});
    auto I3 = Array<double>::random({4, 6});
    
    auto K = Array<double>::random({3, 3});
    auto F = Array<double>::random({3});
    
    auto out1 = I1;
	
    printf("\n\nTesting FilterAdd\n\n");
    
    auto fAdd = shared_ptr<FilterAdd<double> >(new FilterAdd<double>());

    fAdd->apply({&I1, &I2, &I3}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("I2:\n%s\n\n", I2.str().c_str());
    printf("I3:\n%s\n\n", I3.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterUpsample\n\n");
    
    auto fUpsample = shared_ptr<FilterUpsample<double> >(new FilterUpsample<double>());
	
    fUpsample->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterDownsample\n\n");
    
    auto fDownsample = shared_ptr<FilterDownsample<double> >(new FilterDownsample<double>());
	
    fDownsample->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    printf("\n\nTesting FilterFIR\n\n");
    
    auto fFIR = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(K));
	
    fFIR->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("K:\n%s\n\n", K.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterIIR\n\n");
    
    //F.clear(1);
    K.resize({3});
    auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIR<double>(K,F));
	
    fIIR->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("K:\n%s\n\n", K.str().c_str());
    printf("F:\n%s\n\n", F.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterTranspose\n\n");
    
    auto fTranspose = shared_ptr<FilterTranspose<double> >(new FilterTranspose<double>());
	
    fTranspose->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterFlipH\n\n");
    
    auto fFlipH= shared_ptr<FilterFlipH<double> >(new FilterFlipH<double>());
	
    fFlipH->apply({&I1}, out1, temp_images);
    
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
    printf("\n\nTesting FilterDAG\n\n");
    
	auto fTranspose_base = static_pointer_cast<Filter<double> >(fTranspose);
	auto fFlipH_base = static_pointer_cast<Filter<double> >(fFlipH);
	
	FilterDAG<double> fDAG(vector<shared_ptr<Filter<double> > >{fTranspose_base, fFlipH_base, fTranspose_base});
	
	fDAG.apply({&I1}, out1, temp_images);
	
    printf("Will Perform a Tranpose, FlipH, and Transpose, i.e. a flip accross the vertical axis.\n\n");
    printf("I1:\n%s\n\n", I1.str().c_str());
    printf("out1:\n%s\n\n", out1.str().c_str());
    
	printf("\n\n***************Time Testing Filters***************\n\n");
    
	auto fFIR_base = static_pointer_cast<Filter<double> >(fFIR);
	auto fIIR_base = static_pointer_cast<Filter<double> >(fIIR);
	auto fAdd_base = static_pointer_cast<Filter<double> >(fAdd);
	auto fUpsample_base = static_pointer_cast<Filter<double> >(fUpsample);
	auto fDownsample_base = static_pointer_cast<Filter<double> >(fDownsample);
	
	FilterDAG<double> fDAG2(vector<shared_ptr<Filter<double> > >{fFIR_base, fIIR_base, fAdd_base, fTranspose_base, fFlipH_base, fUpsample_base, fDownsample_base});
	
    printf("Will Perform a FilterDAG that contains a FilterFIR (same K as above), FilterIIR (same K and F as above), "
    	   "FilterAdd (will just rewrite input to output since only adding a single image), FilterTranspose, FilterFlipH, "
    	   "FilterUpsample, and FilterDownsample on ImageArray1 100 times per run. There will be 10 runs.\n\n");
    
    printf("image1Array.height(): %d\n", image1Array.height());
    printf("image1Array.width(): %d\n\n", image1Array.width());
    
	clock_t begin, end; 
	
	for(int i = 0; i<10; i++){
		begin = clock();
		for(int j = 0; j<100; j++){
			fDAG2.apply({&image1Array}, out1, temp_images);
		}
		end = clock();
		printf("run[%d]: %f seconds.\n", i, double(end - begin) / CLOCKS_PER_SEC);
	}
}

void compare_main(bool show) {
	
	std::srand( (int) time(NULL) ); // Seed Random
	
	auto imageName = string("before_image"); // Make sure we're not destroying any image files that already exist
	while ( fopen( (imageName+".png").c_str() , "r") != NULL ){
		imageName.push_back('1');
	}
	Array<double> testingImage = (show) ? Array<double>::random({8,10}) :Array<double>::random({800,1000});
    save(testingImage.getHalideImage(), (imageName+".png").c_str() ); // Create new random image
    
	static const vector<string> filterNames = {string("fir"),string("add"),
											   string("iirRight"),string("iirDown"),string("iirLeft"),string("iirUp")};
	static vector<string> resultDirectoryNames = vector<string>(filterNames);
    
    for (int i=0; i < (int) resultDirectoryNames.size(); i++){
	    resultDirectoryNames[i].insert(0,"z");
		while ( opendir( (resultDirectoryNames[i]).c_str() ) != NULL ){ // Make sure we're not destroying any directories that already exist
	    	resultDirectoryNames[i].insert(0,"z");
		}
		std::system( ("mkdir "+resultDirectoryNames[i]).c_str() );
		std::system( ("cd ../compiler && ./master.py compile-single ./samples/compare_"+filterNames[i]+".json ../solver/"+resultDirectoryNames[i]).c_str() );
		std::system( ("cd ../compiler && ./master.py run-single ../solver/"+resultDirectoryNames[i]+"/ ../solver/"+imageName+".png").c_str() );
		//std::cout << "cd ../compiler && ./master.py compile-single ./samples/compare_"+filterNames[i]+".json ../solver/"+resultDirectoryNames[i] + " && cd ../solver" << "\n\n";
		//std::cout << "cd ../compiler && ./master.py run-single ../solver/"+resultDirectoryNames[i]+"/ ../solver/"+imageName+".png" + " && cd ../solver" << "\n\n";
	}
	std::system( "clear && clear" );
	
	auto I_image = load<double>( imageName+".png" ); // I refers to original image
	auto const I_array = Array<double>(I_image);
	
	static vector<Array<double> > temp_images;
	
	Array<double> J_solver_array; // J refers to after image
	Array<double> J_compiler_array;
	
	auto processedImageName = string("after_image"); 
	while ( fopen( (processedImageName+".png").c_str() , "r") != NULL ){ // Make sure we're not destroying any image files that already exist
		processedImageName.push_back('1');
	}
	
	printf("\nComparing Compiler & Solver Filters.\n");
	printf("Original Image Height: %d\n", I_array.height());
	printf("Original Image Width: %d\n", I_array.width());
	printf("Original Image Number of Pixels: %d\n", I_array.nelems);
	
    for (int i=0; i < (int) filterNames.size(); i++){
		// load compiler output
		J_compiler_array = load<double>( string("./"+resultDirectoryNames[i]+"/output.png") );
		
		// load solver output
		printf("\n\nCurrent Filter: %s\n", filterNames[i].c_str());
		if ( filterNames[i] == string("add")) { // Add
			auto fAdd = shared_ptr<FilterAdd<double> >(new FilterAdd<double>()); 
			fAdd->apply({&I_array, &I_array}, J_solver_array, temp_images);
		} else if ( filterNames[i] == string("fir")) {
            // K:
            // [[ 1, 2, 3 ],
            //  [ 4, 5, 6 ],
            //  [ 7, 8, 9 ]]
            auto K = Array<double>::random({3, 3});
            K(0,0) = 1;
            K(0,1) = 2;
            K(0,2) = 3;
            K(1,0) = 4;
            K(1,1) = 5;
            K(1,2) = 6;
            K(2,0) = 7;
            K(2,1) = 8;
            K(2,2) = 9;
            
			auto fFIR = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(K));
			fFIR->apply({&I_array}, J_solver_array, temp_images);
		} else {// else we've got an iir
			
            auto K = Array<double>::random({5});
			K(0) =  1;
			K(1) =  2;
			K(2) =  3;
			K(3) =  4;
			K(4) =  5;
			
			auto F = Array<double>(vector<int>{5});
			F(0) =  1;
			F(1) =  2;
			F(2) =  3;
			F(3) =  4;
			F(4) =  5;
			
			if ( filterNames[i] == string("iirUp")) {
				auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIRMinusY<double>(K,F));
				fIIR->apply({&I_array}, J_solver_array, temp_images);
			} else if ( filterNames[i] == string("iirDown")) {
				auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIRPlusY<double>(K,F));
				fIIR->apply({&I_array}, J_solver_array, temp_images);
			} else if ( filterNames[i] == string("iirLeft")) {
				auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIRMinusX<double>(K,F));
				fIIR->apply({&I_array}, J_solver_array, temp_images);
			} else if ( filterNames[i] == string("iirRight")) {
				auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIR<double>(K,F));
				fIIR->apply({&I_array}, J_solver_array, temp_images);
			}
		}
		
		J_solver_array.normalize(0,1);
		save(J_solver_array.getHalideImage(), string("./"+processedImageName+".png").c_str() );
		J_solver_array = load<double>( string("./"+processedImageName+".png") );
		
		if ( J_solver_array.height() == J_compiler_array.height() && J_solver_array.width() == J_compiler_array.width() ){
			auto difference_array = J_solver_array;
			difference_array.clear();
			int numBadPixels = 0;
			double totalDifference = 0;
			for (int y = 0; y < J_solver_array.height() ; y++){
				for (int x = 0; x < J_solver_array.width() ; x++){
					if ( std::abs(J_solver_array(y,x)-J_compiler_array(y,x)) >= 0.01 ){ // Threshold to tell if pixel difference is noticeable
						difference_array(y,x) = std::abs(J_solver_array(y,x)-J_compiler_array(y,x));
						totalDifference += difference_array(y,x);
						numBadPixels++;
					}
				}
			}
			if (show){
				printf("\nOriginal Image:\n");
				printf(I_array.str().c_str());
				printf("\nSolver Output:\n");
				printf(J_solver_array.str().c_str());
				printf("\nCompiler Output:\n");
				printf(J_compiler_array.str().c_str());
				printf("\nDifference:\n");
				printf(difference_array.str().c_str());
				printf("\n");
			}
			printf("Number of noticably erroneous pixels: %d\n", numBadPixels);
			printf("Number of total pixels: %d\n", J_solver_array.nelems);
			printf("Total Difference: %f\n", totalDifference);
		} else {
			printf("Dimensions of compiler output and solver output don't match.\n");
			printf("Solver Output Height: %d\n", J_solver_array.height());
			printf("Solver Output Width: %d\n", J_solver_array.width());
			printf("Compiler Output Height: %d\n", J_compiler_array.height());
			printf("Compiler Output Height: %d\n", J_compiler_array.width());
		}
	}
	
	
	// Clean up
	//*
	std::system( ("rm -rf "+imageName+".png").c_str() );
	std::system( ("rm "+processedImageName+".png").c_str() );
	for (int i=0; i < (int) filterNames.size(); i++){
		std::system( ("rm -rf "+resultDirectoryNames[i]).c_str() );
    }
    //*/
}

/* -------------------------------------------------------------------------------------------
   Main program
   ------------------------------------------------------------------------------------------- */

void usage() {
    fprintf(stderr, "usage: solve_kernel test                                     -- Run unit tests\n"
                    "       solve_kernel test_solver                              -- Test solver\n"
                    "       solve_kernel test_simplify                            -- Test simplify\n"
                    "       solve_kernel test_mutate                              -- Test mutate\n"
                    "       solve_kernel test_multiscale                          -- Test multiscale\n"
                    "       solve_kernel sample [n=100] [min_size=2] [max_size=5] -- Sample n topologies from min-max size, write pareto_full\n"
                    "       solve_kernel solve target.txt [json.txt] [params]     -- Run solver (use -max_combinatoric, -topology if no json.txt)\n"
                    "       solve_kernel ga target.txt [params]                   -- Run ga\n"
                    "       solve_kernel recalc tgt.txt json.txt [params]         -- Parse pareto_full, recalc time, error. Use -retain_all 1 to\n"
                    "                                                                keep all points in the list and not just the frontier."
                    "       solve_kernel filter [filter_type] [input_image] [output_image]\n"
                    "       solve_kernel apply in.txt json.txt i out.txt [-time_apply 1] -- Apply point i of Pareto frontier in json.txt to in.txt\n"
                    "       solve_kernel compare [options]				          -- Compares filter output with those of compiler, add -s to show difference values\n"
                    "\n"
                    "The filter types are:\n"
                    "       fir:		Finite Impulse Response Filter.\n"
                    "       iir:		Infinite Impulse Response Filter.\n"
                    "       dag:		DAG Filter.\n"
                    "       add:		Adds 2+ images.\n"
                    "       transpose:	Image Transpose.\n"
                    "       fliph:		Horizontal Image Flip.\n"
                    "       upsample:	Image Upsample.\n"
                    "       downsample:	Image Downsample.\n"
                    "\n"
                    "Genetic algorithm parameters:\n"
                    "     -max_combinatoric n     -- Search all DAGs of size <= n\n"
                    "     -ga_downsample b        -- Permit downsample in G.A.\n"
                    "     -ga_always_downsample b -- Require every filter to contain downsample\n"
                    "     -ga_generations n       -- Number of G.A. generations\n"
                    "     -ga_population n\n"
                    "     -ga_frac_elitism f, -ga_frac_crossover f, -ga_frac_mutate f\n"
                    "     -ga_tournament_size n\n"
                    "\n"
                    "Solver parameters:\n"
                    "     -random n               -- Number of re-randomization outer loops\n"
                    "     -contract n             -- Number of contract iterations\n"
                    "     -expand n               -- Number of expand iterations\n"
                    "     -translate n            -- Number of translate iterations\n"
                    "     -iters n                -- Change all of contract, expand, translate iterations\n"
                    "     -filter_w w             -- Maximum size of FIR or IIR filter\n"
                    "     -target_w w             -- If positive crops or pads target to size w x w (default -1)\n"
                    "     -combinatoric n         -- Search combinatorically all topologies up to size n\n"
                    "     -topology name          -- Search only the given topology name (listed in test_solver.h)\n"
                    "     -pareto file.txt        -- Summary Pareto frontier file\n"
                    "     -pareto_full file.txt   -- Full Pareto frontier file\n"
                    "     -trace file.json        -- Write trace file\n"
                    "     -out prefix             -- Write Pareto/Pareto full/trace with given prefix\n"
                    "     -verbose n              -- Verbosity level (0, 1, 2, 3)\n"
                    "     -truncate_row_col b     -- Truncate rows and cols (0 or 1)\n"
                    "     -multires b             -- Multires solver (0 or 1)\n"
                    "     -start_level n          -- Start processing after n levels in multires\n"
                    "     -stop_level n           -- Stop after processing n levels in multires\n"
                    "     -feature_coeff file.txt -- Feature coeffs returned by fit.py, specific to architecture\n"
                    "     -solve_initial_maxiters n, -solve_initial_epsilon e, -solve_add_maxiters n, -solve_add_epsilon e\n"
                    "     -reoptim_add b          -- Run a second optimization when adding a program?\n"
                    "     -reoptim_thresh T\n"
                    "     -contract_is_random b\n"
                    "     -expand_is_random b\n"
                    "     -contract_recheck n\n"
                    "     -expand_recheck n\n"
                    "     -seed n\n"
                    "     -translate_step x\n"
                    "     -translate_all b, -translate_some b, -translate_duplicate b, -translate_duplicate_outside b\n"
                    "     -translate_cross b, -translate_copy b, -translate_border b\n"
                    "     -random_order b\n"
                    "     -visit_one b\n"
                    "     -scale_correct b\n"
                    "     -search_shifts b\n"
                    "     -shift_xmin n, -shift_ymin n, -shift_xmax n, -shift_ymax n\n"
                    "     -check_often b\n"
                    "     -uniform_error b\n"
                    "     -uniform_error_samples n\n"
                    "     -recalc_solve b         -- Run solver when doing 'recalc'\n"
                    "     -symmetry b             -- Enforce symmetry if present in target\n"
                    "     -symmetry_w_min w, -symmetry_w_max w\n"
                    "     -random_init b          -- Whether to seed with random instances\n"
                    "     -random_count n         -- How many random instances to seed with\n"
                    "     -max_of_calls n         -- Stop each topology optimization after n OF calls\n"
                    "     -adaptive b, -adaptive_iters n\n"
                    "     -try_separable b\n"
                    "     -merge_base b\n"
                    "     -multisize b, -multisize_sizes n\n"
                    "     -quantize b, -quantize_count n, -quantize_until iter\n"
                    "     -terminate_early b, -terminate_epsilon f, -terminate_iters n\n"
                    "     -mask_samples n\n"
                    "     -pareto_error_thresh E, -pareto_min_points n, -pareto_min_points_error E\n"
                    "     -antialias b            -- Whether to brute force antialias (default: true)\n"
                    "     -antialias_levels n\n"
                    "     -antialias_subsample n\n"
                    "     -antialias_interval b\n"
                    "     -in s                   -- Load input image (2D matrix text file)\n"
                    "     -preconvolve s          -- Convolve in and target with 2D matrix text file\n"
                    "     -prefilter_sigma x      -- Sigma for prefilter Gaussian\n"
                    "     -prefilter_size n       -- Size of prefilter Gaussian\n"
                    "     -flatten b              -- Whether to flatten DownsamplePrefilter to be FIR + FIR + Downsample\n"
                    "     -feature n              -- Use features n, n=0 (ordinary feature), n=1 (rich descriptor)\n"
                    "     -resume pareto_full.txt -- Resume solve from given Pareto frontier.\n"
                    "     -weights w.txt          -- Use 2D matrix for weights on target.\n"
                    "     -ignore_boundary n      -- Do not evaluate OF on n pixel border of target/out.\n"
                    "     -scale_1norm b          -- Scale 1-norm to 1 for output (default: 0, matches target L2 norm).\n"
                    "     -vh_mode b              -- Use additionally horizontal and vertical FilterFIR\n"
                    "     -seed_convpyr b         -- Seed with convolution pyramids topologies\n"
                    "     -max_time T             -- Stop solver after T minutes (negative for no limit)\n"
                    );

    exit(1);
}

int main(int argc, char *argv[]) {
    argc--;
    argv++;
    
    if (argc == 0) { usage(); }
	
	if (string(argv[0]) == string("test_solver")) { // Test solver
		test_solver();
    } else if (string(argv[0]) == string("test_multiscale")) { // Test solver
        test_multiscale();
    } else if (string(argv[0]) == string("solve")) { // Main solver
        main_solver<adouble>(argc-1, argv+1);
    } else if (string(argv[0]) == string("ga")) { // Main ga
        main_ga<adouble>(argc-1, argv+1);
    } else if (string(argv[0]) == string("apply")) { // Apply filter
        main_apply<adouble>(argc-1, argv+1);
    } else if (string(argv[0]) == string("recalc")) {
        main_recalc<adouble>(argc-1, argv+1);
	} else if (string(argv[0]) == string("test_simplify")) { // Test solver
		test_simplify<adouble>();
	} else if (string(argv[0]) == string("test_mutate")) { // Test mutate
		test_mutate<adouble>(argc-1, argv+1);
	} else if (string(argv[0]) == string("sample")) { // Sample topologies
		main_sample<adouble>(argc-1, argv+1);
	} else if (string(argv[0]) == string("test")) { // Run Unit Tests
		//test_filters();
		test_paul();
	} else if (string(argv[0]) == string("compare")) { // Compares filter outputs to those of compiler
		if (argc == 2) {
			if (string(argv[1]) == string("-s")){
				compare_main(true);
				return 0;
			}
		}
		compare_main(false);
	} else if (string(argv[0]) == string("filter")) { // Filter input image against given filter input
		
		if (argc != 4) { usage(); }
		
    	static vector<Array<double> > temp_images;
		auto filterType = string(argv[1]);
        //Halide::Image<double> inputImage = load<double>(argv[2]);
        auto inputImageArray = Array<double>( load<double>(argv[2]) );
		//auto inputImageArray = Array<double>( CImg<double>( argv[2] ) );
		auto outputImageArray = Array<double>();
		

		if (filterType == string("fir")) {
			
            // K:
            // [[ 1, 2, 3 ],
            //  [ 4, 5, 6 ],
            //  [ 7, 8, 9 ]]
            auto K = Array<double>::random({3, 3});
            K(0,0) = 1;
            K(0,1) = 2;
            K(0,2) = 3;
            K(1,0) = 4;
            K(1,1) = 5;
            K(1,2) = 6;
            K(2,0) = 7;
            K(2,1) = 8;
            K(2,2) = 9;
		    
			auto fFIR = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(K));
			fFIR->apply({&inputImageArray}, outputImageArray, temp_images);
			
			printf("fFIR->str().c_str():\n%s\n",fFIR->str().c_str());
			
		} else if (filterType == string("iir")) {
			
			auto K = Array<double>(vector<int>{1});
			K(0) = 1;
			
			auto F = Array<double>(vector<int>{3});
			F(0) =  0.5;
			F(1) =  0.3;
			F(2) =  0.2;
			
			auto fIIR = shared_ptr<FilterIIR<double> >(new FilterIIR<double>(K,F));
			fIIR->apply({&inputImageArray}, outputImageArray, temp_images);
			
			printf("fIIR->str().c_str():\n%s\n",fIIR->str().c_str());
			
		} else if (filterType == string("dag")) {
			// Transpose, FlipH, Transpose
			// Shoud just flip the image vertically

			auto fTranspose = shared_ptr<FilterTranspose<double> >(new FilterTranspose<double>());
			auto fTranspose_base = static_pointer_cast<Filter<double> >(fTranspose);
			
			auto fFlipH = shared_ptr<FilterFlipH<double> >(new FilterFlipH<double>());
			auto fFlipH_base = static_pointer_cast<Filter<double> >(fFlipH);
			
			FilterDAG<double> fDAG(vector<shared_ptr<Filter<double> > >{fTranspose_base, fFlipH_base, fTranspose_base});
		
			fDAG.apply({&inputImageArray}, outputImageArray, temp_images);
			
			fprintf(stdout, 
					"DAG Elements: \n"
					"    -Transpose\n"
					"    -FlipH\n"
					"    -Transpose\n"
					);
		} else if (filterType == string("add")) {
			// adds image to itself (i.e. doubles values)
			auto fAdd = shared_ptr<FilterAdd<double> >(new FilterAdd<double>());
			fAdd->apply({&inputImageArray,&inputImageArray}, outputImageArray, temp_images);
			printf("Added image to itself.\n");
		} else if (filterType == string("transpose")) {
			auto fTranspose = shared_ptr<FilterTranspose<double> >(new FilterTranspose<double>());
			fTranspose->apply({&inputImageArray}, outputImageArray, temp_images);
		} else if (filterType == string("fliph")) {
			auto fFlipH = shared_ptr<FilterFlipH<double> >(new FilterFlipH<double>());
			fFlipH->apply({&inputImageArray}, outputImageArray, temp_images);
		} else if (filterType == string("upsample")) {
			auto fUpsample = shared_ptr<FilterUpsample<double> >(new FilterUpsample<double>());
			fUpsample->apply({&inputImageArray}, outputImageArray, temp_images);
		} else if (filterType == string("downsample")) {
			auto fDownsample = shared_ptr<FilterDownsample<double> >(new FilterDownsample<double>());
			fDownsample->apply({&inputImageArray}, outputImageArray, temp_images);
		} else {
			usage();
		}
		
		auto outputImage = outputImageArray.getHalideImage();
        save(outputImage, argv[3] );
        
	} else {
		usage();
	}
    
    return 0;
}

