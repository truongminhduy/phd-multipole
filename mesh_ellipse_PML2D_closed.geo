Include "parameters_gmsh_getdp.dat";
SetFactory("OpenCASCADE");

// lc_diff = lambda_m*2e-2;
// lc_mil  = lambda_m*4e-2;
// lc_pml  = lambda_m*1e-1;

lc_diff = lambda_m*1e-2;
lc_mil  = lambda_m*2e-2;
// lc_pml  = lambda_m*5e-2;
lc_pml  = lambda_m*2e-2;

dom_x = R_pml_in*2;
dom_y = R_pml_in*2;
PML_size = R_pml_out - R_pml_in;
// PML geometry
Point(1)  = {-dom_x/2,-dom_y/2, 0};
Point(2)  = {-dom_x/2, dom_y/2, 0};
Point(3)  = { dom_x/2, dom_y/2, 0};
Point(4)  = { dom_x/2,-dom_y/2, 0};

// Point(5)  = {-dom_x/2-PML_size,-dom_y/2-PML_size, 0};
// Point(6)  = {-dom_x/2-PML_size, dom_y/2+PML_size, 0};
// Point(7)  = { dom_x/2+PML_size, dom_y/2+PML_size, 0};
// Point(8)  = { dom_x/2+PML_size,-dom_y/2-PML_size, 0};

// Point(9)   = {-dom_x/2,-dom_y/2-PML_size, 0};
// Point(10)  = {-dom_x/2, dom_y/2+PML_size, 0};
// Point(11)  = { dom_x/2, dom_y/2+PML_size, 0};
// Point(12)  = { dom_x/2,-dom_y/2-PML_size, 0};

// Point(13)  = {-dom_x/2-PML_size,-dom_y/2, 0};
// Point(14)  = {-dom_x/2-PML_size, dom_y/2, 0};
// Point(15)  = { dom_x/2+PML_size, dom_y/2, 0};
// Point(16)  = { dom_x/2+PML_size,-dom_y/2, 0};

// Line(1) = {5, 9};
// Line(2) = {9, 12};
// Line(3) = {12, 8};
// Line(4) = {13, 1};
Line(5) = {1, 4};
// Line(6) = {4, 16};
// Line(7) = {14, 2};
Line(8) = {2, 3};
// Line(9) = {3, 15};
// Line(10) = {6, 10};
// Line(11) = {10, 11};
// Line(12) = {11, 7};
// Line(13) = {5, 13};
// Line(14) = {13, 14};
// Line(15) = {14, 6};
// Line(16) = {9, 1};
Line(17) = {1, 2};
// Line(18) = {2, 10};
// Line(19) = {12, 4};
Line(20) = {4, 3};
// Line(21) = {3, 11};
// Line(22) = {8, 16};
// Line(23) = {16, 15};
// Line(24) = {15, 7};

// Line Loop(1) = {13, 4, -16, -1};
// Plane Surface(1) = {-1};
// Line Loop(2) = {16, 5, -19, -2};
// Plane Surface(2) = {-2};
// Line Loop(3) = {19, 6, -22, -3};
// Plane Surface(3) = {-3};
// Line Loop(4) = {14, 7, -17, -4};
// Plane Surface(4) = {-4};
// Line Loop(5) = {15, 10, -18, -7};
// Plane Surface(5) = {-5};
// Line Loop(6) = {18, 11, -21, -8};
// Plane Surface(6) = {-6};
// Line Loop(7) = {12, -24, -9, 21};
// Plane Surface(7) = {-7};
// Line Loop(8) = {20, 9, -23, -6};
// Plane Surface(8) = {-8};



// ellipse
Ellipse(30) = {0,0,0,r_ellipse_1,r_ellipse_2};
Line Loop(31) = {30};
Plane Surface(10) = {31};

//source
Point(50)  = {xS,yS,0.}; 
// Circle(51) = {xS,yS,0.,Abs[xS]*2e-2}; 
Circle(51) = {xS,yS,0.,r_source}; 
Line Loop(52) = {51};
Plane Surface(11) = {52}; 
Point{50} In Surface{11};

// background
Line Loop(12) = {17, 8, -20, -5};
Plane Surface(9) = {12,31,52}; 

// detector
Point(41)  = {xD1,yD1,0.}; 
// Point(42)  = {xD2,yD2,0.}; 

Physical Point(100) = {4};     	// print   
Physical Point(101) = {50};		// source
Physical Point(102) = {41};		// detector1
Physical Point(103) = {42};		// detector1
Physical Line(200) = {5, 8, 17, 20};		    			// continuité pml
// Physical Line(201) = {1,2,3,10,11,12,13,14,15,22,23,24};	// dirichlet pml

// Physical Surface(301) = {1,3,5,7}; 	// PML corner
// Physical Surface(302) = {4,8}; 		// PML left right
// Physical Surface(303) = {2,6}; 		// PML top bot
Physical Surface(300) = {9}; 		// background
Physical Surface(304) = {10};       // object
Physical Surface(305) = {11};       // source

Mesh.CharacteristicLengthMax = lc_pml;
// Characteristic Length{ PointsOf{ Surface{:}; } } = lc_pml;
Characteristic Length{ PointsOf{ Surface{10}; } } = lc_diff;
// Characteristic Length{ PointsOf{ Surface{11}; } } = lc_diff*3e-1;
Characteristic Length{ PointsOf{ Surface{11}; } } = lc_diff;
