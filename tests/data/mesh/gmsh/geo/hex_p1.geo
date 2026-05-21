SetFactory("OpenCASCADE");

Box(1) = {0,0,0,1,1,1};

// boundary groups
Physical Surface("left")   = {1};
Physical Surface("right")  = {2};

Physical Surface("front")  = {3};
Physical Surface("back")   = {4};

Physical Surface("bottom") = {5};
Physical Surface("top")    = {6};

// domain
Physical Volume("domain") = {1};

// structured 2x2x2 mesh
Transfinite Curve "*" = 3;

Transfinite Surface "*";
Recombine Surface "*";

Transfinite Volume "*";

Mesh.RecombineAll = 1;
Mesh.Algorithm3D = 1;
