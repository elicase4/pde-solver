SetFactory("OpenCASCADE");

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// physical boundaries
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};

Physical Surface("domain") = {1};

// structured mesh
Transfinite Curve {1,3} = 3;
Transfinite Curve {2,4} = 3;

Transfinite Surface {1};
Recombine Surface {1};
