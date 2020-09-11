Resolution {
  { Name Projection;
    System{{ Name M1; NameOfFormulation Ez; Type ComplexValue; }}
    Operation{
      GenerateSeparate[M1];
      EigenSolve[M1, 52,1.20e+01,5.50e+01,EigFilter[],
          { {-1}, {-1,0,0}, { -1.000e+00, 1.764e+02, -6.672e+04, 3.512e+06, -6.735e+08, 1.831e+10, -2.328e+12, 2.804e+13, -2.600e+15, 0.000e+00, 0.000e+00}  } ,
          { { 1}, { 1}    , { 1.000e+00, -4.923e+01, 1.651e+04, -5.523e+05, 9.652e+07, -2.009e+09, 2.387e+11, -2.368e+12, 2.115e+14}  } 
          ];
      SaveSolutions[M1];
    }
  }
}
