/////////////////////////////////////////////////////////////////////////////////////////////////
Group {
  Point_print  = Region[100];
  Point_source = Region[101];
  Point_detect = Region[102];
  backround_nosrc = Region[300];
  src = Region[305];
  dif = Region[304];
  mid = Region[{src,backround_nosrc}];
  Omega_nopml = Region[{dif,mid}];
  Omega_out  = Region[{mid}];
  Omega_in   = Region[{dif}];
  Omega      = Region[{Omega_in,Omega_out}];
  Omega_plot = Region[{Omega,-src}];
  SurfDirichlet = Region[200];
}
/////////////////////////////////////////////////////////////////////////////
Function{
  Freq   = cel/lambda0;
  omega0 = 2.*Pi*cel/lambda0;
  k0     = 2.*Pi/lambda0;   

  epsilon[mid]  = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon[dif]  = Complex[eps_diff_re,eps_diff_im] * TensorDiag[1,1,1];
  // epsilon1[mid] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  // epsilon1[dif] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  mur[mid] = TensorDiag[1,1,1];
  mur[dif] = TensorDiag[1,1,1];

  If (E_or_H == 0)
    chi[]  = CompZZ[epsilon[]];
    // chi1[] = CompZZ[epsilon1[]];
    xsi[]  = Transpose[mur[]]/(CompXX[mur[]]*CompYY[mur[]]-CompXY[mur[]]*Conj[CompYX[mur[]]]);  
  Else
    chi[]  = CompZZ[mur[]];
    xsi[]  = Transpose[epsilon[]]/(CompXX[epsilon[]]*CompYY[epsilon[]]-CompXY[epsilon[]]*Conj[CompYX[epsilon[]]]);  
  EndIf

  // k_mid  = k0*Sqrt[eps_mil_re];
  // alpha_mid = k_mid*Sin[theta];
  // beta_mid  = k_mid*Cos[theta];
  // u_i[pml] = 0.;
  // u_i[mid] = A*Complex[ Cos[alpha_mid*X[]-beta_mid*Y[]] , Sin[alpha_mid*X[]-beta_mid*Y[]] ];
  // u_i[dif] = A*Complex[ Cos[alpha_mid*X[]-beta_mid*Y[]] , Sin[alpha_mid*X[]-beta_mid*Y[]] ];
  // source[] = (omega0/cel)^2*(epsilon1[]-epsilon[])*u_i[];
  source[src]        = 1/(r_source^2*Pi);
  source[Omega_plot] = 0;
  sourcep[Point_source] = 1;
  // source1[] = CompZZ[source[]];
  EigFilter[] = (Norm[$EigenvalueReal] > 1e-20);
}

//// FUNCTION //////////////////////////////////////////////////////////////////
Constraint {
  {Name Dirichlet; Type Assign;
    Case { { Region SurfDirichlet ; Value 0.; }}
  }
}

Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
  // { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
  { Name JLin ; Case { { Region All ; Jacobian Lin ; } } }
}

Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1 ; }
          { GeoElement Line        ; NumberOfPoints  4 ; }
          { GeoElement Triangle    ; NumberOfPoints  6 ; }
        }
      }
    }
  }
}
FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node   ; Support Region[{Omega,Point_source}]; Entity NodesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[{Omega,Point_source}]; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
  
  { Name Hgrad_perp; Type Form1P;
    BasisFunction {
      { Name un;  NameOfCoef un;  Function BF_PerpendicularEdge_1N; Support Region[{Omega,Point_source}]; Entity NodesOf[Omega]; }
      { Name un2; NameOfCoef un2; Function BF_PerpendicularEdge_2E; Support Region[{Omega,Point_source}]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
}
