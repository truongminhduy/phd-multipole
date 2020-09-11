//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
//////////////////////////////////////////////////////////////////////////////
Group {
  Point_print  = Region[100];
  Point_source = Region[101];
  Point_detect = Region[102];
  PML_corner      = Region[301];
  PML_left_right  = Region[302];
  PML_top_bot     = Region[303];
  backround_nosrc = Region[300];
  src = Region[305];
  dif = Region[304];
  mid = Region[{src,backround_nosrc}];
  pml = Region[{PML_corner,PML_left_right,PML_top_bot}];
  Omega_nopml = Region[{dif,mid}];
  Omega_out  = Region[{pml,mid}];
  Omega_in   = Region[{dif}];
  Omega      = Region[{Omega_in,Omega_out}];
  Omega_plot = Region[{Omega,-src}];
  SurfDirichlet = Region[201];
}
///////////////////////////////////////////////////
Function{
  a_pml = Cos[angle];
  b_pml = Sin[angle];
  sx[PML_corner]     = Complex[a_pml,b_pml];
  sx[PML_left_right] = Complex[a_pml,b_pml];
  sx[PML_top_bot]    = 1;
  sx[Omega_nopml]    = 1;
  sy[PML_corner]     = Complex[a_pml,b_pml];
  sy[PML_left_right] = 1;
  sy[PML_top_bot]    = Complex[a_pml,b_pml];
  sy[Omega_nopml]    = 1;
  sz = 1; 
  pmltens[] = TensorDiag[sz*sy[]/sx[],sx[]*sz/sy[],sx[]*sy[]/sz];

  epsilon1[mid] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon1[dif] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon1[pml] = eps_mil_re*pmltens[];
  mur[mid] = TensorDiag[1,1,1];
  mur[dif] = TensorDiag[1,1,1];
  mur[pml] = pmltens[];

  If (E_or_H == 0)
    chi[] = CompZZ[epsilon1[]];
    xsi[] = Transpose[mur[]]/(CompXX[mur[]]*CompYY[mur[]]-CompXY[mur[]]*Conj[CompYX[mur[]]]);
  Else
    chi[] = CompZZ[mur[]];
    xsi[] = Transpose[epsilon1[]]/(CompXX[epsilon1[]]*CompYY[epsilon1[]]-CompXY[epsilon1[]]*Conj[CompYX[epsilon1[]]]);
  EndIf

  EigFilter[] = (Norm[$EigenvalueReal] > 1e-6);
}

Constraint {
  {Name Dirichlet; Type Assign;
    Case { { Region SurfDirichlet ; Value 0.; }}
  }
}

Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
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

/////// DEFINE FUNCTIONSPACE AND WEAK EQUATION /////////////////////////////////
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
      { Name un;  NameOfCoef un;  Function BF_PerpendicularEdge_1N; Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name un2; NameOfCoef un2; Function BF_PerpendicularEdge_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
}

 
////////////////////////////////////////////////////////////////////////////////////////////////
Formulation {
  {Name Ez; Type FemEquation;
    Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad;}}
    Equation {
      If (E_or_H == 0)
        Galerkin { Eig[ cel^2*xsi[]* Dof{d u}, {d u} ]; Rational 1; In Omega    ; Jacobian JVol; Integration Int_1;}
        Galerkin { Eig[ chi[]      * Dof{u}  , {u}   ]; Rational 2; In Omega_out; Jacobian JVol; Integration Int_1;}
        Galerkin { Eig[ chi[]      * Dof{u}  , {u}   ]; Rational 3; In Omega_in ; Jacobian JVol; Integration Int_1;}
      Else
        Galerkin { Eig[ cel^2*xsi[]* Dof{d u}, {d u} ]; Rational 1; In Omega_out; Jacobian JVol; Integration Int_1;}
        Galerkin { Eig[ cel^2*xsi[]* Dof{d u}, {d u} ]; Rational 3; In Omega_in ; Jacobian JVol; Integration Int_1;} 
        Galerkin { Eig[ chi[]      * Dof{u}  , {u}   ]; Rational 2; In Omega    ; Jacobian JVol; Integration Int_1;}
      EndIf
    }
  }
}
Include "resolution.pro";
////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_e; NameOfFormulation Ez;
    Quantity {
      { Name u              ; Value { Local{ [ { u} ] ; In Omega; Jacobian JVol; } } }
      { Name EigenValuesReal; Value { Local{ [$EigenvalueReal]; In Point_print; Jacobian JVol; } } }
      { Name EigenValuesImag; Value { Local{ [$EigenvalueImag]; In Point_print; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name postop_e; NameOfPostProcessing postpro_e;
    Operation {
      Print [EigenValuesReal, OnElementsOf Point_print, Format TimeTable, File "EigenValuesReal.txt"];
      Print [EigenValuesImag, OnElementsOf Point_print, Format TimeTable, File "EigenValuesImag.txt"];
      Print [ u , OnElementsOf Omega, File "u.pos", EigenvalueLegend];
    }
  }
}
