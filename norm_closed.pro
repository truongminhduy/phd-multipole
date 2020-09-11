//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
//////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////
Function{
  epsilon1[mid] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon1[dif] = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  mur[mid] = TensorDiag[1,1,1];
  mur[dif] = TensorDiag[1,1,1];
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
  { Name postpro_norm; NameOfFormulation Ez;
    Quantity {
      If (E_or_H == 0)
        { Name norm1 ; Value { Integral { [ chi[]*{u}*{u} ]; In Omega_out; Integration Int_1; Jacobian JVol; } } }
        { Name norm2 ; Value { Integral { [ {u}*{u} ]      ; In Omega_in; Integration Int_1; Jacobian JVol; } } }
        { Name norm3 ; Value { Integral { [ Inv[mur[]]*{d u}*{d u} ]; In Omega; Integration Int_1; Jacobian JVol; } } }
      Else
        { Name norm1 ; Value { Integral { [ chi[]*{u}*{u} ]; In Omega   ; Integration Int_1; Jacobian JVol; } } }
        { Name norm2 ; Value { Integral { [ {d u}*{d u}   ]; In Omega_in; Integration Int_1; Jacobian JVol; } } }
      EndIf
    }
  }
}
PostOperation {
  { Name postop_norm; NameOfPostProcessing postpro_norm ;
    Operation {
      If (E_or_H == 0)
        Print [ norm1[Omega_out], OnElementsOf Point_print, Format TimeTable, File "norm1.txt"];
        Print [ norm2[Omega_in ], OnElementsOf Point_print, Format TimeTable, File "norm2.txt"];
        Print [  norm3[Omega], OnElementsOf Point_print, Format TimeTable, File "norm3.txt"];
      Else
        Print [ norm1[Omega   ], OnElementsOf Point_print, Format TimeTable, File "norm1.txt"];
        Print [ norm2[Omega_in], OnElementsOf Point_print, Format TimeTable, File "norm2.txt"];
      EndIf
    }
  }
}