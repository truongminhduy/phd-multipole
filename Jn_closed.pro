//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core_closed.pro";

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
  { Name postpro_Jn; NameOfFormulation Ez;
    Quantity {
      // { Name Jn ; Value { Integral { [ CompZZ[source[]]*{u} ]; In Omega; Integration Int_1; Jacobian JVol; } } }
      { Name Jn ; Value { Integral { [ sourcep[]*{u} ]; In Point_source; Integration Int_1; Jacobian JLin; } } }
    }
  }
}

PostOperation {
  { Name postop_Jn; NameOfPostProcessing postpro_Jn ;
    Operation {
      // Print [  Jn[Omega], OnElementsOf Point_print, Format TimeTable, File "Jns.txt"];
      Print [  Jn[Point_source], OnElementsOf Point_print, Format TimeTable, File "Jns.txt"];
    }
  }
}
