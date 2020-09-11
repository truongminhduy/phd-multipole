//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "Pns.dat";
Function{
  For n In {0:neig-1}
    P~{n}[] = Complex[ Pns_re~{n}, Pns_im~{n} ];
  EndFor
}
Include "core.pro";
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
  { Name postpro_modal; NameOfFormulation Ez;
    Quantity {
      { Name up   ;
        Value {
            For i In {0:neig-1}
              Local { [ P~{i}[] * {u}[neig-1-i] ]; In Omega; Jacobian JVol; }
              // Local { [ P~{i}[] * CompZ[{u}[neig-1-i]] ]; In Omega; Jacobian JVol; }
            EndFor
        }
      }
    }
  }
}

PostOperation {
  { Name postop_1; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up1.pos", LastTimeStepOnly ];
    }
  }
  { Name postop_2; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up2.pos", LastTimeStepOnly ];
    }
  }
  { Name postop_3; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up3.pos", LastTimeStepOnly ];
    }
  }
  { Name postop_4; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up4.pos", LastTimeStepOnly ];
    }
  }
  { Name postop_5; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up5.pos", LastTimeStepOnly ];
    }
  }
  { Name postop_0; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega_plot, File "up0.pos", LastTimeStepOnly ];
    }
  }
}
