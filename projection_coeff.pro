//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "Pns.dat";
Include "core.pro";
//////////////////////////////////////////////////////////////////////////////
Formulation {
  {Name Ez; Type FemEquation;
    Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad;}}
    Equation {
      Galerkin { [ Dof{u}       , {u} ]; In Omega; Jacobian JVol; Integration Int_1;}
      Galerkin { [-Field[XYZ[]] , {u} ]; In Omega; Jacobian JVol; Integration Int_1;}
    }
  }
}

Resolution {
  { Name Projection;
    System{
      { Name S; NameOfFormulation Ez ; Type ComplexValue; }
    }
    Operation{
      GmshRead["up0.pos"];
      GmshRead["up1.pos"];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_modal; NameOfFormulation Ez;
    Quantity {            
      { Name uss; Value { Local { [ Field[XYZ[]]{0} + Field[XYZ[]]{1} ] ; In Omega_plot; Jacobian JVol; } } }  
      { Name E2 ; Value { Integral { [ Norm[ Field[XYZ[]]{0} + Field[XYZ[]]{1} ]] ; In Omega_in; Integration Int_1 ; Jacobian JVol ; } } }     
    }
  }
}

PostOperation {
  { Name postop_modal; NameOfPostProcessing postpro_modal ;
    Operation {
      Print[ uss, OnElementsOf Omega_plot, File "uss.pos" ];
      Print[ E2[Omega_in], OnGlobal, File "E2.txt" , Format Table  ];
      Print[ uss, OnPoint {0,0,0}, Format TimeTable, File "E0_0.txt"];
      Print[ uss, OnPoint {xD1,yD1,0}, Format TimeTable, File "E0_1.txt"];
      Print[ uss, OnPoint {xD2,yD2,0}, Format TimeTable, File "E0_2.txt"];
    }
  }
}
