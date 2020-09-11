/////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core_closed.pro";
////////////////////////////////////////////////////////////////////////
Formulation {
  {Name Hz_main; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace Hgrad;}
    }
    Equation {
      Galerkin { [k0^2*chi[]*Dof{u}   , {u}   ]; In Omega; Jacobian JVol; Integration Int_1;  }
      Galerkin { [-xsi[]    *Dof{d u} , {d u} ]; In Omega; Jacobian JVol; Integration Int_1;  }
      // Galerkin { [-source[]           , {u}   ]; In Omega; Jacobian JVol; Integration Int_1;  }
      Galerkin { [-sourcep[]          , {u}   ]; In Point_source; Jacobian JLin; Integration Int_1;  }
    }
  }
  // {Name Hz_main; Type FemEquation;
  //   Quantity {
  //     { Name u; Type Local; NameOfSpace Hgrad_perp;}
  //   }
  //   Equation {
  //     Galerkin { [k0^2*epsilon[]*Dof{u}   , {u}   ]; In Omega; Jacobian JVol; Integration Int_1;  }
  //     Galerkin { [-1/mur[]      *Dof{d u} , {d u} ]; In Omega; Jacobian JVol; Integration Int_1;  }
  //     Galerkin { [-source[]               , {u}   ]; In Omega; Jacobian JVol; Integration Int_1;  }
    
  //   }
  // }
}

Resolution {
  { Name Scattering;
    System {
      { Name S; NameOfFormulation Hz_main; Type ComplexValue; Frequency Freq;}
    }
    Operation {
      Generate[S];
      Solve[S];
      SaveSolutions[S];
    }
  }
}

////// DATA PROCESS ////////////////////////////////////////////////////////////
PostProcessing {
  { Name get_scat; NameOfFormulation Hz_main;
    Quantity {
      // { Name us  ; Value { Local { [ CompZ[{u}]  ] ; In Omega; Jacobian JVol; } } }  
      { Name us  ; Value { Local { [ {u}  ] ; In Omega; Jacobian JVol; } } }    
      { Name E2 ; Value { Integral { [ Norm[{u}]] ; In Omega_in; Integration Int_1 ; Jacobian JVol ; } } }      
    }
  }
}

PostOperation {
  { Name postop_scat; NameOfPostProcessing get_scat ;
    Operation {
      // Print[ us, OnElementsOf Omega, File "us.pos" ];
      Print[ us, OnElementsOf Omega_plot, File "us.pos" ];
      Print[ E2[Omega_in], OnGlobal, File "E2.txt" , Format Table  ];
      Print[ us, OnPoint {0,0,0}, Format TimeTable, File "E0_0.txt"];
      Print[ us, OnPoint {xD1,yD1,0}, Format TimeTable, File "E0_1.txt"];
      Print[ us, OnPoint {xD2,yD2,0}, Format TimeTable, File "E0_2.txt"];
    }
  }
}
///////////////////////////////////////////////////////////// END //////////////
////////////////////////////////////////////////////////////////////////////////
