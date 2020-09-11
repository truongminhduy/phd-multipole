# DQNM expansion in 2D with multipole permittivity

*Application of DQNM expansion of electromagnetic fields in the 2D structure for both bounded and unbounded case: Difficulty encountered in the case of multipole permittivity*

*Python and C++(for getdp and gmsh) code that reproduce the numerical results in the PhD manuscript [arXiv:2009.01307](https://arxiv.org/abs/2009.01307)*

## Affiliation

Minh Duy TRUONG: [minhduy.truong@fresnel.com](minhduy.truong@fresnel.com)

Aix Marseille Univ, CNRS, Centrale Marseille, Institut Fresnel, F-13013 Marseille, France
Athena team

## Instruction

1. Install getdp and gmsh and update the path leading to getdp and gmsh in the main files.
2. Run 4 main files "main_Si_closed_1pole.py", "main_Si_closed_4pole.py", "main_Si_closed_4pole_TM.py" and "main_Si_PML_4pole.py"  in order to obtain data.
3. Run 3 files "plot_epsi.py", "plot_E_1pole.py", "plot_EandH.py" to plot the figures from data.

### Remark
+ Data are saved in 4 folders "pole1_close_TE", "pole4_close_TE", "pole4_close_TM", and "pole4_PML_TE".