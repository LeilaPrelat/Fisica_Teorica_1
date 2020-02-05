## Proyect for Final exam of electrodynamics (Física Teórica 1) 

Here are the codes needed to make the graphics of the mentioned proyect and you can view it in this link: 
https://www.overleaf.com/read/vwdnjwchgqwc

This proyect studies the behavior of a cylinder with 3 or 2 differents mediums. A system of equations (matrix) was create using boundary conditions (the cladding of the cylinder is considerer infinite). The procedure consists on the minimization of the determinant of that matrix and Nelder-Mead was the minimization routine used. The objective of this work was to reproduce the graphs and formulas of two papers: 

1) Amr A. E. Saleh and Jennifer A. Dionne, "Waveguides with a silver lining: Low threshold gain and giant modal gain
in active cylindrical and coaxial plasmonic devices"

2) C. A. Pfeiffer and E. N. Economou, "Surface polaritons in a circularly cylindrical interface: Surface plasmons" **(in construction)**

For ejecuting one of any of the codes you will need the dielectric functions for each medium, so don't forget to download the folder 'dielectric functions'. The code 'determinante.py' is the determinant of the matrix mentioned below and, as an example, the combination of the 3 mediums in 'determinante.py' is: Ag/Si02/Ag. This code only needs 'dielectric_functions.py' to run. 

The code 'disp_relation.py' was used to obtain the dispersion relations and to observe the behavior of the determinant of the matrix. This code uses 'determinante.py' and 'dielectric_functions.py'. The code 'campos.py', like it name says, makes the graphs of the electric field and it needs 'dielectric_functions.py' to run. The code 'campos.py' also uses the archives '.txt' of the dispersion relations which are inside the folder 'Ag/SiO2/Ag'.  
