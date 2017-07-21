#
set terminal epslatex standalone size 4cm,4cm color solid
set output "fig_snapshot_discrete.tex"


set pm3d map
set size ratio -1
set palette defined (0 "white",1 "black")
unset colorbox
unset ytics
unset xtics
unset ztics

set lmargin at screen 0.05
set bmargin at screen 0.05
set rmargin at screen 0.95
set tmargin at screen 0.95

splot "./discrete_l8_a0.5_p0.521_seed1000.dat" matrix with image ti ""
