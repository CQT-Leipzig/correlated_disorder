set terminal epslatex standalone size 4.5cm,4.0cm color solid
set output "fig_snapshot_continous.tex"


set pm3d map
set size ratio -1
load "./matlab.palette"
#set palette model RGB #standard RGB colors
unset ytics
unset xtics
#unset ztics
set palette negative
set cbrange [-4:4]
set format cb '$\makebox[0.4em][r]{\tiny %.0f}$'

set lmargin at screen 0.05
set bmargin at screen 0.05
set tmargin at screen 0.95
set rmargin at screen 0.84

splot "./continuous_l8_a0.5_p0.521_seed1000.dat" u 1:2:($3) matrix with image ti ""
