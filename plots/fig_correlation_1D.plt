set terminal epslatex size 8.5cm,6cm standalone color 8 solid header \
"\\usepackage{amsmath}"
set output "fig_correlation_1D.tex"

set xlabel '$|\mathbf{r}|$'
set ylabel '$C(\mathbf{r})$'
set format xy '$10^{%T}$'

set logscale xy

f_corr(x,a) = (1+x**2)**(-a/2.)

a_str = "0.50"
a(i)=word(a_str,i)+0.
s(i)=word(s_str,i)
size=words(a_str)

set key bottom left

color(i)=word("1 2 3 7", i)
type(i)=word("4 6 8 10", i)+0.

set xrange [1:]
set yrange [:1]

filename(a)=sprintf("../output/correlation_1D_l21_a%.1f.dat",a)

plot \
for [i=1:size] filename(a(i)) u 1:2:3 w e title sprintf('$a = %.2f$',a(i)) pt type(i) lc color(i) ,\
for [i=1:size] f_corr(x,a(i)) w l notitle ls color(i) lw 2

