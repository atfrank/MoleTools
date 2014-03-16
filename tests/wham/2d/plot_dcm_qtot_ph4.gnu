# Write contours of the PES to a (temporary) file

set term table
set out 'contour.dat'

set contour

set cntrparam bspline
set cntrparam levels 120
set cntrparam order 10

#set cntrparam levels discr 0,0.5,1,1.5,2,3,5
#set cntrparam levels incremental  0,1,8
#set cntrparam levels discr  0,0.6,1.2,1.8,3.0,6.0
#set cntrparam levels incremental  0,1,10
#set cntrparam levels discr 0,0.5,1,1.5,2,3,5,10,15,20

#unset surface
#splot "test.2d.pmf.ph4"
set out; set term pop

# Change single blank lines to double blank lines
!awk "NF<2{printf\"\n\"}{print}" <contour.dat>contour1.dat

reset

set size 1,1
# Uncomment the following to line up the axes
set lmargin 1

set size 1,1 
set origin 0,0

set yrange [11.3717:77.6992]
set xrange [52:372]

set zrange [0.0:15.0]
set cbrange [0.0:8.0]

# Draw the plo
#set palette rgbformulae 33,13,10

set style line 1 lt 1 lw 2 pal

set pm3d map
set pm3d explicit

#set xlabel "Asp356 CA-CYS403 CA Distance (Angs.)"
#set ylabel "Trp 354 Chi2,1 Torsion (degrees)"

set key off

set size square
set origin 0,0

splot 'test.2d.pmf.ph4' notitle with pm3d,\
      'contour1.dat' notitle w l lt -1 lw 1

#l '<./cont.sh contour1.dat 0 15 0'
#splot 'wham.pmf' notitle with pm3d,\
#      'contour1.dat' notitle with lines ls 1, \
#       'test' u 2:3:4 w l lw 3,'new_path.dat' u 2:3:4 w l lw 3, \
#   'b' u 2:3:4 w l lw 4


set terminal postscript eps color lw 2 "Helvetica" 18

set out "test.2d.pmf.ph4.pdf"
replot
set encoding iso_8859_1
set term pop

# remove all customization
reset

pause -1

