# Set the terminal to PNG for saving plots
set terminal pngcairo size 800,600

# Plot dNd3p vs pT
set output './plotting/dNd3p_vs_pT.png'
set title 'dNd3p vs pT'
set xlabel 'pT'
set ylabel 'dNd3p'
set grid
plot './output/pi_plus.dat' using ($3 != 0 ? $1 : 1/0):($3 != 0 ? $4 : 1/0) with linespoints title 'dNd3p vs pT'

# Plot dNd3p vs phi_p
set output './plotting/dNd3p_vs_phi_p.png'
set title 'dNd3p vs phi_p'
set xlabel 'phi_p'
plot './output/pi_plus.dat' using 2:4 with linespoints title 'dNd3p vs phi_p'

# Plot dNd3p vs y_p
set output './plotting/dNd3p_vs_y_p.png'
set title 'dNd3p vs y_p'
set xlabel 'y_p'
plot './output/pi_plus.dat' using 3:4 with linespoints title 'dNd3p vs y_p'

# Reset output
set output
