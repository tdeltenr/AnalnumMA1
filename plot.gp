reset

fname = "data.txt" # Trouver le fichier des datas

set title "Convergence of the two-grid method"

# Supprimer les réglages des échelles automatiques
unset autoscale x
unset autoscale y

# Définir une échelle semi-logarithmique pour l'axe des x
set logscale x

set xlabel "Number of iterations"
set ylabel "Residual"

plot fname u 1:2 with points title "Residual vs. Number of iterations"

reset 
fname = "data.txt" # Trouver le fichier des datas

set title "Température en °C dans la pièce"

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set xlabel "L [m]"
set ylabel "l [m]"

set zlabel "Température [°C]" # Ajoutez une étiquette pour l'axe z

set palette defined(0 "blue", 10 "yellow", 20 "red")

set pm3d # Active le tracé en relief en 3D
splot fname u 1:2:3 with pm3d # Utilisez "fname" pour spécifier le fichier de données
