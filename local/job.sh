python3 generate_matrix.py

echo -n "Ingrese número de iteraciones: "
read n_iteraciones

echo -n "Ingrese número de procesadores: "
read n_procesadores

mpic++ main.cpp -std=c++11
time mpirun -np $n_procesadores ./a.out $n_iteraciones