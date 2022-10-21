#include <iostream>
#include <fstream>
#include <mpi.h>
#include <stdio.h>
#include <cmath>

using namespace std;

int main(int argc, const char* argv[])
{
    int nrows, ncols;
    double *my_matrix;
    double tmp;

    ifstream file;

    file.open("matrix.txt");

    if (file.is_open())
    {
        file >> nrows;

        file >> ncols;

        // MPI_Init: inicialización MPI

        MPI_Init(NULL, NULL);

        int world_size, world_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        // my_nrows: número de filas de cada procesador

        int my_nrows;
        int first_nrows = nrows/world_size;
        int last_nrows = nrows/world_size + nrows%world_size;

        if (world_rank == world_size - 1)
        {
            my_nrows = last_nrows;
        }

        else
        {
            my_nrows = first_nrows;
        }

        // my_matrix: partes de la matriz de cada procesador

        int my_firstrow = world_rank * first_nrows;

        cout << "proceso " << world_rank + 1 << "/" << world_size << ": filas " << my_firstrow << " -> " << my_firstrow + my_nrows - 1 << endl;

        my_matrix = new double [my_nrows * ncols];

        for (int i=0; i<(my_firstrow)*ncols; i++) 
        {
            file >> tmp;
        }

        for (int i=0; i<my_nrows*ncols; i++) 
        {
            file >> my_matrix[i];
        }

        file.close();

        // my_b_0: parte del vector b_0 de cada procesador

        double my_b_0[my_nrows], mu_k;

        for (int i = 0; i < my_nrows; ++i)
        {
            my_b_0[i] = 1;
        }

        // counts: cantidades de elementos a enviar por procesador

        int counts[world_size];
        for (int i=0; i<world_size; i++)
        {
            if (i == world_size - 1)
            {
                counts[i] = last_nrows;
            }

            else
            {
                counts[i] = first_nrows;
            }
        }
        
        // displacements: desplazamientos de las posiciones de los elementos a enviar por procesador

        int displacements[world_size];
        for (int i=0; i<world_size; i++)
        {
            displacements[i] = i*first_nrows;
        }

        // buffer: dirección donde recibir todos los datos

        double buffer[nrows] = {0};

        // my_values: valores de b_k a enviar por procesador
        // my_values_count: cantidad de valores a enviar por procesador

        double* my_values;
        int my_values_count  = my_nrows;

        my_values = (double*)malloc(sizeof(double) * my_values_count);

        for (int i = 0; i < my_nrows; ++i)
        {
            my_values[i] = my_b_0[i];
        }

        // n_iter: número de iteraciones a efectuar

        int n_iter;

        if ( argc > 1 ) 
        {
            n_iter = atoi( argv[1] );
        }        
        
        // MPI_Allgatherv: comunicación entre los procesos

        for (int iter = 0; iter < n_iter; ++iter)
        {
            MPI_Allgatherv(my_values, my_values_count, MPI_DOUBLE, buffer, counts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);

            // my_Ab_k: parte del vector Ab_k por procesador

            double my_Ab_k[my_nrows];

            for (int i = 0; i < my_nrows; ++i)
            {
                my_Ab_k[i] = 0;
            }

            for (int i = 0; i < nrows*my_nrows; ++i)
            {   
                my_Ab_k[i/nrows] += my_matrix[i]*buffer[i%(nrows)];
            }

            // norm: norma del vector Ab_k completo

            double norm = 0;

            double buffer_[nrows] = {0};

            double* my_values_;

            my_values_ = (double*)malloc(sizeof(double) * my_values_count);

            for (int i = 0; i < my_nrows; ++i)
            {
                my_values_[i] = my_Ab_k[i];
            }

            MPI_Allgatherv(my_values_, my_values_count, MPI_DOUBLE, buffer_, counts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);

            for (int i = 0; i < nrows; ++i)
            {   
                norm += buffer_[i]*buffer_[i];
            }
                
            norm = sqrt(norm);

            // my_values: valores de b_k+1 a enviar por procesador

            for (int i = 0; i < my_nrows; ++i)
            {
                my_values[i] = my_Ab_k[i]/norm;
            }

            if (iter == n_iter - 1)
            {
                mu_k = norm;
            }
        }

        // MPI_Finalize: finalización MPI

        MPI_Finalize();

        if (world_rank == world_size - 1)
        {
            // b_k

            for (int i = 0; i < my_nrows; ++i)
            {   
                //cout << "proceso " << world_rank + 1 << "/" << world_size << ": b_k[" << i << "] = " << my_values[i] << endl;
            }

            // mu_k
            
            cout << "mu_k " << mu_k << endl;

            // error

            double mu_teorico = 10;
            double error = abs(mu_k - mu_teorico);
            
            cout << "error método: " << error << endl;
        }
    }    
    
    else
    {
        cout << "Unable to open file." << endl;
    }

    delete[] my_matrix;

    return 0;
}