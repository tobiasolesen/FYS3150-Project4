//Program to solve the two-dimensional Ising model.
//The coupling constant J = 1
//Boltzmann's constant k = 1, temperature thus has dimension energy
//Metropolis sampling is used. Periodic boundary conditions (PBC).
//Temperature T is in range [1, 3].
//The array w[17] contains values of dE spanning from -8J to 8J, and it is precalculated in the main part for every new temp.
//The lattice is square.

//In the Metropolis function I loop over all spins, but choose the lattice positions x and y randomly.
//If the move is accepted after performing the Metropolis test, I update the energy and the magnetization.
//The new values are then used to update the averages computed in the main function.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
//#include "lib.h"
using namespace std;
ofstream ofile;

//Inline function for periodic boundary conditions:
inline int periodic(int i, int limit, int add){
    return (i + limit + add) % (limit);
}

//Prototyping functions:

// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&);

void initialize(int n_spins, double temp, int **spin_matrix, double &E, double &M);

void Metropolis(int, long&, int**, double&, double&, double *);

// prints to file the results of the calculations
void output(int, int, double, double *);

random_device rd;
mt19937 engine(rd());
uniform_real_distribution <double> randd(0., 1.0);

//Random numbers:
double rrandom(){
    return randd(engine);
}

double r = rrandom();

//main program:
int main(int argc, char* argv[]){
    cout << r << endl;
    char *outfilename;
    long idum;
    double temp;
    int ** spin_matrix, n_spins, mcs; //Matrise med alle spins (verdier +1 eller -1), antall spins i én retning, antall monte carlo cycles
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

    // Read in output file, abort if there are too few command-line arguments
      if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] <<
          " read also output file on same line" << endl;
        exit(1);
      }
      else{
        outfilename=argv[1];
      }
      ofile.open(outfilename);
      //    Read in initial values such as size of lattice, temp and cycles
      read_input(n_spins, mcs, initial_temp, final_temp, temp_step);

    //spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int)); //Gir feil
    idum = -1; //random starting point
    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){
        E = M = 0; //Initialize energy and magnetization
    }

    //Set up array for possible energy changes:
    for (int de = -8; de <= 8; de++) w[de+8] = 0;
    for (int de = -8; de <= 8; de += 4) w[de+8] = exp(-de/temp);
    //Initialize array for expectation values:
    for (int i = 0; i < 5; i++) average[i] = 0.; //Setter alle elementer lik null. [E, E^2, M, M^2, abs(M)]

    initialize(n_spins, temp, spin_matrix, E, M); //Kaller initialize

    //Start Monte Carlo computation:
    for (int cycles = 1; cycles <= mcs; cycles++){
        //Metropolis(n_spins, idum, spin_matrix, E, M, w);
        //Update expectation values:
        //average[0] += E; average[1] += E*E;
        //average[2] += M; average[3] += M*M; average[4] += fabs(M);
    }

    //Print results:
    output(n_spins, mcs, temp, average); //temp er kalt temperature i Mortens program

    //free_matrix((void **) spin_matrix ); //free memory
    ofile.close(); //close output file
    return 0;

}//End main


//All the needed functions:

//Function to initialize energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M){
    //Setup for spin matrix and initial magnetization
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){
            if (temp < 1.5) spin_matrix[y][x] = 1; //Spin orientation for the ground state. Sets ll spins = 1 if temp is low (T<1.5).
            M += (double) spin_matrix[y][x]; //Initialize magnetization M
            E -= (double) spin_matrix[y][x] * (spin_matrix[periodic(y, n_spins, -1)][x] + spin_matrix[y][periodic(x, n_spins, -1)] );
        }
    }

}//End function initialize

// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
        double& final_temp, double& temp_step)
{
  cout << "Number of Monte Carlo trials =";
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal) =";
  cin >> n_spins;
  cout << "Initial temperature with dimension energy=";
  cin >> initial_temp;
  cout << "Final temperature with dimension energy=";
  cin >> final_temp;
  cout << "Temperature step with dimension energy=";
  cin >> temp_step;
} // end of function read_input

void output(int n_spins, int mcs, double temperature, double *average)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;
} // end output function

/*
//Function that performs the Metropolis algo:
void Metropolis(int n_spins, long & idum, int **spin_matrix, double & E, double & M, double *w){
    //Loop over all spins:
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){
            //find random position:
            int ix = (int) (ran1&idum)*(double) n_spins;
            int iy = (int) (ran1&idum)*(double) n_spins;
            //Calculate energy difference:
            int deltaE = 2*spin_matrix[iy][ix] *
            (spin_matrix[iy][periodic(ix, n_spins, -1)] +
            spin_matrix[periodic(iy, n_spins, -1)][ix] +
            spin_matrix[iy][periodic(ix, n_spins, 1)] +
            spin_matrix[periodic(iy, n_spins, 1)][ix] );

            //Performing the Metropolis test:
            if (ran1 (&idum) <= w[deltaE+8] ){
                spin_matrix[iy][ix] *= -1; //Flip one spin and accept new spin config
                //Update energy and magnetization:
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
            }
        }
    }
}//End Metropolis function (Metropolis sampling over spins)

*/

