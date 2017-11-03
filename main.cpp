//Program to solve the two-dimensional Ising model.
//The coupling constant J = 1
//Boltzmann's constant k = 1, temperature thus has dimension energy
//Metropolis sampling is used. Periodic boundary conditions (PBC).
//Temperature T is in range [1, 3].
//The array w[17] contains values of dE spanning from -8J to 8J, and it is precalculated in the main part for every new temp.
//The lattice is square.

//The program takes as input the initial temperature, final temperature, a temperature step, the number of spins in one direction
//and the number of monte carlo cycles

//In the Metropolis function I loop over all spins, but choose the lattice positions x and y randomly.
//If the move is accepted after performing the Metropolis test, I update the energy and the magnetization.
//The new values are then used to update the averages computed in the main function.

//argc (argument count) is the number of variables pointed to by argv (argument vector).
//This will (in practice) be 1 plus the number of arguments. argv inneholder argumentene.

//Results with mcs = 100000: Eavg = -1.99918, mcs = 1000000: Eavg = -1.99932 og Mavg = 0.999846, mcs = 10000000: Eavg = -1.99934, mcs = 100000000: Eavg = -1.99933
//mcs = 1000000000: Eavg = -1.99933.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
//#include "lib.h"
using namespace std;
ofstream ofile;

//Inline function for periodic boundary conditions:
inline int periodic(int i, int limit, int add){
    return (i + limit + add) % (limit);
}

//Prototyping functions:

double analyticEavg(double, double);
double analyticE2avg(double, double);
double analyticMavg(double, double);
double analyticM2avg(double, double);
double analyticSusceptibility(double, double, double);
double analyticHeatCapacity(double, double, double);

// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&); //Bruker ikke denne

void initialize(int n_spins, double temp, int **spin_matrix, double &E, double &M);

void Metropolis(int, int**, double&, double&, double *);

// prints to file the results of the calculations
void output(int, int, double, double *, vector <double>);

random_device rd;
mt19937 engine(rd());
uniform_real_distribution <double> randd(0.0, 1.0);

//Random numbers:
double rrandom(){
    return randd(engine);
}

//main program:
int main(){

    //char *outfilename;
    double temp = 1.;
    double initial_temp = 1.;
    double final_temp = 3.;
    double temp_step = 0.2;
    double beta = 1./temp;  //1/kT med k=1
    //int ** spin_matrix, n_spins, mcs; //Matrise med alle spins (verdier +1 eller -1), antall spins i én retning, antall monte carlo cycles
    int n_spins = 2; //Antall spins i én retning
    int mcs = 10000000; //antall monte carlo cycles should be a billion?
    double w[17], average[5], E, M;
    double Z = 12 + 2*exp(8*beta) + 2*exp(-8*beta); //For 2x2-lattice
    vector <double> E_vec(mcs);

    int** spin_matrix = new int*[n_spins]; //Deklarerer matrisen spin_matrix
    for(int i = 0; i < n_spins; i++)
        spin_matrix[i] = new int[n_spins];

    //Analytisk forventningsverdi av energi og magnetisering:
    double analytic_Eavg = analyticEavg(Z, beta); //function call
    double analytic_E2avg = analyticE2avg(Z, beta);
    double analytic_Mavg = analyticMavg(Z, beta);
    double analytic_M2avg = analyticM2avg(Z, beta);
    double analytic_Susceptibility = analyticSusceptibility(analytic_Mavg, analytic_M2avg, temp);
    double analytic_HeatCapacity = analyticHeatCapacity(analytic_Eavg, analytic_E2avg, temp);

    cout << "Analytical expectation value of energy: " << analytic_Eavg/n_spins/n_spins << endl; //Eavg/n_spins/n_spins
    cout << "Analytical expectation value of energy^2: " << analytic_E2avg/n_spins/n_spins << endl; //Eavg/n_spins/n_spins
    cout << "Analytical expectation value of magnetization: " << analytic_Mavg/n_spins/n_spins << endl; //Mavg/n_spins/n_spins
    cout << "Analytical expectation value of magnetization^2: " << analytic_M2avg/n_spins/n_spins << endl; //Mavg/n_spins/n_spins
    cout << "Analytical susceptibility: " << analytic_Susceptibility << endl;
    cout << "Analytical heat capacity: " << analytic_HeatCapacity << endl;

    /*
    // Read in output file, abort if there are too few command-line arguments
      if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] <<
          " read also output file on same line" << endl;
        exit(1);
      }
      else{
        //outfilename=argv[1];
        outfilename = "output.txt";
      }
      */
      //outfilename = "output.txt";
      ofile.open("output.txt");
      //Read in initial values such as size of lattice, temp and cycles
      //read_input(n_spins, mcs, initial_temp, final_temp, temp_step);  //Command line arguments

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){
        E = M = 0; //Initialize energy and magnetization
    }

    //Set up array for possible energy changes:
    for (int de = -8; de <= 8; de++) {
        w[de+8] = 0;
    }
    for (int de = -8; de <= 8; de += 4) w[de+8] = exp(-de/temp);
    //Initialize array for expectation values:
    for (int i = 0; i < 5; i++) average[i] = 0.; //Setter alle elementer lik null. [E, E^2, M, M^2, abs(M)]



    initialize(n_spins, temp, spin_matrix, E, M); //Kaller initialize
    //Start Monte Carlo computation:
    for (int cycles = 1; cycles <= mcs; cycles++){
        Metropolis(n_spins, spin_matrix, E, M, w);
        //Update expectation values:
        average[0] += E; average[1] += E*E;  //average is expectation values. I c) vil jeg plotte E jeg faar etter hver sweep gjennom lattice mot mcs
        average[2] += M; average[3] += M*M; average[4] += fabs(M);
        E_vec[cycles] = E;  //Inneholder E-verdiene jeg skal plotte. Antall E-verdier vil vaere lik mcs.
        //cout << E_vec[cycles] << endl;
    }

    //Print results:
    output(n_spins, mcs, temp, average, E_vec); //temp er kalt temperature i Mortens program

    //free_matrix((void **) spin_matrix ); //free memory
    ofile.close(); //close output file

    //Skriver ut:
    cout << "Temperatur: " << temp << endl;
    cout << "Antall Monte Carlo cycles: " << mcs << endl;

    return 0;

}//End main

//All the needed functions:

//Function to initialize energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M){
    //Setup for spin matrix and initial magnetization
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){
            //if (temp < 1.5) spin_matrix[y][x] = 1; //Spin orientation for the ground state. Sets ll spins = 1 if temp is low (T<1.5).

            //Random start congif:
            double r = rrandom();
            if (r <= 0.5) spin_matrix[y][x] = -1;
            else{
                spin_matrix[y][x] = 1;
            }

            M += (double) spin_matrix[y][x]; //Initialize magnetization M
        }
    }

    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){
            E -= (double) spin_matrix[y][x] * (spin_matrix[periodic(y, n_spins, -1)][x] + spin_matrix[y][periodic(x, n_spins, -1)] );
        }
    }
}//End function initialize
/*
// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
        double& final_temp, double& temp_step)
{
  cin >> mcs;
  cout << "Number of Monte Carlo trials =";
  //cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal) =";
  cin >> n_spins;
  cout << "Initial temperature with dimension energy=";
  cin >> initial_temp;
  cout << "Final temperature with dimension energy=";
  cin >> final_temp;
  cout << "Temperature step with dimension energy=";
  cin >> temp_step;
} // end of function read_input
*/

//Calculates and prints numerical values:
void output(int n_spins, int mcs, double temp, double *average, vector <double> E_vec)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm; //average of E^2
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  double num_HeatCapacity = (E2average - Eaverage*Eaverage)/(temp*temp);
  double num_Susceptibility = (M2average - Mabsaverage*Mabsaverage)/((double) temp);
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temp; //1 kolonne i filen output.txt
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins; //2 kolonne
  ofile << setw(15) << setprecision(8) << Evariance/temp/temp; //3 kolonne
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins; //4 kolonne
  ofile << setw(15) << setprecision(8) << Mvariance/temp; //5 kolonne
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl; //6 kolonne
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;

  cout << "Numerical E average: " << Eaverage/n_spins/n_spins << endl;
  cout << "Numerical abs M average: " << Mabsaverage/n_spins/n_spins << endl;
  cout << "Numerical E^2 average: " << E2average/n_spins/n_spins << endl;
  cout << "Numerical M^2 average: " << M2average/n_spins/n_spins << endl;
  cout << "Numerical heat capacity: " << num_HeatCapacity << endl;
  cout << "Numerical susceptibility: " << num_Susceptibility << endl;

} // end output function

//Function that performs the Metropolis algo:
void Metropolis(int n_spins, int **spin_matrix, double & E, double & M, double *w){
    //Loop over all spins:
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){
            //find random position:
            int ix = (int) (rrandom() * n_spins);
            int iy = (int) (rrandom() * n_spins);
            //Calculate energy difference: (likn 13.6 i forel notater)
            int deltaE = 2*spin_matrix[iy][ix] *
            (spin_matrix[iy][periodic(ix, n_spins, -1)] +
            spin_matrix[periodic(iy, n_spins, -1)][ix] +
            spin_matrix[iy][periodic(ix, n_spins, 1)] +
            spin_matrix[periodic(iy, n_spins, 1)][ix] );

            //Performing the Metropolis test:
            if (rrandom() <= w[deltaE+8] ){
                spin_matrix[iy][ix] *= -1; //Flip one spin and accept new spin config
                //Update energy and magnetization:
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
            }
        }
    }
}//End Metropolis function (Metropolis sampling over spins)

//Function that calculates analytic average energy
double analyticEavg(double Z, double beta){
    double analytic_Eavg = (16*exp(-8*beta) - 16*exp(8*beta))/Z;
    return analytic_Eavg; //Skal bli ca -2
}

//Function that calculates analytic average magnetization
double analyticMavg(double Z, double beta){
    double analytic_Mavg = (8*exp(8*beta) + 16)/Z;
    return analytic_Mavg; //Skal bli ca
}

//Function that calculates analytic average energy squared
double analyticE2avg(double Z, double beta){
    double analytic_E2avg = 256*cosh(8*beta)/Z;
    return analytic_E2avg;
}

//Function that calculates analytic average magnetization squared:
double analyticM2avg(double Z, double beta){
    double analytic_M2avg = (32 + 32*exp(8*beta))/Z;
    return analytic_M2avg;
}

double analyticSusceptibility(double analytic_Mavg, double analytic_M2avg, double temp){
    //double analytical_Susceptibility = 1/temp*( (32 + +32*exp(8*beta))/Z - pow( 16 + 8*exp(8*beta), 2 )/(Z*Z) );
    double analytic_Susceptibility = (analytic_M2avg - analytic_Mavg*analytic_Mavg)/temp;
    return analytic_Susceptibility;
}

double analyticHeatCapacity(double analytic_Eavg, double analytic_E2avg, double temp){
    double analytic_Heat_Capacity = (analytic_E2avg - analytic_Eavg*analytic_Eavg)/(temp*temp);
    return analytic_Heat_Capacity;
}

