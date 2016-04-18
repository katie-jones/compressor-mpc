#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

#include "compressor.h"
#include "compressor_simulation.h"

using namespace std;

typedef Comp Cp;

CompressorSimulation comp;
ofstream statefile;

void Callback(const Cp::CompressorState x, const double t) {
  // cout << comp.u(0) << endl;
  for (int i=0; i<Cp::n_states; i++)
    statefile << x[i] << "\t";
  statefile << "\n";
}

int main(int argc, char *argv[]) {
  statefile.open("states.txt");
  IntegrateCompressor(comp,0,10,0.05,Callback);

  statefile.close();

  return 0;
}

