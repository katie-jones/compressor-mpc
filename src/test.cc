#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

#include "compressor.h"

using namespace std;

typedef Compressor Cp;

Cp::CompressorInput uinit =
    ((Cp::CompressorInput() << 0.304, 0.405, 0.393, 0, 1).finished());

Cp::CompressorState xinit =
    ((Cp::CompressorState() << 0.898, 1.126, 0.15, 440, 0).finished());

Compressor comp(xinit, uinit);
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

