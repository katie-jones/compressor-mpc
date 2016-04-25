#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

#include "compressor.h"
#include "parallel_compressors.h"
#include "compressor_simulation.h"

using namespace std;

typedef Comp Cp;
typedef ParallelCompressors PC;

CompressorSimulation comp;
ParallelCompressors compsys;
ofstream statefile;
ofstream sys_states;

void Callback(const Cp::CompressorState x, const double t) {
  // cout << comp.u(0) << endl;
  for (int i = 0; i < Cp::n_states; i++) statefile << x[i] << "\t";
  statefile << "\n";
}

void Callback2(const PC::SysState x, const double t) {
  // cout << comp.u(0) << endl;
  for (int i = 0; i < PC::n_states; i++) sys_states << x[i] << "\t";
  sys_states << "\n";
}

int main(int argc, char *argv[]) {
  statefile.open("states.txt");
  sys_states.open("sysstates.txt");
  // IntegrateCompressor(comp,0,10,0.05,Callback);
  IntegrateSystem(compsys, 0, 0.10, 0.05, Callback2);

  statefile.close();
  sys_states.close();

  return 0;
}
