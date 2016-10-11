// Centralized parallel simulation
#include <exception>
#include <sstream>
#include <string>

// Define two macros to use correct typedefs etc. in "common-variables.h"
#define CONTROLLER_TYPE_CENTRALIZED
#define SYSTEM_TYPE_PARALLEL

// Defines the following types:
// SimSystem: simulation system used for numerical integration
// NvCtr: nerve center of controllers
// AugmentedSystem(1,2): (sub)systems with augmented state matrices
// Obsv(1,2): observers of (sub)systems
// Controller(1,2): (sub)controllers
#include "common-variables.h"

// -------------------------------------------------------------------------- //
// ---------------------------- GLOBAL VARIABLES ---------------------------- //
// -------------------------------------------------------------------------- //

int n_solver_iterations;  // Number of solver iterations
int n_timing_iterations;  // Number of solver iterations to time

SimSystem *p_sim_compressor;        // Pointer to the simulation
CompressorSystem *p_compressor;  // Pointer to the parallel compressor sys.
NvCtr *p_controller;                // Pointer to controller nerve center
std::ofstream output_file;          // File where timing data is saved

boost::timer::nanosecond_type time_total;  // Solver time (ns) per iteration

// -------------------------------------------------------------------------- //
// -------------------------------- CALLBACK -------------------------------- //
// -------------------------------------------------------------------------- //

void Callback(CompressorSystem::State x, double t) {
  // Get output from system
  CompressorSystem::Output y = p_compressor->GetOutput(x);

  // Get and apply next input
  NvCtr::ControlInput u = p_controller->GetNextInputWithTiming(
      p_compressor->GetOutput(x), n_timing_iterations, &time_total);

  // Output information to file: t, x, y, u, time_total
  output_file << t << std::endl
              << x.transpose() << std::endl
              << u.transpose() << std::endl
              << y.transpose() << std::endl
              << time_total << std::endl
              << std::endl;

  // Apply next input
  p_sim_compressor->SetInput(u);
}

// -------------------------------------------------------------------------- //
// ----------------------------- DEFAULT SETUP ------------------------------ //
// -------------------------------------------------------------------------- //

// Simulation variables controlled by setup file
struct SimulationVariables {
  int n_solver_iterations = 1;
  int n_timing_iterations = 1;

  std::string folder_name = "parallel/";
  std::string output_fname = "cent_output.dat";

  CompressorSystem::Input u_init = CompressorSystem::GetDefaultInput();
  CompressorSystem::ControlInput u_control_init =
      CompressorSystem::ControlInput::Zero();
  CompressorSystem::State x_init = CompressorSystem::GetDefaultState();

  double sampling_time = 0.05;

  Obsv::ObserverMatrix M =
      (Obsv::ObserverMatrix()
           << Eigen::Matrix<double, CompressorSystem::n_states,
                            CompressorSystem::n_outputs>::Zero(),
       Eigen::Matrix<double, n_disturbance_states,
                     CompressorSystem::n_outputs>::Identity())
          .finished();

  // Weights
  NvCtr::UWeightType uwt = NvCtr::UWeightType::Identity();
  Controller::YWeightType ywt = Controller::YWeightType::Identity();

  CompressorSystem::Output y_ref = CompressorSystem::Output::Zero();

  // Input constraints
  InputConstraints<Controller::n_control_inputs> constraints;

  SimulationVariables() { constraints.use_rate_constraints = true; }
};

// -------------------------------------------------------------------------- //
// ----------------------------- READ VARIABLES ----------------------------- //
// -------------------------------------------------------------------------- //

SimulationVariables ReadSimulationVariables(std::ifstream &setup_file,
                                            int &line_number) {
  SimulationVariables x = SimulationVariables();
  std::string line;
  std::istringstream line_stream;

  // Read file line by line
  while (std::getline(setup_file, line)) {
    line_number++;

    // Skip commented and empty lines
    if ((line[0] == '#') || (line.empty()) ||
        (std::all_of(line.begin(), line.end(), isspace)))
      continue;

    if (line.find("n-iterations") != std::string::npos) {
      ReadNumbers(&x.n_solver_iterations, 1, setup_file, line_number);
    } else if (line.find("n-timing-iterations") != std::string::npos) {
      ReadNumbers(&x.n_timing_iterations, 1, setup_file, line_number);
    } else if (line.find("folder-name") != std::string::npos) {
      x.folder_name = ReadString(setup_file, line_number);
      // Always end with /
      if (x.folder_name.back() != '/') x.folder_name.push_back('/');
    } else if (line.find("output-filename") != std::string::npos) {
      x.output_fname = ReadString(setup_file, line_number);
    } else if (line.find("u-init") != std::string::npos) {
      ReadNumbers(x.u_init.data(), x.u_init.size(), setup_file, line_number);
    } else if (line.find("x-init") != std::string::npos) {
      ReadNumbers(x.x_init.data(), x.x_init.size(), setup_file, line_number);
    } else if ((line.find("Ts") != std::string::npos) ||
               (line.find("sampling-time") != std::string::npos)) {
      ReadNumbers(&x.sampling_time, 1, setup_file, line_number);
    } else if (line.find("observer-matrix") != std::string::npos) {
      ReadNumbers(x.M.data(), x.M.size(), setup_file, line_number);
    } else if (line.find("uwt") != std::string::npos) {
      ReadNumbers(x.uwt.data(), x.uwt.size(), setup_file, line_number);
    } else if (line.find("ywt") != std::string::npos) {
      ReadNumbers(x.ywt.data(), x.ywt.size(), setup_file, line_number);
    } else if (line.find("yref") != std::string::npos) {
      ReadNumbers(x.y_ref.data(), x.y_ref.size(), setup_file, line_number);
    } else if (line.find("constraints-lower") != std::string::npos) {
      ReadNumbers(x.constraints.lower_bound.data(),
                  x.constraints.lower_bound.size(), setup_file, line_number);
    } else if (line.find("constraints-upper") != std::string::npos) {
      ReadNumbers(x.constraints.upper_bound.data(),
                  x.constraints.upper_bound.size(), setup_file, line_number);
    } else if (line.find("constraints-rate-lower") != std::string::npos) {
      ReadNumbers(x.constraints.lower_rate_bound.data(),
                  x.constraints.lower_rate_bound.size(), setup_file,
                  line_number);
    } else if (line.find("constraints-rate-upper") != std::string::npos) {
      ReadNumbers(x.constraints.upper_rate_bound.data(),
                  x.constraints.upper_rate_bound.size(), setup_file,
                  line_number);
    } else if (line.find("simulation") != std::string::npos) {
      break;
    } else {
      std::cerr << "Option \"" << line << "\" on line " << line_number
                << " not recognized." << std::endl;
    }
  }

  return x;
}

int main(int argc, char **argv) {
  std::string setup_fname = "setup";
  std::ifstream setup_file;
  if (argc > 1) {
    std::istringstream ss(argv[1]);
    if (!(ss >> setup_fname)) {
      setup_fname = "setup";
      std::cerr << "Warning: Invalid filename: " << argv[1]
                << ". Using default name " << setup_fname << "." << std::endl;
    }
  }

  setup_file.open(setup_fname);

  int line_number = 0;
  SimulationVariables x = ReadSimulationVariables(setup_file, line_number);
  std::cout << "Folder name: " << x.folder_name << std::endl;
  std::cout << "Output filename: " << x.output_fname << std::endl;
  std::cout << "Initial state: " << x.x_init << std::endl;
  std::cout << "Initial input: " << x.u_init << std::endl;
  std::cout << "Initial control input: " << x.u_control_init << std::endl;
  std::cout << "Observer matrix: " << std::endl << x.M << std::endl;
  std::cout << "Constraints: " << x.constraints.lower_bound << std::endl
            << x.constraints.upper_bound << std::endl
            << x.constraints.lower_rate_bound << std::endl
            << x.constraints.upper_rate_bound << std::endl;

  n_solver_iterations = x.n_solver_iterations;
  n_timing_iterations = x.n_timing_iterations;

  if (n_solver_iterations < 0 || n_solver_iterations > 1) {
    std::cout << "Warning: Number of solver iterations for centralized "
                 "controller should be 0 or 1. Using 1 solver iteration."
              << std::endl;
    n_solver_iterations = 1;
  }

  time_total = 0;
  output_file.open(x.folder_name + x.output_fname);

  CompressorSystem compressor;
  p_compressor = &compressor;

  SimSystem sim_comp(p_compressor, x.u_init, x.x_init);
  p_sim_compressor = &sim_comp;

  // Setup controller
  AugmentedSystem sys(compressor, x.sampling_time);
  Controller ctrl(sys, x.constraints, x.M);

  // Create a nerve center
  std::tuple<Controller> ctrl_tuple(ctrl);
  NvCtr nerve_center(ctrl_tuple, n_solver_iterations);
  p_controller = &nerve_center;

  NvCtr::SubYWeightType ywts(x.ywt);
  nerve_center.SetWeights(x.uwt, ywts);
  nerve_center.SetOutputReference(x.y_ref.replicate<Controller::p, 1>());

  nerve_center.Initialize(x.x_init, x.u_control_init, x.u_init,
                          compressor.GetOutput(x.x_init));

  // Initialize disturbance
  CompressorSystem::Input u_disturbance;
  double t_past = -x.sampling_time;
  double t_next;

  // Time entire simulation
  std::cout << "Running centralized simulation... ";
  std::cout.flush();
  boost::timer::cpu_timer simulation_timer;

  // Read inputs and times from file and simulate
  while (ReadNumbers(u_disturbance.data(), u_disturbance.size(), setup_file,
                     line_number, false) == u_disturbance.size()) {
    if (!(ReadNumbers(&t_next, 1, setup_file, line_number, false))) {
      std::cerr << "Simulation time could not be read." << std::endl;
      break;
    }
    std::cout << "Simulating from time " << t_past << " to time " << t_next
              << " with offset:" << std::endl;
    for (int i = 0; i < u_disturbance.size(); i++) {
      std::cout << u_disturbance(i) << "\t";
    }
    std::cout << std::endl;

    sim_comp.SetOffset(x.u_init + u_disturbance);

    sim_comp.Integrate(t_past + x.sampling_time, t_next, x.sampling_time,
                       &Callback);

    t_past = t_next;
  }

  boost::timer::cpu_times simulation_cpu_time = simulation_timer.elapsed();
  boost::timer::nanosecond_type simulation_ns(simulation_cpu_time.system +
                                              simulation_cpu_time.user);

  std::cout << "Finished." << std::endl
            << "Total time required:\t"
            << static_cast<double>(simulation_ns) / 1.0e6 << " ms." << std::endl
            << std::endl;

  setup_file.close();
  output_file.close();

  return 0;
}
