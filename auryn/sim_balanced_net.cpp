#include "auryn.h"

using namespace auryn;
int main(int ac, char* av[]) {
  auryn_init(ac, av, "plastic");

  int nb_exc_neurons = 20000;
  int nb_inh_neurons = nb_exc_neurons / 4;
  IFGroup* neurons_exc = new IFGroup(nb_exc_neurons);
  neurons_exc->set_name("exc neurons");
  neurons_exc->get_state_vector("g_nmda")->set_random();
  IFGroup* neurons_inh = new IFGroup(nb_inh_neurons);
  neurons_inh->set_tau_mem(5e-3);
  neurons_inh->set_name("inh neurons");

  int nb_input_neurons = 5000;
  float poisson_rate = 2.0;
  PoissonGroup* poisson = new PoissonGroup(nb_input_neurons, poisson_rate);

  float weight = 0.2;       // conductance amplitude in units of leak conductance
  float sparseness = 0.05;  // probability of connection
  SparseConnection* con_ext_exc = new SparseConnection(poisson, neurons_exc, weight, sparseness, GLUT);

  float gamma = 4.0;
  SparseConnection* con_ee = new SparseConnection(neurons_exc, neurons_exc, weight, sparseness, GLUT);
  SparseConnection* con_ei = new SparseConnection(neurons_exc, neurons_inh, weight, sparseness, GLUT);
  SparseConnection* con_ie = new SparseConnection(neurons_inh, neurons_exc, gamma * weight, sparseness, GABA);
  SparseConnection* con_ii = new SparseConnection(neurons_inh, neurons_inh, gamma * weight, sparseness, GABA);

  SpikeMonitor* exc_spike_mon = new SpikeMonitor(neurons_exc, sys->fn("exc", "ras"));
  VoltageMonitor* voltage_mon = new VoltageMonitor(neurons_exc, 0, sys->fn("neuron", "mem"));
  voltage_mon->record_for(2);

  sys->run(10);

  con_ee->set_block(0, 500, 0, 500, 5 * weight);
  sys->run(2);

  auryn_free();
}