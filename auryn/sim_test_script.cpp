#include "auryn.h"

using namespace auryn;
int main(int ac, char* av[]) {
  auryn_init(ac, av, "test");

  int n_input_neurons = 100;
  float poisson_rate = 1.0;

  PoissonGroup* poisson = new PoissonGroup(n_input_neurons, poisson_rate);
  IzhikevichGroup* neuron = new IzhikevichGroup(1);

  float weight = 0.2;
  AllToAllConnection* conn = new AllToAllConnection(poisson, neuron, weight);
  conn->set_transmitter(GLUT);

  SpikeMonitor* input_spike_mon = new SpikeMonitor(poisson, sys->fn("input", "ras"));
  SpikeMonitor* output_spike_mon = new SpikeMonitor(neuron, sys->fn("output", "ras"));
  VoltageMonitor* output_voltage_mon = new VoltageMonitor(neuron, 0, sys->fn("output", "mem"));

  sys->run(2);

  auryn_free();
}