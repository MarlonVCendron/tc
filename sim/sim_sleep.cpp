#include "startup.h"

int main(int ac, char *av[]) {
  int ret = startup(ac, av);
  if (ret != 0) return ret;

  std::cout << (sleeps ? "SLEEPING" : "NOT SLEEPING") << '\n';
  std::cout << sleep_duration << '\n';
  std::cout << awake_duration << '\n';

  AIF2Group *neurons_e = new AIF2Group(size);
  neurons_e->dg_adapt1 = 0.1;
  neurons_e->dg_adapt2 = adapt;
  neurons_e->set_tau_ampa(5e-3);
  neurons_e->set_tau_gaba(10e-3);
  neurons_e->set_tau_nmda(100e-3);
  neurons_e->set_ampa_nmda_ratio(0.2);

  IFGroup *neurons_i = new IFGroup(size / 4);
  neurons_i->set_tau_ampa(5e-3);
  neurons_i->set_tau_gaba(10e-3);
  neurons_i->set_tau_nmda(100e-3);
  neurons_i->set_ampa_nmda_ratio(0.3);

  SleepGroup *sleepgroup;
  if (sleeps) {
    sleepgroup = new SleepGroup(size / 16, sleep_duration, awake_duration);
  }

  SleepyStimulusGroup *stimgroup;
  sprintf(strbuf, "%s/%s.%d.stimtimes", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  string stimtimefile = strbuf;
  stimgroup = new SleepyStimulusGroup(size, stimtimefile);
  stimgroup->set_mean_on_period(ontime);
  stimgroup->set_mean_off_period(offtime);
  stimgroup->binary_patterns = true;
  stimgroup->scale = scale;
  stimgroup->background_rate = bgrate;
  stimgroup->background_during_stimulus = true;
  stimgroup->sleep_duration = sleep_duration;
  stimgroup->awake_duration = sleeps ? awake_duration : simtime; // Just a quick way to make the net not sleep
  // stimgroup->randomintensities = true;
  if (seed != 1) stimgroup->seed(seed);

  P11Connection *con_ee = new P11Connection(neurons_e, neurons_e, wee, sparseness, eta, kappa, wmax);
  con_ee->set_transmitter(AMPA);
  con_ee->set_name("E->E");
  con_ee->set_weight_a(weight_a);
  con_ee->set_weight_c(weight_c);
  con_ee->consolidation_active = consolidation;
  con_ee->pot_strength = pot_strength / normalization_factor;
  logger->parameter("normalized pot_strength", con_ee->pot_strength);
  if (noisy_initial_weights) con_ee->random_data(wee, wee);
  if (consolidate_initial_weights) con_ee->consolidate();
  // STP parameters
  con_ee->set_tau_d(taud);
  con_ee->set_tau_f(tauf);
  con_ee->set_ujump(ujump);
  con_ee->set_beta(beta);
  con_ee->delta = raw_delta * eta;
  con_ee->set_min_weight(wmin);
  con_ee->set_tau_hom(tauh);

  STPConnection *con_ei = new STPConnection(neurons_e, neurons_i, 3 * wei, sparseness, GLUT);
  con_ei->set_tau_d(taud);
  con_ei->set_tau_f(0.6);
  con_ei->set_ujump(0.2);

  double geta = -eta * 1e-4;
  SparseConnection *con_ii = new SparseConnection(neurons_i, neurons_i, wii, sparseness, GABA);
  con_ii->set_name("I->I");

  GlobalPFConnection *con_ie;
  con_ie = new GlobalPFConnection(neurons_i, neurons_e, wie, sparseness, 10.0, eta / 50, alpha, wmaxi, GABA);
  con_ie->set_name("I->E");

  SparseConnection *con_sleep_e;
  if (sleeps) {
    con_sleep_e = new SparseConnection(sleepgroup, neurons_e, 0.04, 0.2, GLUT);
    con_sleep_e->set_name("Sleep->E");
  }

  P11Connection *con_stim_e = new P11Connection(stimgroup, neurons_e, wext, sparseness_ext, eta, kappa, wmax, GLUT);
  con_stim_e->set_weight_a(weight_a);
  con_stim_e->set_weight_c(weight_c);
  con_stim_e->set_tau_d(taud);
  con_stim_e->set_tau_f(tauf);
  con_stim_e->set_ujump(ujump);
  con_stim_e->set_beta(beta);
  con_stim_e->delta = raw_delta * eta;
  con_stim_e->set_min_weight(wmin);
  if (noisy_initial_weights) con_stim_e->random_data(wext, wext);
  con_stim_e->set_name("Stim->E");
  con_stim_e->consolidation_active = consolidation;
  con_stim_e->pot_strength = pot_strength / normalization_factor;
  con_stim_e->set_tau_hom(tauh);
  if (consolidate_initial_weights) con_stim_e->consolidate();

  if (!stimfile.empty()) {
    logger->msg("Setting up stimulus ...", PROGRESS, true);
    stimgroup->load_patterns(stimfile.c_str());
    stimgroup->set_next_action_time(10); // let network settle for some time
    // stimgroup->set_next_action_time(2); // let network settle for some time

    sprintf(strbuf, "%s/%s.%d.s.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    BinarySpikeMonitor *smon_s = new BinarySpikeMonitor(stimgroup, string(strbuf), size);

    sprintf(strbuf, "%s/%s.%d.wse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    new WeightStatsMonitor(con_stim_e, string(strbuf));
  }

  // load if necessary
  if (!infilename.empty()) {
    logger->msg("Loading from file ...", PROGRESS, true);
    sys->load_network_state(infilename.c_str());
  }

  if (!prefile.empty() && chi > 0.0) {
    con_ee->patterns_ignore_gamma = true;
    con_ee->load_patterns(prefile, chi);
    // con_ee->consolidate();
  }

  if (!prefile.empty() && xi > 0.0) {
    con_stim_e->patterns_ignore_gamma = true;
    con_stim_e->load_patterns(prefile, xi, false);
    if (consolidate_initial_weights) con_stim_e->consolidate();
  }

  sprintf(strbuf, "%s/%s.%d.wse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  new WeightStatsMonitor(con_stim_e, string(strbuf));

  if (!recfile.empty() && xi > 0.0) {
    con_stim_e->load_fragile_matrix(recfile); // TODO
    con_stim_e->scale_all(xi);
  }

  sprintf(strbuf, "%s/%s.%d.sse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  WeightMonitor *wmon_s = new WeightMonitor(con_stim_e, string(strbuf), 1.0);
  wmon_s->add_equally_spaced(50);
  if (!monfile.empty()) {
    sprintf(strbuf, "%s/%s.%d.pact", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    PatternMonitor *patmon = new PatternMonitor(neurons_e, string(strbuf), monfile.c_str(), 100);

    if (!stimfile.empty())
      wmon_s->load_pattern_connections(stimfile, monfile, 20, 20, ASSEMBLIES_ONLY);
    else
      wmon_s->load_pattern_connections(monfile, 20, 20, ASSEMBLIES_ONLY);
  }

  sprintf(strbuf, "%s/%s.%d.see", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  WeightMonitor *wmon = new WeightMonitor(con_ee, string(strbuf), 1.0);
  wmon->add_equally_spaced(50);

  if (!monfile.empty()) wmon->load_pattern_connections(monfile, 10, 10, ASSEMBLIES_ONLY);

  sprintf(strbuf, "%s/%s.%d.hom", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  StateMonitor *stmon_hom = new StateMonitor(con_ee->hom, 12, string(strbuf), 1);

  sprintf(strbuf, "%s/%s.%d.mem", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  VoltageMonitor *stmon_mem = new VoltageMonitor(neurons_e, 3, string(strbuf));
  stmon_mem->record_for(10); // stops recording after 10s

  sprintf(strbuf, "%s/%s.%d.imem", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  VoltageMonitor *stmon_imem = new VoltageMonitor(neurons_i, 3, string(strbuf));
  stmon_imem->record_for(10); // stops recording after 10s

  sprintf(strbuf, "%s/%s.%d.sie", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  WeightMonitor *wmon_ie = new WeightMonitor(con_ie, string(strbuf));
  wmon_ie->add_equally_spaced(50);

  sprintf(strbuf, "%s/%s.%d.wee", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  new WeightStatsMonitor(con_ee, string(strbuf));

  sprintf(strbuf, "%s/%s.%d.wie", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  new WeightStatsMonitor(con_ie, string(strbuf));

  if (!monfile.empty()) {
    sprintf(strbuf, "%s/%s.%d.wprec", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    WeightPatternMonitor *wpmon = new WeightPatternMonitor(con_ee, string(strbuf), 60);
    wpmon->load_patterns(monfile);
  }

  if (!premonfile.empty() && !monfile.empty()) {
    sprintf(strbuf, "%s/%s.%d.wpin", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    WeightPatternMonitor *wpmon = new WeightPatternMonitor(con_stim_e, string(strbuf), 60);
    wpmon->load_pre_patterns(premonfile);
    wpmon->load_post_patterns(monfile);
  }

  BinarySpikeMonitor *smon_sleep;
  if (sleeps) {
    sprintf(strbuf, "%s/%s.%d.sleep.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
    smon_sleep = new BinarySpikeMonitor(sleepgroup, string(strbuf), size);
  }

  sprintf(strbuf, "%s/%s.%d.e.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  BinarySpikeMonitor *smon_e = new BinarySpikeMonitor(neurons_e, string(strbuf), size);

  sprintf(strbuf, "%s/%s.%d.i.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  BinarySpikeMonitor *smon_i = new BinarySpikeMonitor(neurons_i, string(strbuf), size);

  sprintf(strbuf, "%s/%s.%d.e.prate", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  PopulationRateMonitor *pmon_e = new PopulationRateMonitor(neurons_e, string(strbuf), 0.1);

  sprintf(strbuf, "%s/%s.%d.i.prate", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  PopulationRateMonitor *pmon_i = new PopulationRateMonitor(neurons_i, string(strbuf), 0.1);

  RateChecker *chk = new RateChecker(neurons_e, -1, 20., 0.1);

  stimgroup->set_mean_on_period(ontime);
  stimgroup->set_mean_off_period(offtime);

  if (eta > 0) {
    con_ee->stdp_active = true;
    con_stim_e->stdp_active = true;
    con_ie->stdp_active = true;
  } else {
    con_ee->stdp_active = false;
    con_stim_e->stdp_active = false;
    con_ie->stdp_active = false;
  }

  con_ie->stdp_active = isp_active;

  logger->msg("Main simtime ...", PROGRESS, true);
  if (!sys->run(simtime, false)) errcode = 1;

  logger->msg("Saving ...", PROGRESS, true);
  if (save) {
    sys->set_output_dir(dir);
    sys->save_network_state(file_prefix);
  }

  logger->msg("Writing connectivity matrices ...", PROGRESS, true);
  sprintf(strbuf, "%s/%s.%d.ee.wmat", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  con_ee->write_to_file(strbuf);

  sprintf(strbuf, "%s/%s.%d.ext.wmat", dir.c_str(), file_prefix.c_str(), sys->mpi_rank());
  con_stim_e->write_to_file(strbuf);

  if (errcode) auryn_abort(errcode);

  logger->msg("Freeing ...", PROGRESS, true);
  auryn_free();

  return errcode;
}
