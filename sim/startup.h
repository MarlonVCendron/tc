#pragma once

#include "GlobalPFConnection.h"
#include "P11Connection.h"
#include "SleepGroup.h"
#include "SleepyStimulusGroup.h"
#include "auryn.h"

using namespace auryn;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

string dir = "./data/sim";
string file_prefix = "rc";
string infilename = "";

char strbuf[255];
string msg;

bool sleeps = false;
bool save = false;
bool consolidation = true;
bool isp_active = true;

bool noisy_initial_weights = false;
bool consolidate_initial_weights = false;
bool quiet = false;

NeuronID size = 4096;
NeuronID seed = 1;
double alpha = 3;
double kappa = 10;
double tauf = 200e-3;
double ujump = 0.2;
double taud = 200e-3;
double eta = 1e-3;

double tauh = 600.0;

double beta = 5.0e-2;
double delta = 2.0e-2;
double weight_a = 0.1;
double weight_c = 0.5;
double adapt = 0.0;

double raw_delta = delta * eta / 1e-3;

double pot_strength = 0.1;

double ontime = 1.0;
double offtime = 5.0;

double scale = 35;

double wmax = 5.0;
double wmin = 0.0;
double wmaxi = 5.0;
double wmini = 0.0;

double wtmax = 1.0 / 4 * (weight_c - weight_a);
double normalization_factor = (wtmax - weight_a) * (wtmax - (weight_a + weight_c) / 2) * (wtmax - weight_c);

double bgrate = 10.0;

double sleep_duration = 100;
double awake_duration = 200;

string stimfile = "";   // stimulus patterns file
string prefile = "";    // preload patters file
string recfile = "";    // file with receptive fields
string monfile = "";    // patternsto monitor file
string premonfile = ""; // patternsto monitor file

AurynWeight wee = 0.1;
AurynWeight wei = 0.2;
AurynWeight wie = 0.2;
AurynWeight wii = 0.2;

AurynWeight chi = 1.0;
AurynWeight xi = 0.5;

double sparseness = 0.1;
double sparseness_ext = 0.05;
double wext = 0.2;

double simtime = 3600.;

int errcode = 0;

int startup(int ac, char *av[]) {
  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("load", po::value<string>(), "input weight matrix")("prefix", po::value<string>(), "set file prefix")(
        "sleeps", "whether the network should sleep")("save", "save network state at end of sim")("noconsolidation", "switches off consolidation")(
        "noisp", "switches off isp")("noisyweights", "switches noisy initial weights on")("consolidateweights", "initialize weights as consolidated")(
        "quiet", "quiet mode")("alpha", po::value<double>(), "exc input rate")("kappa", po::value<double>(), "hom parameter")(
        "taud", po::value<double>(), "time constant of synaptic depression")("tauf", po::value<double>(), "time constant of synaptic facilitation")(
        "tauh", po::value<double>(), "time constant of homeostasis")("ujump", po::value<double>(), "u jump STP constant")(
        "chi", po::value<double>(), "chi factor - pattern preload strength")("xi", po::value<double>(), "xi factor - stimulation strength")(
        "wext", po::value<double>(), "recurrent weight (ext)")("wee", po::value<double>(), "recurrent weight (wee)")(
        "wei", po::value<double>(), "recurrent weight (wei)")("wii", po::value<double>(), "recurrent weight (wii)")(
        "wie", po::value<double>(), "recurrent weight (wie)")("extsparse", po::value<double>(), "external sparseness")("intsparse", po::value<double>(),
                                                                                                                       "external sparseness")(
        "simtime", po::value<double>(), "simulation time")("ontime", po::value<double>(), "simulation time")("offtime", po::value<double>(), "simulation time")(
        "dir", po::value<string>(), "output dir")("eta", po::value<double>(), "the learning rate")("beta", po::value<double>(), "decay parameter")(
        "potstrength", po::value<double>(), "potential strength parameter")("delta", po::value<double>(), "growth parameter")(
        "weight_a", po::value<double>(), "weight_a")("weight_c", po::value<double>(), "weight_c")("size", po::value<int>(), "simulation size")(
        "seed", po::value<int>(), "random seed ")("stimfile", po::value<string>(), "stimulus file")("prefile", po::value<string>(), "preload file")(
        "recfile", po::value<string>(), "receptive field file")("scale", po::value<double>(), "stimulus strength")(
        "adapt", po::value<double>(), "adaptation jump size for long time constant")("bgrate", po::value<double>(), "background rate of input")(
        "monfile", po::value<string>(), "monitor file")("premonfile", po::value<string>(), "presynaptic monitor file");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    if (vm.count("load")) infilename = vm["load"].as<string>();

    if (vm.count("prefix")) file_prefix = vm["prefix"].as<string>();

    if (vm.count("sleeps")) sleeps = true;

    if (vm.count("save")) save = true;

    if (vm.count("noconsolidation")) consolidation = false;

    if (vm.count("noisp")) isp_active = false;

    if (vm.count("noisyweights")) noisy_initial_weights = true;

    if (vm.count("consolidateweights")) consolidate_initial_weights = true;

    if (vm.count("quiet")) quiet = true;

    if (vm.count("alpha")) alpha = vm["alpha"].as<double>();

    if (vm.count("kappa")) kappa = vm["kappa"].as<double>();

    if (vm.count("taud")) taud = vm["taud"].as<double>();

    if (vm.count("tauf")) tauf = vm["tauf"].as<double>();

    if (vm.count("tauh")) tauh = vm["tauh"].as<double>();

    if (vm.count("ujump")) ujump = vm["ujump"].as<double>();

    if (vm.count("simtime")) simtime = vm["simtime"].as<double>();

    if (vm.count("ontime")) ontime = vm["ontime"].as<double>();

    if (vm.count("offtime")) offtime = vm["offtime"].as<double>();

    if (vm.count("dir")) dir = vm["dir"].as<string>();

    if (vm.count("chi")) chi = vm["chi"].as<double>();

    if (vm.count("xi")) xi = vm["xi"].as<double>();

    if (vm.count("wext")) wext = vm["wext"].as<double>();

    if (vm.count("wee")) wee = vm["wee"].as<double>();

    if (vm.count("wei")) wei = vm["wei"].as<double>();

    if (vm.count("wii")) wii = vm["wii"].as<double>();

    if (vm.count("wie")) wie = vm["wie"].as<double>();

    if (vm.count("extsparse")) sparseness_ext = vm["extsparse"].as<double>();

    if (vm.count("intsparse")) sparseness = vm["intsparse"].as<double>();

    if (vm.count("eta")) eta = vm["eta"].as<double>();

    if (vm.count("beta")) beta = vm["beta"].as<double>();

    if (vm.count("potstrength")) pot_strength = vm["potstrength"].as<double>();

    if (vm.count("delta")) delta = vm["delta"].as<double>();

    if (vm.count("weight_a")) weight_a = vm["weight_a"].as<double>();

    if (vm.count("weight_c")) weight_c = vm["weight_c"].as<double>();

    if (vm.count("size")) size = vm["size"].as<int>();

    if (vm.count("stimfile")) {
      stimfile = vm["stimfile"].as<string>();
      monfile = stimfile;
      premonfile = stimfile;
    }

    if (vm.count("prefile")) prefile = vm["prefile"].as<string>();

    if (vm.count("recfile")) recfile = vm["recfile"].as<string>();

    if (vm.count("scale")) scale = vm["scale"].as<double>();

    if (vm.count("adapt")) adapt = vm["adapt"].as<double>();

    if (vm.count("bgrate")) bgrate = vm["bgrate"].as<double>();

    if (vm.count("monfile")) monfile = vm["monfile"].as<string>();

    if (vm.count("premonfile")) premonfile = vm["premonfile"].as<string>();

    if (vm.count("seed")) seed = vm["seed"].as<int>();
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  raw_delta = delta * eta / 1e-3;
  wtmax = 1.0 / 4 * (weight_c - weight_a);
  normalization_factor = (wtmax - weight_a) * (wtmax - (weight_a + weight_c) / 2) * (wtmax - weight_c);

  auryn_init(ac, av, dir, "sim_sleep", file_prefix);
  sys->set_master_seed(42);

  logger->set_logfile_loglevel(VERBOSE);
  logger->parameter("alpha", alpha);
  logger->parameter("beta", beta);
  logger->parameter("delta", delta);
  logger->parameter("eta", eta);
  logger->parameter("wee", wee);
  logger->parameter("wext", wext);
  logger->parameter("chi", chi);
  logger->parameter("xi", xi);
  logger->parameter("stimfile", stimfile);
  logger->parameter("monfile", monfile);
  logger->parameter("offtime", offtime);
  logger->parameter("ontime", ontime);
  logger->parameter("taud", taud);
  logger->parameter("tauf", tauf);
  logger->parameter("ujump", ujump);
  logger->parameter("tauh", tauh);

  return 0;
}