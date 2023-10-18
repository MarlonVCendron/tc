/*
 * Copyright 2014-2018 Friedemann Zenke
 *
 * This file is part of Auryn, a simulation package for plastic
 * spiking neural networks.
 *
 * Auryn is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Auryn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
 *
 * If you are using Auryn or parts of it for your work please cite:
 * Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations
 * of spiking neural networks using general-purpose computers.
 * Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
 */

#include "SleepGroup.h"

using namespace auryn;

boost::mt19937 SleepGroup::gen = boost::mt19937();

void SleepGroup::init(AurynFloat sleep_duration_seconds, AurynFloat awake_duration_seconds) {
  sleep_duration = (AurynTime)(sleep_duration_seconds / auryn_timestep);
  awake_duration = (AurynTime)(awake_duration_seconds / auryn_timestep);
  // std::cout << sleep_duration << std::endl;
  current_stage = AWAKE;
  next_stage_change_time = awake_duration;
  next_sleep_awake_time = awake_duration;
  sws_duration = (AurynTime)(sleep_duration / 8);
  rem_duration = (AurynTime)(sleep_duration / 8);
  poisson_rate = 25;
  sws_rate = 1;
  rem_rate = 16;

  auryn::sys->register_spiking_group(this);
  if (evolve_locally()) {


    dist = new boost::exponential_distribution<>(1.0);
    die = new boost::variate_generator<boost::mt19937 &, boost::exponential_distribution<> >(gen, *dist);
    salt = sys->get_seed();
    seed(sys->get_seed());
    neuronId = 0;
    set_rate(poisson_rate);
  }
}

SleepGroup::SleepGroup(NeuronID n, AurynFloat sleep_duration_seconds, AurynFloat awake_duration_seconds) : SpikingGroup(n) {
  init(sleep_duration_seconds, awake_duration_seconds);
}

SleepGroup::~SleepGroup() {
  if (evolve_locally()) {
    delete dist;
    delete die;
  }
}

void SleepGroup::set_rate(AurynDouble rate) {
  lambda = 1.0 / (1.0 / rate - auryn_timestep);
  if (evolve_locally()) {
    if (rate > 0.0) {
      AurynDouble r = (*die)() / lambda;
      neuronId = (NeuronID)(r / auryn_timestep + 0.5);
    } else {
      neuronId = std::numeric_limits<NeuronID>::max();
    }
  }
}

AurynDouble SleepGroup::get_rate() { return lambda; }

void SleepGroup::control_sleep() {
  if (auryn::sys->get_clock() >= next_sleep_awake_time) {
    if (current_stage != AWAKE) {
      current_stage = AWAKE;
      next_stage_change_time = auryn::sys->get_clock() + awake_duration;
      next_sleep_awake_time = auryn::sys->get_clock() + awake_duration;
    } 
  }

  if (auryn::sys->get_clock() >= next_stage_change_time) {
    switch (current_stage) {
    case AWAKE:
      current_stage = SWS;
      next_stage_change_time = auryn::sys->get_clock() + sws_duration;
      next_sleep_awake_time = auryn::sys->get_clock() + sleep_duration;
      break;
    case SWS:
      current_stage = REM;
      next_stage_change_time = auryn::sys->get_clock() + rem_duration;
      break;
    case REM:
      current_stage = SWS;
      next_stage_change_time = auryn::sys->get_clock() + sws_duration;
      break;
    }
    // std::cout << current_stage << std::endl;
  }
}

void SleepGroup::generate_spikes() {
  const AurynDouble desired_period = current_stage == SWS ? 1/sws_rate : 1/rem_rate;
  const AurynDouble omega = 2 * M_PI / (desired_period / auryn_timestep);
  const AurynDouble oscillation_factor = sin(omega * auryn::sys->get_clock());

  while (neuronId < get_rank_size()) {
    push_spike(neuronId);

    AurynDouble modulated_lambda = lambda * (1.01 + oscillation_factor);

    if (modulated_lambda <= 0) {
      neuronId++;
      continue;
    }

    AurynDouble r = (*die)() / modulated_lambda;
    neuronId += (NeuronID)(r / auryn_timestep + 1.5);
  }

  neuronId -= get_rank_size();
}

void SleepGroup::evolve() {
  control_sleep();
  if (current_stage != AWAKE) {
    generate_spikes();
  }
}
void SleepGroup::seed(unsigned int s) {
  std::stringstream oss;
  oss << "SleepGroup:: Seeding with " << s << " and " << salt << " salt";
  auryn::logger->msg(oss.str(), NOTIFICATION);

  gen.seed(s + salt);
}
