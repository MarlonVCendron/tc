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

#ifndef SLEEPGROUP_H_
#define SLEEPGROUP_H_

#include "auryn.h"

#include <boost/random/exponential_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

namespace auryn {

enum SleepStage { AWAKE, SWS, REM };

/*! \brief A SpikingGroup that creates poissonian spikes with a given rate.
 *
 * This is the standard Poisson spike generator of Auryn. It implements a
 * group of given size of Poisson neurons all firing at the same rate.
 * The implementation is very efficient if the rate is constant throughout.
 *
 * The random number generator will be seeded identically every time. Use
 * the seed function to seed it randomly if needed. Note that all SleepGroups
 * in a simulation share the same random number generator. Therefore it
 * sufficed to seed one of them.
 */
class SleepGroup : public SpikingGroup {
private:
  AurynTime *clk;
  AurynDouble lambda;
  static boost::mt19937 gen;
  boost::exponential_distribution<> *dist;
  boost::variate_generator<boost::mt19937 &, boost::exponential_distribution<> > *die;

  unsigned int salt;

  AurynTime next_stage_change_time;
  AurynTime next_sleep_awake_time;
  SleepStage current_stage;
  AurynTime rem_duration;
  AurynTime sws_duration;
  AurynTime sleep_duration;
  AurynTime awake_duration;

  void init(AurynFloat sleep_duration_seconds, AurynFloat awake_duration_seconds);

protected:
  NeuronID neuronId;

public:
  /*! Standard constructor.
   * @param n is the size of the SpikingGroup, i.e. the number of Poisson neurons.
   * @param rate is the mean firing rate of all poisson neurons in the group.
   */
  SleepGroup(NeuronID n, AurynFloat sleep_duration_seconds, AurynFloat awake_duration_seconds);
  virtual ~SleepGroup();
  virtual void evolve();

  void control_sleep();
  void generate_spikes();
  void set_rate(AurynDouble rate);
  AurynDouble get_rate();
  void seed(unsigned int s);

  float poisson_rate;
  float sws_rate;
  float rem_rate;

};

} // namespace auryn

#endif /*SLEEPGROUP_H_*/
