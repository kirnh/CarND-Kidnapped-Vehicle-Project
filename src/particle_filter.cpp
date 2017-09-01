/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <chrono>

#include "map.h"
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//  x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  // Setting number of particles
  num_particles = 1;

  // Set the first positions of the particles based on estimates from GPS as mean and uncertainty of the GPS as Gaussian noise
  // construct a trivial random generator engine from a time-based seed:
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);

  normal_distribution<double> x_distribution (x, std[0]);
  normal_distribution<double> y_distribution (y, std[1]);
  normal_distribution<double> theta_distribution (theta, std[2]);

  for (int i=0; i<num_particles; i++){
  Particle particle;
  particle.x = x_distribution(generator);
  particle.y = y_distribution(generator);
  particle.theta = theta_distribution(generator);
  particle.weight = 1;

  particles.push_back(particle);

  }

  // Set the flag to indicate the initialization
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Instantiating a new vector to hold the predicted particle states
	vector<Particle> predicted_particles;

  // Looping over all the particles (i.e., num_particles number of times)
	for (int i=0; i< num_particles; i++) {
		Particle particle = particles[i];
		double x = particle.x;
		double y = particle.y;
		double theta = particle.theta;

    // Adding the predicted measurements obtained from previous velocity, yaw_rate and delta_t
	  if (yaw_rate == 0) {
			x = x + velocity*delta_t*cos(theta);
			y = y + velocity*delta_t*sin(theta);
			theta = theta;
		}

		else if (yaw_rate != 0) {
			x = x + (velocity/yaw_rate)*(sin(theta+yaw_rate*delta_t)-sin(theta));
			y = y + (velocity/yaw_rate)*(cos(theta)-cos(theta+yaw_rate*delta_t));
			theta = theta + yaw_rate*delta_t;
		} 

		// Generating random seed to get the random Gaussian noise from the standard deviation values
	  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	  default_random_engine generator (seed);

	  // Adding random noise to introduce uncertainties when measuring state variables
	  normal_distribution<double> x_distribution (x, std_pos[0]);
	  normal_distribution<double> y_distribution (y, std_pos[1]);
	  normal_distribution<double> theta_distribution (theta, std_pos[2]);

	  // Generating random state variables after considering mean values and the standard deviations
	  particle.x = x_distribution(generator);
	  particle.y = y_distribution(generator);
	  particle.theta = theta_distribution(generator);

	  // Push the predicted_particle to the new vector holding all predicted particle states
	  predicted_particles.push_back(particle);

	  // Assign this new vector to represent the particles vector
	  particles = predicted_particles; 
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	// New vector to hold the updated observations list
	vector<LandmarkObs> updated_observations;
  // Looping over all the observations
  for (int i=0; i<observations.size(); i++){
  	LandmarkObs observation = observations[i];
  	// Least distance
  	double least_distance = 9999999999999;
  	// Looping over all the predicted landmark positions
  	for (int j=0; j<predicted.size(); j++){
  		LandmarkObs landmark = predicted[j]; 
  		// Compare the observation to the predicted landmark and associate it with the closest one.
  		// Using dist(x1, y1, x2, y2) function from helper_functions.h to get the Euclidean distance 
  		// between the predicted landmark and the observation measurements
  		double distance = dist(landmark.x, landmark.y, observation.x, observation.y);
  		// store the state if its the least distance
  		if (distance < least_distance){
  			least_distance = distance;
  			observation.id = j;
  		}
  	}
  	updated_observations.push_back(observation);
  }
  observations = updated_observations;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// standard deviations in x and y direction
	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];

  // Looping through all the particles
	for (int i=0; i<num_particles; i++){
		// current Particle in the loop
		Particle particle = particles[i];
		// Ground truth of current particle's state in Map's coordinate system
		ground_truth predicted;
		predicted.x = particle.x;
		predicted.y = particle.y;
		predicted.theta = particle.theta;
		
		// Vector to store the observation in map's coordinate system after transformation from car's system
		vector<LandmarkObs> observations_m;
		// Converting all the observations in car's coordinate system to the map's coordinate system
    for (int j=0; j<observations.size(); j++){
    	LandmarkObs observation = observations[j];
    	LandmarkObs observation_m;

      // observation variables
			double x_o = observation.x;
			double y_o = observation.y;
			// predicted variables
			double x_p = predicted.x;
			double y_p = predicted.y;
			double theta = predicted.theta;	
			// Transformation function
			double x_m = x_p + cos(theta)*x_o + (-sin(theta)*y_o);
			double y_m = y_p + sin(theta)*x_o + cos(theta)*y_o;
			// Observation in map's coordinate system
			observation_m.x = x_m;
			observation_m.y = y_m;

    	// Push this observation to the vector storing all obs in map's system
    	observations_m.push_back(observation_m);
    }
    
    // Using the predicted state variables and the Map to calculate a vector of 
    // predicted landmark positions after translation.
    // Landmarks list from the map
    vector<Map::single_landmark_s> landmark_list = map_landmarks.landmark_list;
    // Predicted landmark list from the predicted state variables
    vector<LandmarkObs> predicted_landmark_list;
   	// Looping over the landmark list to get the predicted landmark list
   	for (int j=0; j<landmark_list.size(); j++){
   		// current landmark from map
   		Map::single_landmark_s landmark = landmark_list[j];
      double landmark_x = landmark.x_f;
      double landmark_y = landmark.y_f;
      int landmark_id = landmark.id_i;
      // Predicted landmark observation
      LandmarkObs predicted_landmark;
      predicted_landmark.x = predicted.x + landmark_x;
      predicted_landmark.y = predicted.y + landmark_y;
      predicted_landmark.id = landmark_id;
      // Push this to the list storing all the predicted landmark obs
      predicted_landmark_list.push_back(predicted_landmark);
   	}

    // Data association of this measured observation (in map's coordinate system) with the 
    // predicted landmark positions (obtained from the current state of the particle and the map)
    dataAssociation(predicted_landmark_list, observations_m);

    // Update weights for each particle depending on P(Z|X) where Z is the observation and X
    // is the particle's position
    // Looping through all the observations and capturing its multivariate gaussian probability
    double weight = 1; 
    for (int j=0; j<observations_m.size(); j++){
    	// current observation
    	LandmarkObs observation = observations[j];
    	int landmark_id = observation.id;
    	// associated predicted landmark
    	LandmarkObs predicted_landmark = predicted_landmark_list[landmark_id];
    	// get variables for better readability
    	double x_obs = observation.x;
    	double y_obs = observation.y;
    	double mu_x = predicted_landmark.x;
    	double mu_y =  predicted_landmark.y;
      // calculate normalization term
			double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
			// calculate exponent
			double exponent = (pow(x_obs - mu_x,2))/(2 * sig_x * sig_x) + (pow(y_obs - mu_y,2))/(2 * sig_y * sig_y);
      // calculate probability using normalization terms and exponent
      double probability = gauss_norm * exp(-exponent);
      // the final weight is the product of the probabilities of all the observations
      weight = weight * probability;
    }

    // Update the weight of the particle
    particle.weight = weight;

    // Update the weights vector of the particle filter as well
    weights.push_back(weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
  // Discrete distribution of particle weights used to resample
  random_device rd;
  mt19937 gen(rd());
  discrete_distribution<> d(weights.begin(), weights.end());

  // A vector of particles to store the new particels while resampling is undergoing
  vector<Particle> particles_new;
  // Looping over num_particles number of times and resampling from particles vector according to weights
  // from the weights vector 
  for(int n=0; n<num_particles; n++) {
      Particle particle_resampled = particles[d(gen)];
      particles_new.push_back(particle_resampled);
  }

  // The new resampled particles now represents the particles
  particles = particles_new;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
