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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	
	// Creating a noramal distribution for x y and theta
	normal_distribution <double> dist_x(x,std[0]);
	normal_distribution <double> dist_y(y,std[1]);
	normal_distribution <double> dist_theta(theta,std[2]);
	
	
	num_particles = 20;
	for(int i =0;i<num_particles;i++){
		double sample_x,sample_y,sample_theta;
		
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		
		// Creating a particle
		Particle p = {};
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1.0;
		// Adding partcle to particle list
		particles.push_back(p);
		weights.push_back(1.0);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	// having a 0 centred normal distribution for addition of motion uncertainty
	normal_distribution <double> dist_x(0,std_pos[0]);
	normal_distribution <double> dist_y(0,std_pos[1]);
	normal_distribution <double> dist_theta(0,std_pos[2]);
	
	for(int i = 0;i<= num_particles;i++){
		Particle &pt = particles[i];
		double delta_x,delta_y,delta_theta;
		
		// Change in x y and theta calculation
		if(fabs(yaw_rate) > 0.0001){
			delta_x = (velocity/yaw_rate)*( sin(pt.theta + yaw_rate * delta_t) - sin(pt.theta) );
			delta_y = (velocity/yaw_rate)*( cos(pt.theta) - cos(pt.theta + yaw_rate * delta_t)  );
			delta_theta = yaw_rate * delta_t ;
		}
		else{
			delta_x = velocity * delta_t * cos(pt.theta);
			delta_y = velocity * delta_t * sin(pt.theta);
			delta_theta = 0.0;
		}
		
		// Adding delta change and corresponding noise to x y and theta of a particle
		pt.x += (delta_x + dist_x(gen));
		pt.y += (delta_y + dist_y(gen));
		pt.theta += (delta_theta + dist_theta(gen));
		pt.weight = 1.0;
		weights[i] = 1;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	//-----------------------------------------------------------------------------------------------------
	// getting landmarks in the range of vehicle
	
	
	for(int j=0;j < particles.size();j++){
		
		// to store map landmarks in the range of each particle
		vector<LandmarkObs> landmarks_in_range;
		
		for(int k = 0;k < map_landmarks.landmark_list.size();k++){
			float dist_pt_lmrk = sqrt((map_landmarks.landmark_list[k].x_f - particles[j].x) * (map_landmarks.landmark_list[k].x_f - particles[j].x) + (map_landmarks.landmark_list[k].y_f - particles[j].y) * (map_landmarks.landmark_list[k].y_f - particles[j].y));
			if(dist_pt_lmrk <= sensor_range){
				LandmarkObs temp;
				temp.id = map_landmarks.landmark_list[k].id_i;
				temp.x = map_landmarks.landmark_list[k].x_f;
				temp.y = map_landmarks.landmark_list[k].y_f;
				
				landmarks_in_range.push_back(temp);
			}
		}
		// Translating observation landmark from vehicle co-ordinate system to map co-ordinate system
		// Comparing map landmarks for a particle with observation landmarks
		double x_par,y_par,theta_par,meas_prob = 1;
		x_par = particles[j].x;
		y_par = particles[j].y;
		theta_par = particles[j].theta;
		for(int p = 0;p<observations.size();p++){
			double x_mes,y_mes,x_map,y_map;
			x_mes = observations[p].x;
			y_mes = observations[p].y;
			
			x_map = x_par + (cos(theta_par)*x_mes ) - (sin(theta_par)*y_mes);
			y_map = y_par + (sin(theta_par)*x_mes ) + (cos(theta_par)*y_mes);
			
			// find associated landmark on map corresponding to x_map y_map
			// the one landmark with min distance from observation pt
			double land_idx = 0;
			float dist_land = sensor_range;
			for(int t =0;t< landmarks_in_range.size();t++){
				float td = sqrt((x_map - landmarks_in_range[t].x)*(x_map - landmarks_in_range[t].x) + (y_map - landmarks_in_range[t].y)*(y_map - landmarks_in_range[t].y));
				if(td < dist_land){
					land_idx = t;
					dist_land = td;
				}
			}
			
			// Now we have a observation(x1,y1) near to land mark (xl,yl)
			// will find the probablity to get observation from land mark using 2D gaussian and given sensor deviation in error
			float sig_x = std_landmark[0];
			float sig_y = std_landmark[1];
			//calculate normalization term
			float gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
			// calculate exponent
			float exponent= ((x_map - landmarks_in_range[land_idx].x)*(x_map - landmarks_in_range[land_idx].x))/(2 * sig_x*sig_x) + ((y_map - landmarks_in_range[land_idx].y)*(y_map - landmarks_in_range[land_idx].y))/(2 * sig_y*sig_y);
			// calculate weight using normalization terms and exponent
			float weight= gauss_norm * exp(-exponent);
			
			meas_prob *= weight;
		}
		particles[j].weight = meas_prob;
		weights[j] = meas_prob;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> new_particles;

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> distribution(weights.begin(), weights.end());

    for(int i = 0; i < num_particles; i++){
        Particle p = particles[distribution(gen)];
        new_particles.push_back(p);
    }
    particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
