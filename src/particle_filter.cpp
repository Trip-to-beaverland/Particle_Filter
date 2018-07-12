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

//0. Include random value
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    if (is_initialized==true){
        return;
    }
	//1. Initialize number of particles
	num_particles = 100;

	//2. Extract standard deviations
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	//3. Creating normal distribution
	normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);

    //4. Initialize particles
    for (int i = 0; i<num_particles; i++){
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        particles.push_back(particle);
    }

    //5. Initialize filter
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    //1. Extract standard deviations
	double std_pos_x = std_pos[0];
	double std_pos_y = std_pos[1];
	double std_pos_theta = std_pos[2];

	//2. Creating normal distribution
	normal_distribution<double> dist_x(0, std_pos_x);
    normal_distribution<double> dist_y(0, std_pos_y);
    normal_distribution<double> dist_theta(0, std_pos_theta);

    //3. Adding measurements - just by using our bicycle motion model
    for (int i = 0; i<num_particles; i++){
        double theta = particles[i].theta;
        if (fabs(yaw_rate)<0.0001){
            particles[i].x += velocity*delta_t*cos(theta);
            particles[i].y += velocity*delta_t*sin(theta);
    }
        else{
            particles[i].x += (velocity/yaw_rate)*(sin(theta+yaw_rate*delta_t)-sin(theta));
            particles[i].y += (velocity/yaw_rate)*(cos(theta)-cos(theta+yaw_rate*delta_t));
            particles[i].theta += yaw_rate*delta_t;
    }
    //4. Adding noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
    unsigned int n_obs = observations.size();
    unsigned int n_preds = predicted.size();

    for (unsigned int i = 0; i < n_obs; i++){
        double minimum_dist = 1000000.0;
        int obs_id = -1;

        for (unsigned int j = 0; j < n_preds; j++ ){
                double distances = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

                if (distances < minimum_dist ){
                        minimum_dist = distances;
                        obs_id = predicted[j].id;
                        }
        }
    observations[i].id = obs_id;
    }
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

	//1. Pre-calculate Gaussian elements
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];
	double multiplier = 1 / (2 * M_PI * std_x * std_y);
	double x_den = 2 * pow(std_x,2);
	double y_den = 2 * pow(std_y,2);

	//2. Iterate over all the particles for each observation
	for (int i = 0; i<num_particles; i++){
	    double x_p = particles[i].x;
        double y_p = particles[i].y;
        double theta = particles[i].theta;
        double sensor_range2 = sensor_range*sensor_range;

    //3. Check, whether landmark is under sensor range, if truth - store it
        vector<LandmarkObs> RangedLandmarks;
        for (unsigned int k = 0; k<map_landmarks.landmark_list.size(); k++){
            float lmrk_x = map_landmarks.landmark_list[k].x_f;
            float lmrk_y = map_landmarks.landmark_list[k].y_f;
            int lmrk_id = map_landmarks.landmark_list[k].id_i;
            double lmrk_dist = pow(x_p-lmrk_x,2)+pow(y_p-lmrk_y,2);//it should be faster here, than dist.
            if (lmrk_dist<=sensor_range2){
                    RangedLandmarks.push_back(LandmarkObs{lmrk_id, lmrk_x, lmrk_y});
                }
            }
    //4. Transform map observations from car to map reference frame
        vector<LandmarkObs> MappedObservations;
        for(unsigned int j = 0; j<observations.size(); j++){
            double x_m = x_p + observations[j].x*cos(theta) - observations[j].y*sin(theta);
            double y_m = y_p + observations[j].x*sin(theta) + observations[j].y*cos(theta);
            MappedObservations.push_back(LandmarkObs{observations[j].id, x_m, y_m});
        }
    //5. Finally use dataAssociation method
        dataAssociation(RangedLandmarks, MappedObservations);

    //6. Reset weights of particles
        particles[i].weight = 1.0;

    //7. Run over all the observations, transformed to map reference frame.
        for (unsigned int l=0; l<MappedObservations.size(); l++){
            double obs_x = MappedObservations[l].x;
            double obs_y = MappedObservations[l].y;
            int LM_id = MappedObservations[l].id;
            double LM_X;
            double LM_Y;
            unsigned int m = 0;
            unsigned int n_lmrks = RangedLandmarks.size();
            bool found = false;
            while(!found && m<n_lmrks){
                if(RangedLandmarks[m].id==LM_id){
                    found = true;
                    LM_X = RangedLandmarks[m].x;
                    LM_Y = RangedLandmarks[m].y;
                }
                m++;
            }
    //8. Calculating weights
            double diff_x = obs_x - LM_X;
            double diff_y = obs_y - LM_Y;
            double gauss = multiplier*exp(-(diff_x*diff_x/x_den + diff_y*diff_y/y_den));
            if (gauss==0){
                particles[i].weight *= 0.0001;
            }else{
                particles[i].weight *= gauss;
            }
        }
    }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//Coding in C++ Sebastian's python class on wheels
	vector<double> new_weights;
	double max_weight = 0.0001;
	for(int i = 0; i < num_particles; i++) {
            new_weights.push_back(particles[i].weight);
        if ( particles[i].weight > max_weight ) {
        max_weight = particles[i].weight;
        }
    }

    uniform_real_distribution<double> distDouble(0.0, max_weight);
    uniform_int_distribution<int> distInt(0, num_particles - 1);
    int index = distInt(gen);

    double beta = 0.0;

    //The wheel itself
    vector<Particle> NewParticles;
    for(int i = 0; i < num_particles; i++) {
            beta += distDouble(gen) * 2.0;
            while( beta > new_weights[index]) {
                    beta -= new_weights[index];
                    index = (index + 1) % num_particles;
                    }
            NewParticles.push_back(particles[index]);
            }
    particles = NewParticles;
}


Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Until I added these 3 rows of code, the simulator ran slowly and returned me the message:
    //"You ran out of time" with the same errors for x~0.1, y~0.09 and yaw~0.004.
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
