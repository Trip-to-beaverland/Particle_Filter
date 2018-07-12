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
    if (is_initialized){
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
    /*std::vector<double> distances;
    std::vector<double> min_dist;
    distances.resize(predicted.size());
    for (int j = 0; j<predicted.size(); j++){
        for(int i = 0; i<predicted.size(); i++){
            double x1 = observations[0];
            double y1 = observations[1];
            double x2 = predicted[0];
            double y2 = predicted[1];
            distances[i] = dist(x1,y1,x2,y2);
        }
        //double minimum_dist = 0.0;
        min_dist[j] = std::min_element(distances.begin(), distances.end());
        int map_id =
    }*/
    unsigned int n_obs = observations.size();
    unsigned int n_preds = predicted.size();

    for (unsigned int i = 0; i < n_obs; i++){
        //double minimum_dist = 1000000.0;
        double minimum_dist = numeric_limits<double>::max();
        int obs_id = -1;

        for (unsigned int j = 0; j < n_preds; j++ ){
                double xdist = observations[i].x - predicted[j].x;
                double ydist = observations[i].y - predicted[j].y;
                double distances = xdist*xdist + ydist*ydist;
                //double distances = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

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
	//double gauss = 1.0;
	//2. Iterate over all the particles for each observation to transform coordinates
	//for (unsigned int i = 0; i<num_particles; i++){
	for (int i = 0; i<num_particles; i++){
	    double x_p = particles[i].x;
        double y_p = particles[i].y;
        double theta = particles[i].theta;
        double sensor_range2 = sensor_range*sensor_range;
        vector<LandmarkObs> RangedLandmarks;
            for (unsigned int k = 0; k<map_landmarks.landmark_list.size(); k++){
                float lmrk_x = map_landmarks.landmark_list[k].x_f;
                float lmrk_y = map_landmarks.landmark_list[k].y_f;
                int lmrk_id = map_landmarks.landmark_list[k].id_i;
                double lmrk_dist = pow(x_p-lmrk_x,2)+pow(y_p-lmrk_y,2);
                if (lmrk_dist<=sensor_range2){
                    //landmark_dist[k] = sqrt(pow(x_m-lmrk_x,2)+pow(y_m-lmrk_y,2));
                    RangedLandmarks.push_back(LandmarkObs{lmrk_id, lmrk_x, lmrk_y});
                }
            }
        //Transform map observations
        vector<LandmarkObs> MappedObservations;
        for(unsigned int j = 0; j<observations.size(); j++){
            double x_c = observations[j].x;
            double y_c = observations[j].y;

            double x_m = x_p + x_c*cos(theta) - y_c*sin(theta);
            double y_m = y_p + x_c*sin(theta) + y_c*cos(theta);
            MappedObservations.push_back(LandmarkObs{observations[j].id, x_m, y_m});
        }

        dataAssociation(RangedLandmarks, MappedObservations);

        //Reset weights
        particles[i].weight = 1.0;

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
            //Calculating weights
            double diff_x = obs_x - LM_X;
            double diff_y = obs_y - LM_Y;
            double gauss = multiplier*exp(-(diff_x*diff_x/x_den + diff_y*diff_y/y_den));
            if (gauss==0){
                particles[i].weight *= 0.0001;
            }else{
                particles[i].weight *= gauss;
            }
        }
    //3. Iterate over all the particles for each landmark to find the nearest landmark

    //4. Associate all the points with their nearest landmarks; calculate Gaussian

            //int min_index = distance(landmark_dist.begin(),min_element(landmark_dist.begin(),landmark_dist.end()));



            //gauss = gauss*multiplier*exp(-1*((pow(diff_x,2)/x_den)+(pow(diff_y,2)/y_den)));
    //5. Update weights;
    }
            //particles[i].weight = gauss;
            //weights[i] = particles[i].weight;
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<double> weights;
	double maxWeight = numeric_limits<double>::min();
	for(int i = 0; i < num_particles; i++) {
            weights.push_back(particles[i].weight);
    if ( particles[i].weight > maxWeight ) {
        maxWeight = particles[i].weight;
        }
    }

  // Creating distributions.
    uniform_real_distribution<double> distDouble(0.0, maxWeight);
    uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generating index.
    int index = distInt(gen);

    double beta = 0.0;

  // the wheel
    vector<Particle> resampledParticles;
    for(int i = 0; i < num_particles; i++) {
            beta += distDouble(gen) * 2.0;
            while( beta > weights[index]) {
                    beta -= weights[index];
                    index = (index + 1) % num_particles;
                    }
            resampledParticles.push_back(particles[index]);
            }
    particles = resampledParticles;
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
