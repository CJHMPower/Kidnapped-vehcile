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
#include "helper_functions.h"


using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
        num_particles=8;
        random_device r;
        default_random_engine gen(r());
        normal_distribution<> x_axis(x, std[0]);
        normal_distribution<> y_axis(y, std[1]);
        normal_distribution<> theta_axis(theta, std[2]);

        for(int i=0;i<num_particles;i++){
           Particle * part=(Particle*)malloc(sizeof(Particle));
           part->id=i+1;
           part->x=x_axis(gen);
           part->y=y_axis(gen);
           part->theta=theta_axis(gen);
           part->weight=1.0;
           weights.push_back(part->weight);
           part->associations={};
           part->sense_x={};
           part->sense_y={};
           particles.push_back(*part);

        }
        is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        random_device r;
        default_random_engine gen(r());
        normal_distribution<> x_noise(0, std_pos[0]);
        normal_distribution<> y_noise(0, std_pos[1]);
        normal_distribution<> theta_noise(0, std_pos[2]);
        for(int i=0;i<num_particles;i++){
           double new_theta=particles[i].theta+yaw_rate*delta_t;
           if(fabs(yaw_rate)>0.000001){
           particles[i].x+=velocity/yaw_rate*(sin(new_theta)-sin(particles[i].theta));
           particles[i].y+=velocity/yaw_rate*(cos(particles[i].theta)-cos(new_theta));
           particles[i].theta+=yaw_rate*delta_t;
           
          
           }
           else{ 
            particles[i].x+=velocity*delta_t*cos(particles[i].theta);
            particles[i].y+=velocity*delta_t*sin(particles[i].theta);
           }
           particles[i].x+=x_noise(gen);
           particles[i].y+=y_noise(gen);
           particles[i].theta+=theta_noise(gen);
        }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        for(int i=0;i<observations.size();i++){
           double min_dist=pow(2,31);
           int min_index=0;
           for(int j=0;j<predicted.size();j++){
                double distance=dist(observations[i].x, observations[i].y,predicted[j].x,
                predicted [j].y);
                if(distance<min_dist){
                   min_dist=distance;
                   min_index=j;
                }
           }
           observations[i].id=min_index;
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
        
       
       
        for(int i=0;i<num_particles;i++){
            vector<LandmarkObs> predicted={};

         for (auto l = map_landmarks.landmark_list.begin();
         l != map_landmarks.landmark_list.end(); ++l) {
            if (dist(l->x_f, l->y_f, particles[i].x, particles[i].y) < sensor_range) {
            LandmarkObs *obs=(LandmarkObs*)malloc(sizeof(LandmarkObs));
            obs->x = l->x_f;
            obs->y = l->y_f;
            obs->id = l->id_i;
            predicted.push_back(*obs);
            }
         }
            double cos_theta=cos(particles[i].theta);
            double sin_theta=sin(particles[i].theta);
            
            vector<LandmarkObs> transform_observations(observations.size());
            for(int j=0;j<observations.size();j++){
            double xm=particles[i].x+cos_theta*observations[j].x-sin_theta*observations[j].y;
            double ym=particles[i].y+sin_theta*observations[j].x+cos_theta*observations[j].y;
            LandmarkObs * transform=(LandmarkObs*)malloc(sizeof(LandmarkObs));
            transform->id=observations[j].id;
            transform->x=xm;
            transform->y=ym;
            transform_observations[j]=*transform;
            }
            dataAssociation(predicted, transform_observations);
            long double prob=1.0;
            

            for(int j=0;j<transform_observations.size();j++){
               
                double coef=1.0/(2*M_PI*std_landmark[0]*std_landmark[1]);
                double term=pow(transform_observations[j].x-predicted[transform_observations[j].id].x,2)/pow(std_landmark[0],2)+pow(transform_observations[j].y-predicted[transform_observations[j].id].y,2)/pow(std_landmark[1],2);
               
                double weight=coef*exp(-0.5*term);
                prob*=weight;
                
            }
            weights[i]=prob;
            particles[i].weight=prob;
            
        }
        

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
       
       default_random_engine gen;
       discrete_distribution<> d(weights.begin(), weights.end());
       vector<Particle> new_particles;
     for (size_t i = 0; i < particles.size(); ++i) {
    const Particle &src = particles[d(gen)];
    new_particles.push_back(src);
     }
      particles.clear();
      particles.insert(particles.end(), new_particles.begin(), new_particles.end());
       /*uniform_int_distribution<int> uindex(0, num_particles-1);
      
       for(int i=0;i<particles.size();i++){
          weight_sum+=particles[i].weight;
          if(particles[i].weight>wmax)
             wmax=particles[i].weight;
       }
       cout<<"weight sum: "<<weight_sum<<" "<<wmax<<endl;
      
       for(int i=0;i<particles.size();i++){
          particles[i].weight=particles[i].weight/weight_sum; 
          weights[i]=particles[i].weight;
       }
       wmax=wmax/weight_sum;
       
       uniform_real_distribution<double> weight_distribution(0, 2*wmax);
       vector<Particle> new_particles={};
       double beta=0;
       int index=uindex(gen);
       for(int i=0;i<particles.size();i++){
           beta+=weight_distribution(gen);
           while(beta>particles[index].weight){
               beta-=particles[index].weight;
               index=(index+1)%num_particles;
           }
           new_particles.push_back(particles[index]);
       }
       particles=new_particles;
       return;*/
       
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
