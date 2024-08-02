#ifndef PARTICLE_CATALOGUE_H
#define PARTICLE_CATALOGUE_H

#include <unordered_map>
#include <vector>
#include <iostream>
#include <memory>
#include <iomanip>
#include <map>
#include <algorithm>
#include "particle.h" 

// Using the template prevents this file from being split into interface and implementation.
template<typename T>
class ParticleCatalogue
{
private:
  // Map from particle type to a list of particles of that type
  std::unordered_map<std::string, std::vector<std::shared_ptr<T>>> particles_by_type;

public:
  void add_particle(std::shared_ptr<T> particle)
  {
    particles_by_type[particle->get_type()].push_back(particle);
  }

  void remove_particle(const std::string& type, std::shared_ptr<T> particle)
  {
    auto& particles = particles_by_type[type];
    auto new_end = std::remove_if(particles.begin(), particles.end(),
                                  [&particle](const std::shared_ptr<T>& p) { return p == particle; });
    particles.erase(new_end, particles.end());
  }

  std::vector<std::shared_ptr<T>> get_particles_of_type(const std::string& type) const
  {
    auto it = particles_by_type.find(type);
    if(it != particles_by_type.end())
    {
      return it->second;
    }
    return {}; // Return empty vector if not found
  }

  void print_catalogue_by_type(const std::string& type) const
  {
    auto particles = get_particles_of_type(type);
    std::cout<<"Printing "<<particles.size()<<" particles of type "<<type<<" and its decay products:\n";
    for(const auto& particle : particles) {
      particle->print();
      std::cout<<"Total number of decay products for "<<particle->get_type()<<" (including subsequent decays): " 
                   <<particle->total_decay_products()<<"\n\n";
    }
  }

  void number_of_type(const std::string& type) const
  {
    std::cout<<"Number of particles of type "<<type<<": "<<get_particles_of_type(type).size()<<std::endl;
  }

  void total_number() const
  {
    size_t total_particles = 0;
    for (const auto& entry : particles_by_type)
    {
      total_particles += entry.second.size();
    }
    std::cout<<"Total number of particles in catalogue: "<<total_particles<<std::endl;
  }

  void print_all() const
  {
    size_t total_particles = 0;
    size_t decay_particles = 0;
    std::cout<<"Printing all particles in the catalogue:"<<std::endl;
    for(const auto& entry : particles_by_type)
    {
      total_particles += entry.second.size();
      for(const auto& particle : entry.second)
      {
        particle->print();
        std::cout<<"Total number of decay products for "<<particle->get_type()<<" (including subsequent decays): " 
                     <<particle->total_decay_products()<<"\n\n";
        decay_particles += particle->total_decay_products();
      }
    }
    std::cout<<"Total number of base particles printed: "<<total_particles<<"\n";
    std::cout<<"Total number of decay particles printed: "<<decay_particles<<"\n";
  }

  void sum_all() const
  {
    FourMomentum total_momentum(0, 0, 0, 0); // Initialize with zero four-momentum
    FourMomentum totaldecay(0, 0, 0, 0);

    for(const auto& entry : particles_by_type)
    {
      for(const auto& particle : entry.second)
      {
        auto [energy, px, py, pz] = particle->get_momentum();
        // Create a temporary FourMomentum object to sum
        FourMomentum particle_momentum(energy, px, py, pz);
        total_momentum = total_momentum + particle_momentum; // Sum up four-momenta
        totaldecay = totaldecay + particle->sum_decay_products_fourmomentum();
      }
    }
    std::cout<<"Total sum of Four-Momentum of all base particles in the catalogue: ("<<total_momentum.get_e()<<", "<<total_momentum.get_px()
             <<", "<<total_momentum.get_py()<<", "<<total_momentum.get_pz()<<") MeV/c\n";
    double total_invariant_mass = total_momentum.invariant_mass(); // Calculate the total invariant mass
    std::cout<<"Total Invariant Mass of all base particles in the catalogue: "<<total_invariant_mass<<" MeV/c^2\n";
    
    std::cout<<"Total Four-Momentum of all decay particles in the catalogue: ("<<totaldecay.get_e()<<", "<<totaldecay.get_px()
             <<", "<<totaldecay.get_py()<<", "<<totaldecay.get_pz()<<") MeV/c\n";
    double decay_invariant_mass = totaldecay.invariant_mass(); // Calculate the total invariant mass
    std::cout<<"Total Invariant Mass of all decay particles in the catalogue: "<<decay_invariant_mass<<" MeV/c^2\n";
  }


void print_particle_types()
{
  std::cout<<"Available particle types:\n";
  std::cout<<std::left<<std::setw(27)<<"Type"<<std::setw(5)<<"Number"<<std::endl;
  for(const auto& entry : particles_by_type)
  {
    std::cout<<std::left<<std::setw(25)<<entry.first
             <<std::setw(5)<<std::internal<<std::setfill(' ')<<get_particles_of_type(entry.first).size()<<std::endl;
  }
}

std::vector<std::string> get_particle_types() const
{
  std::vector<std::string> types;
  for (const auto& entry : particles_by_type)
  {
    types.push_back(entry.first);
  }
  return types;
}


};

#endif //PARTICLE_CATALOGUE_H
