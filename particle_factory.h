#ifndef PARTICLE_FACTORY_H
#define PARTICLE_FACTORY_H

#include "particle_catalogue.h"
#include "particle.h"
#include <memory>
#include <iostream>

template<typename ParticleType, typename... Args>
std::shared_ptr<ParticleType> create_add_particle(ParticleCatalogue<Particle>& catalogue, Args&&... args)
{
  try
  {
    auto particle = std::make_shared<ParticleType>(std::forward<Args>(args)...);
    catalogue.add_particle(particle);
    return particle;
  }
  catch(const std::exception& e) // If there's an error, eg input momentum out of bounds (too big), then particle is not added to catalogue
  {
    std::cerr<<"Error creating particle: "<<e.what()<<std::endl;
    return nullptr;
  }
}

#endif // PARTICLE_FACTORY_H
