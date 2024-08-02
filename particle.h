#ifndef PARTICLE_H
#define PARTICLE_H

#include "fourmom.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

class Particle
{
protected:
  std::string particle_type;
  std::unique_ptr<FourMomentum> four_momentum;
  double mass;
  double charge;
  double spin;
  bool is_antiparticle;
  std::vector<std::shared_ptr<Particle>> decay_products;

public:
  Particle(double mass, double charge, double spin, double E, double px, double py, double pz, const std::string& type, bool is_anti);
  Particle(const Particle& other, bool copy_decay_products = true); // Copy constructor
  Particle& operator=(const Particle& other); // Copy assignment operator
  Particle(Particle&& other) noexcept;
  Particle& operator=(Particle&& other) noexcept;
  virtual ~Particle();

  virtual void decay() = 0; // Pure virtual function for decay mechanisms
  virtual void print() const;
  virtual std::shared_ptr<Particle> clone() const = 0;

  double get_mass() const;
  double get_charge() const;
  double get_spin() const;
  std::string get_type() const;
  bool get_is_antiparticle() const;

  double get_e() const;
  double get_px() const;
  double get_py() const;
  double get_pz() const;

  virtual int get_electron_lepton_number() const;
  virtual int get_muon_lepton_number() const;
  virtual int get_tau_lepton_number() const;
  virtual double get_baryon_number() const;

  void set_momentum(double E, double px, double py, double pz);
  std::tuple<double, double, double, double> get_momentum() const;

  void copying_decay_products(const Particle& source);

  int total_decay_products() const;
  void add_decay_product(std::shared_ptr<Particle> product);
  const std::vector<std::shared_ptr<Particle>>& get_decay_products() const;
  void clear_decay_products();

  FourMomentum sum_decay_products_fourmomentum() const;

  bool check_lepton_number_conservation(int initial_electron_number, int initial_muon_number, int initial_tau_number, 
                                   const std::vector<std::shared_ptr<Particle>>& decay_products);
    
  bool check_baryon_number_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products);
  bool check_charge_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products) const;
  void distribute_energy_momentum(std::vector<std::shared_ptr<Particle>>& decay_products, double total_energy, double initial_px,
     double initial_py, double initial_pz, double borrowed_energy);
  bool check_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products, double initial_energy, double initial_px, double initial_py, double initial_pz);
  bool check_invariant_mass(const std::vector<std::shared_ptr<Particle>>& decay_products, double borrowed_energy) const;

};

#endif // PARTICLE_H