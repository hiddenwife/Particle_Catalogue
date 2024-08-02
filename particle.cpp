#include "particle.h"
#include "quark.h"
#include "lepton.h"
#include "fourmom.h"
#include <iostream>
#include <iomanip>
#include <random>

Particle::Particle(double mass, double charge, double spin, double E, double px, double py, double pz, const std::string& type, bool is_anti)
  : particle_type(type),
    four_momentum(std::make_unique<FourMomentum>(E, px, py, pz)),
    mass(mass), 
    charge(charge), 
    spin(spin),
    is_antiparticle(is_anti) {}

// Copy constructor
Particle::Particle(const Particle& other, bool copy_decay_products)
  : particle_type(other.particle_type),
    four_momentum(other.four_momentum ? std::make_unique<FourMomentum>(*other.four_momentum) : nullptr),
    mass(other.mass),
    charge(other.charge),
    spin(other.spin),
    is_antiparticle(other.is_antiparticle)
{
  if(copy_decay_products)
  {
    copying_decay_products(other);  // Recursive copy
  }
}

void Particle::copying_decay_products(const Particle& source)
{
  decay_products.clear();  // Clear existing decay products if any
  for(const auto& dp : source.decay_products)
  {
    decay_products.push_back(dp->clone());
    dp->copying_decay_products(*dp);  // Recursive copy
  }
}

// Move constructor
Particle::Particle(Particle&& other) noexcept
  : particle_type(std::move(other.particle_type)),
    four_momentum(std::move(other.four_momentum)),
    mass(other.mass),
    charge(other.charge),
    spin(other.spin),
    is_antiparticle(other.is_antiparticle),
    decay_products(std::move(other.decay_products)) {}

// Move assignment operator
Particle& Particle::operator=(Particle&& other) noexcept
{
  if(this != &other)
  {
    particle_type = std::move(other.particle_type);
    four_momentum = std::move(other.four_momentum);
    mass = other.mass;
    charge = other.charge;
    spin = other.spin;
    is_antiparticle = other.is_antiparticle;
    decay_products = std::move(other.decay_products);
  }
    return *this;
}

// Copy assignment operator
Particle& Particle::operator=(const Particle& other)
{
  if(this != &other)
  {
    particle_type = other.particle_type;
    mass = other.mass;
    charge = other.charge;
    spin = other.spin;
    is_antiparticle = other.is_antiparticle;
    four_momentum = other.four_momentum ? std::make_unique<FourMomentum>(*other.four_momentum) : nullptr;
    decay_products.clear();
    decay_products.reserve(other.decay_products.size());
    for(const auto& dp : other.decay_products)
    {
      decay_products.push_back(dp->clone());;
    }
  }
  return *this;
}

// Virtual destructor
Particle::~Particle() {}

void Particle::print() const
{
  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"Type: "<<particle_type<<"\n"
           <<"  Mass: "<<mass<<" MeV/c^2\n"
           <<"  Invariant mass: "<<four_momentum->invariant_mass()<<" MeV/c^2\n"
           <<"  Charge: "<<charge<<"\n"
           <<"  Spin: "<<spin<<"\n"
           <<"  Four-Momentum: ("<<four_momentum->get_e()<<", "<<four_momentum->get_px()
           <<", "<<four_momentum->get_py()<<", "<<four_momentum->get_pz()<<") MeV/c\n";
}

double Particle::get_mass() const { return mass; }
double Particle::get_charge() const { return charge; }
double Particle::get_spin() const { return spin; }
std::string Particle::get_type() const { return particle_type; }
bool Particle::get_is_antiparticle() const { return is_antiparticle; }

double Particle::get_e() const { return four_momentum->get_e(); }
double Particle::get_px() const { return four_momentum->get_px(); }
double Particle::get_py() const { return four_momentum->get_py(); }
double Particle::get_pz() const { return four_momentum->get_pz(); }

int Particle::get_electron_lepton_number() const { return 0; }
int Particle::get_muon_lepton_number() const { return 0; }
int Particle::get_tau_lepton_number() const { return 0;}

double Particle::get_baryon_number() const { return 0; }
void Particle::add_decay_product(std::shared_ptr<Particle> product) { decay_products.push_back(product); }
const std::vector<std::shared_ptr<Particle>>& Particle::get_decay_products() const { return decay_products; }
void Particle::clear_decay_products() { decay_products.clear(); }

void Particle::set_momentum(double E, double px, double py, double pz)
{
  if(!four_momentum)
  {
    four_momentum = std::make_unique<FourMomentum>(E, px, py, pz);
  }
  else
  {
    four_momentum->set_e(E);
    four_momentum->set_px(px);
    four_momentum->set_py(py);
    four_momentum->set_pz(pz);
  }
}

std::tuple<double, double, double, double> Particle::get_momentum() const
{
  return std::make_tuple(four_momentum->get_e(), four_momentum->get_px(), 
                         four_momentum->get_py(), four_momentum->get_pz());
}

int Particle::total_decay_products() const
{
  int total = decay_products.size(); // Start with direct decay products
  for(const auto& product : decay_products)
  {
    total += product->total_decay_products(); // Recursively count decay products of each decay product
  }
  return total;
}

FourMomentum Particle::sum_decay_products_fourmomentum() const
{
  FourMomentum totalMomentum(0, 0, 0, 0); // Initialize with zero four-momentum for decay products
  for(const auto& decayProduct : decay_products)
  {
    // Add the four-momentum of the decay product
    auto [e, px, py, pz] = decayProduct->get_momentum();
    totalMomentum = totalMomentum + FourMomentum(e, px, py, pz);

    // Recursively add the four-momentum of this decay product's decay products
    totalMomentum = totalMomentum + decayProduct->sum_decay_products_fourmomentum();
  }
  return totalMomentum;
}

void Particle::distribute_energy_momentum(std::vector<std::shared_ptr<Particle>>& decay_products, double total_energy, double initial_px, double initial_py,
                                          double initial_pz, double borrowed_energy)
{ // Borrowed energy for virtual particles (eg in Higgs decay to W-W+ or ZZ)
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for(int iteration = 0; iteration < 5000000; ++iteration)
  { // Limit iterations to prevent infinite loop
    double redistributed_energy = total_energy;
    double remaining_px = initial_px, remaining_py = initial_py, remaining_pz = initial_pz;

    for(size_t i = 0; i < decay_products.size(); ++i)
    {
      double mass = decay_products[i]->get_mass() - borrowed_energy;
      double energy_fraction, px, py, pz;
            
      if(iteration == 0 || dis(gen) < 0.5)
      { // First iteration or 50% chance to explore
        energy_fraction = redistributed_energy / decay_products.size();
        double phi = 2.0 * M_PI * dis(gen);
        double theta = M_PI * dis(gen);
        double p = std::sqrt(std::max(energy_fraction * energy_fraction - mass * mass, 0.0));
        px = p * std::sin(theta) * std::cos(phi);
        py = p * std::sin(theta) * std::sin(phi);
        pz = p * std::cos(theta);
      }
      else
      { // 50% chance to refine based on the best difference
        auto [energy, pX, pY, pZ] = decay_products[i]->get_momentum();
        // Small adjustments towards conservation
        px = pX * (1 + dis(gen) * 0.1 - 0.05);
        py = pY * (1 + dis(gen) * 0.1 - 0.05);
        pz = pZ * (1 + dis(gen) * 0.1 - 0.05);
        energy_fraction = std::sqrt(px * px + py * py + pz * pz + mass * mass);
      }

      // Adjust the last particle's momentum to ensure momentum conservation
      if(i == decay_products.size() - 1)
      {
        px = remaining_px;
        py = remaining_py;
        pz = remaining_pz;
        energy_fraction = std::sqrt(px * px + py * py + pz * pz + mass * mass);
      }
      else
      {
        remaining_px -= px;
        remaining_py -= py;
        remaining_pz -= pz;
      }

      decay_products[i]->set_momentum(energy_fraction, px, py, pz);
      redistributed_energy -= energy_fraction;
    }

    // Check if this iteration is better
    if(check_conservation(decay_products, total_energy, initial_px, initial_py, initial_pz))
    {
      std::cout<<"Energy and momentum conservation for "<<particle_type<<" decay achieved in "<<iteration + 1<<" iterations."<<std::endl;
      return;
    }
  }

  std::cout<<"Failed to achieve conservation for "<<particle_type<<" decay within iteration limit."<<std::endl;
}

bool Particle::check_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products, double initial_energy, double initial_px, double initial_py, double initial_pz)
{
  double total_energy = 0.0, total_px = 0.0, total_py = 0.0, total_pz = 0.0;
  for(const auto& particle : decay_products)
  {
    auto [energy, px, py, pz] = particle->get_momentum(); 
    total_energy += energy;
    total_px += px;
    total_py += py;
    total_pz += pz;
  }

  double tolerance_e = total_energy * 0.001; // Accepting conservation at 0.1% energy tolerance
  return std::abs(total_energy - initial_energy) < tolerance_e &&
         std::abs(total_px - initial_px) < tolerance_e &&
         std::abs(total_py - initial_py) < tolerance_e &&
         std::abs(total_pz - initial_pz) < tolerance_e;
}

bool Particle::check_lepton_number_conservation(int initial_electron_number, int initial_muon_number, int initial_tau_number, 
                                                 const std::vector<std::shared_ptr<Particle>>& decay_products)
{ // Check there is lepton number conservation for the decay products
  int final_electron_number = 0, final_muon_number = 0, final_tau_number = 0;

  // Calculate final lepton numbers
  for(const auto& product : decay_products)
  {
    final_electron_number += product->get_electron_lepton_number();
    final_muon_number += product->get_muon_lepton_number();
    final_tau_number += product->get_tau_lepton_number();
  }

  // Check conservation for each lepton family
  return (initial_electron_number == final_electron_number) &&
         (initial_muon_number == final_muon_number) &&
         (initial_tau_number == final_tau_number);
}

bool Particle::check_baryon_number_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products)
{ // Check there is baryon number conservation for the decay products
  double initial_baryon_number = this->get_baryon_number();
  double final_baryon_number = 0.0;

  // Calculate final baryon number
  for(const auto& product : decay_products)
  {
    final_baryon_number += product->get_baryon_number();
  }

  // Check conservation of baryon number
  return initial_baryon_number == final_baryon_number;
}

bool Particle::check_charge_conservation(const std::vector<std::shared_ptr<Particle>>& decay_products) const
{ // Check there is charge conservation for the decay products
  double initial_charge = this->get_charge(); // Assume getCharge() returns the charge in integer format
  double final_charge = 0;

  // Sum the charge of all decay products
  for(const auto& product : decay_products)
  {
    final_charge += product->get_charge();
  }

  // Check if the total charge is conserved
  return initial_charge == final_charge;
}

bool Particle::check_invariant_mass(const std::vector<std::shared_ptr<Particle>>& decay_products, double borrowed_energy) const
{ // Check the invariant mass of decay products equal rest mass (for non-virtual particles)
  int n = 0;
  for(const auto& product : decay_products)
  {
    double calc_invariant_mass = product->four_momentum->invariant_mass();
    double actual_mass = product->get_mass();
    double tolerance = 1e-2; // Same tolerance as given for energy and momentum conservation
    if((std::abs(calc_invariant_mass - actual_mass) > tolerance) && (borrowed_energy == 0)) // Borrowed energy included to ignore virtual particles
    {
      std::cout<<"Invariant mass of "<<product->get_type()<<" does not equal rest mass for decay of "<<this->get_type()<<"\n";
      n +=1;
    }
  }
  return (n > 0 ? false : true); // All decay products' invariant masses match their actual masses
}
