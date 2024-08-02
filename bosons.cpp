#include "bosons.h"
#include "quark.h"
#include "lepton.h"
#include "particle.h"
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <iostream>

constexpr double planck_constant = 4.135667696e-21; // Planck constant in MeV·s
constexpr double speed_of_light = 299792458; // Speed of light in m/s
constexpr double eV_to_joules = 1.602176634e-19;

Boson::Boson(double mass, double charge, double spin, double px, double py, double pz, const std::string& type)
  : Particle(mass, charge, spin, std::sqrt(px*px + py*py + pz*pz + mass*mass), px, py, pz, type, false) {}

Boson::Boson(const Boson& other, bool copy_decay_products)
  : Particle(other, copy_decay_products) {}

Boson& Boson::operator=(const Boson& other)
{
  if(this != &other)
  {
    Particle::operator=(other);
  }
  return *this;
}

// Move constructor
Boson::Boson(Boson&& other) noexcept
  : Particle(std::move(other)) {}

// Move assignment operator
Boson& Boson::operator=(Boson&& other) noexcept
{
  if(this != &other)
  {
    Particle::operator=(std::move(other));
  }
  return *this;
}

void Boson::print() const
{
  Particle::print(); // Call the base class print function
}

void Boson::decay() {}

// Photon
Photon::Photon(double px, double py, double pz) : Boson(0, 0, 1, px, py, pz, "Photon") {}

Photon::Photon(const Photon& other)
  : Boson(other) {}

Photon& Photon::operator=(const Photon& other)
{
  if(this != &other)
  {
    Boson::operator=(other);
  }
  return *this;
}

Photon::Photon(Photon&& other) noexcept
  : Boson(std::move(other)) {}

Photon& Photon::operator=(Photon&& other) noexcept
{
  if(this != &other)
  {
    Boson::operator=(std::move(other));
  }
  return *this;
}

void Photon::decay() {}

void Photon::print() const
{
  double energy_MeV = this->get_e(); // Energy in MeV
  double energy_joules = energy_MeV * 1e6 * eV_to_joules; // Convert energy from MeV to Joules

  double frequency = energy_joules / planck_constant; // Frequency in Hz
  double wavelength = speed_of_light / frequency; // Wavelength in meters

  wavelength *= 1e9; // nm
  frequency *= 1e-9; // GHz
  Boson::print();
  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"  Frequency: "<<frequency<<" GHz\n";
  std::cout<<"  Wavelength: "<<wavelength<<" nm\n";
}

std::shared_ptr<Particle> Photon::clone() const
{
  return std::make_shared<Photon>(*this);
}

// WBoson
WBoson::WBoson(int charge, double px, double py, double pz, double borrowed_energy)
  : Boson(W_mass, charge, 1, px, py, pz, charge > 0 ? "W+" : "W-"), borrowed_energy(borrowed_energy) {}

WBoson::WBoson(const WBoson& other, bool copy_decay_products)
  : Boson(other, copy_decay_products), borrowed_energy(other.borrowed_energy), decay_type(other.decay_type) {}

WBoson& WBoson::operator=(const WBoson& other)
{
  if(this != &other)
  {
    Boson::operator=(other);
    borrowed_energy = other.borrowed_energy;
    decay_type = other.decay_type;
  }
  return *this;
}

WBoson::WBoson(WBoson&& other) noexcept
  : Boson(std::move(other)), borrowed_energy(other.borrowed_energy), decay_type(other.decay_type) {}

// Move assignment operator for WBoson
WBoson& WBoson::operator=(WBoson&& other) noexcept
{
  if(this != &other)
  {
    Boson::operator=(std::move(other));
    borrowed_energy = other.borrowed_energy;
    decay_type = other.decay_type;
  }
  return *this;
}

std::shared_ptr<Particle> WBoson::clone() const
{
  auto new_W = std::make_shared<WBoson>(*this);  // Uses copy constructor
  new_W->clear_decay_products();  // Clear any existing decay products
  for(const auto& dp : decay_products)
  {
    new_W->add_decay_product(dp->clone());  // Recursively clone decay products
  }
  return new_W;
}

void WBoson::print() const
{
  if(!(borrowed_energy==0))
  {
    std::cout<<"Virtual WBoson with borrowed energy: "<<borrowed_energy<<" MeV\n";
  }
  Boson::print();
  std::cout<<"Decay Type: "<<(decay_type)<<"\n";
  std::cout<<"Decay Products:\n";
  for(const auto& product : decay_products)
  {
    product->print();
  }
}

constexpr double WBoson::get_W_mass() { return W_mass; }

void WBoson::decay()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0); 
  double borrowed_energy = 0;

  if(dis(gen) < 0.33) // Leptonic decay
  {
    decay_type = "Leptonic";
    if(dis(gen) < 1.0/3.0) // Electronic decay, 1/3 probability
    {

      auto electron_W = std::make_shared<Electron>(0, 0, 0, std::vector<double>{0.511,0,0,0}, charge > 0 ? true : false); // Initial momenta set to 0
      auto electronNeutrino_W = std::make_shared<ElectronNeutrino>(0, 0, 0, false, charge > 0 ? false : true);

      // Instead of using a local vector, directly add to this->decay_products
      this->add_decay_product(electron_W);
      this->add_decay_product(electronNeutrino_W);

      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
      electron_W->adjust_calorimeter_deposits();        

    }
    else if(dis(gen) < 2.0/3.0)
    { // Muonic decay, 1/3 possibility

      auto muon_W = std::make_shared<Muon>(0, 0, 0, false, charge > 0 ? true : false); // Initial momenta set to 0
      auto muonNeutrino_W = std::make_shared<MuonNeutrino>(0, 0, 0, false, charge > 0 ? false : true);

      // Instead of using a local vector, directly add to this->decay_products
      this->add_decay_product(muon_W);
      this->add_decay_product(muonNeutrino_W);

      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
    }
    else
    {
      auto tau_W = std::make_shared<Tau>(0, 0, 0, charge > 0 ? true : false); // Initial momenta set to 0
      auto tauNeutrino_W = std::make_shared<TauNeutrino>(0, 0, 0, false, charge > 0 ? false : true);

      // Instead of using a local vector, directly add to this->decay_products
      
      this->add_decay_product(tauNeutrino_W);
      this->add_decay_product(tau_W);
      
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
      tau_W->decay(); // Tau decays after W decay
    }

  }
  else
  { // Hadronic
    decay_type = "Hadronic";
    if(dis(gen) < 0.5)
    { // Up quark + other quark decay

      auto up_W = std::make_shared<UpQuark>(0, 0, 0, charge > 0 ? ColourCharge::Green : ColourCharge::AntiGreen, charge > 0 ? false : true); // Initial momenta set to 0
      this->add_decay_product(up_W);

      if(dis(gen) < 1.0/3.0)
      { // Down quark
        auto down_W = std::make_shared<DownQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiGreen : ColourCharge::Green, charge > 0 ? true : false);
       this->add_decay_product(down_W);
      }
      else if(dis(gen) < 2.0/3.0)
      { // Strange
        auto strange_W = std::make_shared<StrangeQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiGreen : ColourCharge::Green, charge > 0 ? true : false);
        this->add_decay_product(strange_W);
      }
      else
      { // Bottom
        auto bottom_W = std::make_shared<BottomQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiGreen : ColourCharge::Green, charge > 0 ? true : false);
        this->add_decay_product(bottom_W);
      }

      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);

    }
    else
    { // Charm quark + other quark decay

      auto charm_W = std::make_shared<CharmQuark>(0, 0, 0, charge > 0 ? ColourCharge::Blue : ColourCharge::AntiBlue, charge > 0 ? false : true); // Initial momenta set to 0
      this->add_decay_product(charm_W);

      if(dis(gen) < 1.0/3.0)
      { // Down quark
        auto down_W = std::make_shared<DownQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiBlue : ColourCharge::Blue, charge > 0 ? true : false);
        this->add_decay_product(down_W);
      }
      else if(dis(gen) < 2.0/3.0)
      { // Strange
        auto strange_W = std::make_shared<StrangeQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiBlue : ColourCharge::Blue, charge > 0 ? true : false);
        this->add_decay_product(strange_W);
      }
      else
      { // Bottom
        auto bottom_W = std::make_shared<BottomQuark>(0, 0, 0, charge > 0 ? ColourCharge::AntiBlue : ColourCharge::Blue, charge > 0 ? true : false);
        this->add_decay_product(bottom_W);
      }

      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);

    }
  }

  int initial_electron_number = this->get_electron_lepton_number();
  int initial_muon_number = this->get_muon_lepton_number();
  int initial_tau_number = this->get_tau_lepton_number();

  if(!(check_lepton_number_conservation(initial_electron_number, initial_muon_number, initial_tau_number, this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: lepton number conservation violated."<<std::endl;
  }
  if(!(check_baryon_number_conservation(this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: baryon number conservation violated."<<std::endl;
  }
  if(!(check_charge_conservation(this->decay_products))) {
    std::cerr<<"Invalid particle decay: charge conservation violated."<<std::endl;
  }
  if(!(check_invariant_mass(this->decay_products, borrowed_energy)))
  {
    std::cerr<<"Invalid particle decay: invariant mass violated."<<std::endl;
  }   
}

// ZBoson
ZBoson::ZBoson(double px, double py, double pz, double borrowed_energy)
  : Boson(Z_mass, 0, 1, px, py, pz, "ZBoson"), borrowed_energy(borrowed_energy) {}

void ZBoson::print() const
{
  if(!(borrowed_energy == 0))
  {
    std::cout<<"Virtual ZBoson with borrowed energy: "<<borrowed_energy<<" MeV\n";
  }
  Boson::print();
  std::cout<<"Decay Type: "<<(decay_type)<<"\n";
  std::cout<<"Decay Products:\n";
  for(const auto &product : decay_products)
  {
    product->print();
  }
}

ZBoson::ZBoson(const ZBoson &other, bool copy_decay_products)
  : Boson(other, copy_decay_products), borrowed_energy(other.borrowed_energy), decay_type(other.decay_type) {}

ZBoson &ZBoson::operator=(const ZBoson &other)
{
  if(this != &other)
  {
    Boson::operator=(other);
    borrowed_energy = other.borrowed_energy;
    decay_type = other.decay_type;
  }
  return *this;
}

// Move constructor for ZBoson
ZBoson::ZBoson(ZBoson&& other) noexcept
  : Boson(std::move(other)), borrowed_energy(other.borrowed_energy), decay_type(other.decay_type) {}

// Move assignment operator for ZBoson
ZBoson& ZBoson::operator=(ZBoson&& other) noexcept
{
  if(this != &other)
  {
    Boson::operator=(std::move(other));
    borrowed_energy = other.borrowed_energy;
    decay_type = other.decay_type;
  }
  return *this;
}

std::shared_ptr<Particle> ZBoson::clone() const
{
  auto new_Z = std::make_shared<ZBoson>(*this);  // Uses copy constructor
  new_Z->clear_decay_products();  // Clear any existing decay products
  for(const auto& dp : decay_products)
  {
    new_Z->add_decay_product(dp->clone());  // Recursively clone decay products
  }
  return new_Z;
}

constexpr double ZBoson::get_Z_mass() { return Z_mass; }

void ZBoson::decay()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  double borrowed_energy = 0;

  if(dis(gen) < 1.0 / 3.0)
  { // Leptonic decay with 1/3 probability
    decay_type = "Leptonic";
    if(dis(gen) < 1.0 / 6.0)
    { // Electronic decay, 1/6 probability

      auto electron_Z = std::make_shared<Electron>(0, 0, 0, std::vector<double>{0.511, 0, 0, 0}, false); // Initial momenta set to 0
      auto Antielectron_Z = std::make_shared<Electron>(0, 0, 0, std::vector<double>{0.511, 0, 0, 0}, true);
      this->add_decay_product(electron_Z);
      this->add_decay_product(Antielectron_Z);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
      electron_Z->adjust_calorimeter_deposits();
      Antielectron_Z->adjust_calorimeter_deposits();

    }
    else if(dis(gen) < 1.0 / 3.0)
    { // Muonic decay, 1/6 possibility
      auto muon_Z = std::make_shared<Muon>(0, 0, 0, false, false); // Initial momenta set to 0
      auto Antimuon_Z = std::make_shared<Muon>(0, 0, 0, false, true);

      // Instead of using a local vector, directly add to this->decay_products
      this->add_decay_product(muon_Z);
      this->add_decay_product(Antimuon_Z);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
    }
    else if(dis(gen) < 1.0 / 2.0)
    {
      auto tau_Z = std::make_shared<Tau>(0, 0, 0, false); // Initial momenta set to 0
      auto Antitau_Z = std::make_shared<Tau>(0, 0, 0, true);
      this->add_decay_product(tau_Z);
      this->add_decay_product(Antitau_Z);

      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
      tau_Z->decay(); // Tau decays after W decay
      Antitau_Z->decay();
    }
    else if(dis(gen) < 2.0 / 3.0)
    {
      auto ElectronNeutrino_Z = std::make_shared<ElectronNeutrino>(0, 0, 0, false, false);
      auto AntiElectronNeutrino_Z = std::make_shared<ElectronNeutrino>(0, 0, 0, false, true);
      this->add_decay_product(ElectronNeutrino_Z);
      this->add_decay_product(AntiElectronNeutrino_Z);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
    }
    else if(dis(gen) < 5.0 / 6.0)
    {
      auto MuonNeutrino_Z = std::make_shared<MuonNeutrino>(0, 0, 0, false, false);
      auto AntiMuonNeutrino_Z = std::make_shared<MuonNeutrino>(0, 0, 0, false, true);
      this->add_decay_product(MuonNeutrino_Z);
      this->add_decay_product(AntiMuonNeutrino_Z);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);

    }
    else
    {
      auto TauNeutrino_Z = std::make_shared<TauNeutrino>(0, 0, 0, false, false);
      auto AntiTauNeutrino_Z = std::make_shared<TauNeutrino>(0, 0, 0, false, true);
      this->add_decay_product(TauNeutrino_Z);
      this->add_decay_product(AntiTauNeutrino_Z);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
    }
  }
  else
  { // Hadronic
    decay_type = "Hadronic";
    if (dis(gen) < 1.0 / 5.0)
    { // Up quark
      auto up_Z = std::make_shared<UpQuark>(0, 0, 0, ColourCharge::Green, false);
      auto Antiup_Z = std::make_shared<UpQuark>(0, 0, 0, ColourCharge::AntiGreen, true); // Initial momenta set to 0
      this->add_decay_product(up_Z);
      this->add_decay_product(Antiup_Z);
    }
    else if(dis(gen) < 2.0 / 5.0)
    { // Down quark
      auto down_Z = std::make_shared<DownQuark>(0, 0, 0, ColourCharge::Red, false);
      auto Antidown_Z = std::make_shared<DownQuark>(0, 0, 0, ColourCharge::AntiRed, true);
      this->add_decay_product(down_Z);
      this->add_decay_product(Antidown_Z);
    }
    else if(dis(gen) < 3.0 / 5.0)
    { // Charm
      auto charm_Z = std::make_shared<CharmQuark>(0, 0, 0, ColourCharge::Blue, false); // Initial momenta set to 0
      auto Anticharm_Z = std::make_shared<CharmQuark>(0, 0, 0, ColourCharge::AntiBlue, true);
      this->add_decay_product(charm_Z);
      this->add_decay_product(Anticharm_Z);
    }
    else if(dis(gen) < 4.0 / 5.0)
    { // Strange
      auto strange_Z = std::make_shared<StrangeQuark>(0, 0, 0, ColourCharge::Green, false);
      auto Antistrange_Z = std::make_shared<StrangeQuark>(0, 0, 0, ColourCharge::AntiGreen, true);
      this->add_decay_product(strange_Z);
      this->add_decay_product(Antistrange_Z);
    }
    else
    { // Bottom
      auto bottom_Z = std::make_shared<BottomQuark>(0, 0, 0, ColourCharge::Red, false);
      auto Antibottom_Z = std::make_shared<BottomQuark>(0, 0, 0, ColourCharge::AntiRed, true);
      this->add_decay_product(bottom_Z);
      this->add_decay_product(Antibottom_Z);
    }

    distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
  }

  int initial_electron_number = this->get_electron_lepton_number();
  int initial_muon_number = this->get_muon_lepton_number();
  int initial_tau_number = this->get_tau_lepton_number();

  if(!(check_lepton_number_conservation(initial_electron_number, initial_muon_number, initial_tau_number, this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: lepton number conservation violated."<<std::endl;
  }
  if(!(check_baryon_number_conservation(this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: baryon number conservation violated."<<std::endl;
  }
  if(!(check_charge_conservation(this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: charge conservation violated."<<std::endl;
  }
  if(!(check_invariant_mass(this->decay_products, borrowed_energy)))
  {
    std::cerr<<"Invalid particle decay: invariant mass violated."<<std::endl;
  }
}

// HiggsBoson
HiggsBoson::HiggsBoson(double px, double py, double pz)
  : Boson(higgs_mass, 0, 0, px, py, pz, "HiggsBoson") {}

HiggsBoson::HiggsBoson(const HiggsBoson &other, bool copy_decay_products)
  : Boson(other, copy_decay_products), decay_type(other.decay_type) {}

HiggsBoson &HiggsBoson::operator=(const HiggsBoson &other)
{
  if(this != &other)
  {
    Boson::operator=(other);
    decay_type = other.decay_type;
  }
  return *this;
}

// Move constructor for HiggsBoson
HiggsBoson::HiggsBoson(HiggsBoson&& other) noexcept
  : Boson(std::move(other)), decay_type(other.decay_type) {}

// Move assignment operator for HiggsBoson
HiggsBoson& HiggsBoson::operator=(HiggsBoson&& other) noexcept
{
  if(this != &other)
  {
    Boson::operator=(std::move(other));
    decay_type = other.decay_type;
  }
  return *this;
}

std::shared_ptr<Particle> HiggsBoson::clone() const
{
  auto new_H = std::make_shared<HiggsBoson>(*this);  // Uses copy constructor
  new_H->clear_decay_products();  // Clear any existing decay products
  for (const auto& dp : decay_products)
  {
    new_H->add_decay_product(dp->clone());  // Recursively clone decay products
  }
  return new_H;
}

void HiggsBoson::print() const
{
  Boson::print();
  std::cout<<"Decay Type: "<<(decay_type)<<"\n";
  std::cout<<"Decay Products:\n";
  for(const auto &product : decay_products)
  {
    product->print();
  }
}

void HiggsBoson::decay()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0); 
  double borrowed_energy = 0;

  if(dis(gen) < 1.0/4.0)
  { // Virtual Z-Boson decay
    decay_type = "Virtual ZZ";
    borrowed_energy = (2 * ZBoson::get_Z_mass() - this->get_mass())/2; // Amount of energy borrowed
    auto Z1_H = std::make_shared<ZBoson>(0, 0, 0, borrowed_energy);
    auto Z2_H = std::make_shared<ZBoson>(0, 0, 0, borrowed_energy);
    this->add_decay_product(Z1_H);
    this->add_decay_product(Z2_H);
    distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), borrowed_energy);
    Z1_H->decay();
    Z2_H->decay();
  }
  else if(dis(gen) < 2.0/4.0)
  { // Virtual W-Boson decay
    decay_type = "Virtual W-W+";
    borrowed_energy = (2 * WBoson::get_W_mass() - this->get_mass())/2; // Amount of energy borrowed
    auto W_minus_H = std::make_shared<WBoson>(-1, 0, 0, 0, borrowed_energy);
    auto W_plus_H = std::make_shared<WBoson>(1, 0, 0, 0, borrowed_energy);
    this->add_decay_product(W_minus_H);
    this->add_decay_product(W_plus_H);
    distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), borrowed_energy);
    W_minus_H->decay();
    W_plus_H->decay();
  }
  else if(dis(gen) < 3.0/4.0)
  { // Photon decay
    decay_type = "Photon-Photon";
    auto photon1_H = std::make_shared<Photon>(0, 0, 0);
    auto photon2_H = std::make_shared<Photon>(0, 0, 0);
    this->add_decay_product(photon1_H);
    this->add_decay_product(photon2_H);
    distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
  }
  else
  { // Bottom quark decay
    decay_type = "Hadronic";
    auto bottom_H = std::make_shared<BottomQuark>(0, 0, 0, ColourCharge::Red, false); // Initial momenta set to 0
    auto Antibottom_H = std::make_shared<BottomQuark>(0, 0, 0, ColourCharge::AntiRed, true);
    this->add_decay_product(bottom_H);
    this->add_decay_product(Antibottom_H);
    distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
  }

  int initial_electron_number = this->get_electron_lepton_number();
  int initial_muon_number = this->get_muon_lepton_number();
  int initial_tau_number = this->get_tau_lepton_number();

  if(!(check_lepton_number_conservation(initial_electron_number, initial_muon_number, initial_tau_number, this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: lepton number conservation violated."<<std::endl;
  }
  if(!(check_baryon_number_conservation(this->decay_products))) 
  {
    std::cerr<<"Invalid particle decay: baryon number conservation violated."<<std::endl;
  }
  if(!(check_charge_conservation(this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: charge conservation violated."<<std::endl;
  }   
  if(!(check_invariant_mass(this->decay_products, borrowed_energy)))
  {
    std::cerr<<"Invalid particle decay: invariant mass violated."<<std::endl;
  }   
}

// Gluon
Gluon::Gluon(ColourCharge colour1, ColourCharge colour2, double px, double py, double pz)
  : Boson(0, 0, 1, px, py, pz, "Gluon"), colour1(colour1), colour2(colour2)
{
  check_colour_consistency();
}

Gluon::Gluon(const Gluon& other) 
  : Boson(other), colour1(other.colour1), colour2(other.colour2) {}

// Gluon deep copy assignment operator
Gluon& Gluon::operator=(const Gluon& other)
{
  if(this != &other)
  {
    Boson::operator=(other);
    colour1 = other.colour1;
    colour2 = other.colour2;
  }
  return *this;
}

// Move constructor
Gluon::Gluon(Gluon&& other) noexcept
  : Boson(std::move(other)), colour1(std::move(other.colour1)), colour2(std::move(other.colour2)) {
    other.colour1 = ColourCharge::Neutral;
    other.colour2 = ColourCharge::Neutral;
}

// Move assignment operator
Gluon& Gluon::operator=(Gluon&& other) noexcept
{
  if(this != &other) {
    Boson::operator=(std::move(other));
    colour1 = std::move(other.colour1);
    colour2 = std::move(other.colour2);
    other.colour1 = ColourCharge::Neutral;
    other.colour2 = ColourCharge::Neutral;
  }
  return *this;
}


void Gluon::check_colour_consistency()
{
  bool is_first_colour = colour1 == ColourCharge::Red || colour1 == ColourCharge::Green || colour1 == ColourCharge::Blue;
  bool is_second_anticolour = colour2 == ColourCharge::AntiRed || colour2 == ColourCharge::AntiGreen || colour2 == ColourCharge::AntiBlue;

  if(!(is_first_colour && is_second_anticolour))
  {
    throw std::invalid_argument("Gluons must have one colour and one anticolour.");
  }
}

void Gluon::print() const
{
  Boson::print(); // Print the base class information
  std::cout<<"  Colour 1: "<<colour_charge_to_string(colour1)<<"\n"
           <<"  Colour 2: "<<colour_charge_to_string(colour2)<< std::endl;
}

void Gluon::decay() {}

// Clones
std::shared_ptr<Particle> Gluon::clone() const
{
  return std::make_shared<Gluon>(*this);
}
