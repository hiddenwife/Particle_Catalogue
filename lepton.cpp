#include "lepton.h"
#include "particle.h"
#include "fourmom.h"
#include "quark.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <numeric>

// Lepton implementation
Lepton::Lepton(double mass, double charge, double px, double py, double pz,
               const std::string& type, bool is_anti, int electron_number, int muon_number, int tau_number)
  : Particle(mass, charge, 0.5, std::sqrt(px*px + py*py + pz*pz + mass*mass), px, py, pz, type, is_anti),
    electron_lepton_number(electron_number),
    muon_lepton_number(muon_number),
    tau_lepton_number(tau_number) {}

Lepton::Lepton(const Lepton& other, bool copy_decay_products) : Particle(other, copy_decay_products)
{
  electron_lepton_number = other.electron_lepton_number;
  muon_lepton_number = other.muon_lepton_number;
  tau_lepton_number = other.tau_lepton_number;
}

Lepton& Lepton::operator=(const Lepton& other)
{
  if(this != &other)
  {
    Particle::operator=(other);
    electron_lepton_number = other.electron_lepton_number;
    muon_lepton_number = other.muon_lepton_number;
    tau_lepton_number = other.tau_lepton_number;
  }
  return *this;
}

// Move constructor
Lepton::Lepton(Lepton&& other) noexcept
  : Particle(std::move(other)), // Invoke the move constructor of the base class
    electron_lepton_number(other.electron_lepton_number),
    muon_lepton_number(other.muon_lepton_number),
    tau_lepton_number(other.tau_lepton_number)
{
  // Clear the moved-from object's data members
  other.electron_lepton_number = 0;
  other.muon_lepton_number = 0;
  other.tau_lepton_number = 0;
}

// Move assignment operator
Lepton& Lepton::operator=(Lepton&& other) noexcept
{
  if(this != &other)
  {
    Particle::operator=(std::move(other)); // Invoke the move assignment operator of the base class
    electron_lepton_number = other.electron_lepton_number;
    muon_lepton_number = other.muon_lepton_number;
    tau_lepton_number = other.tau_lepton_number;

    // Clear the moved-from object's data members
    other.electron_lepton_number = 0;
    other.muon_lepton_number = 0;
    other.tau_lepton_number = 0;
  }
  return *this;
}

void Lepton::print() const
{
  Particle::print(); // Call the base class print function
  std::cout<<"  Electron Lepton Number: "<<electron_lepton_number<<"\n"
           <<"  Muon Lepton Number: "<<muon_lepton_number<<"\n"
           <<"  Tau Lepton Number: "<<tau_lepton_number<<"\n"
           <<"  Antiparticle: "<<(is_antiparticle ? "Yes" : "No")<<"\n";
}

// Electron implementations
Electron::Electron(double px, double py, double pz, std::vector<double> deposits, bool is_anti)
  : Lepton(electron_mass, is_anti ? 1 : -1, px, py, pz, is_anti ? "AntiElectron" : "Electron", is_anti,
           is_anti ? -1 : 1, 0, 0)
{
  if(electron_mass <= 0)
  {
    throw std::invalid_argument("Electron mass must be positive.");
  }
  calorimeter_deposits = std::move(deposits);
  adjust_calorimeter_deposits(); // Ensure calorimeter deposits match electron's energy
}

Electron::Electron(const Electron& other) : Lepton(other),
  calorimeter_deposits(other.calorimeter_deposits) {}

Electron& Electron::operator=(const Electron& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    calorimeter_deposits = other.calorimeter_deposits;
  }
  return *this;
}

Electron::Electron(Electron&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    calorimeter_deposits(std::move(other.calorimeter_deposits)) {} // Move the calorimeter deposits

Electron& Electron::operator=(Electron&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    calorimeter_deposits = std::move(other.calorimeter_deposits); // Move the calorimeter deposits
  }
  return *this;
}

void Electron::adjust_calorimeter_deposits()
{
  double total_deposited_energy = std::accumulate(calorimeter_deposits.begin(), calorimeter_deposits.end(), 0.0);
  if(std::abs(total_deposited_energy - get_e()) > 0.05)
  {
    std::cerr<<"Warning: Total energy deposited in calorimeter does not match the " 
             <<(is_antiparticle ? "AntiElectron" : "Electron")
             <<"'s energy ("<<total_deposited_energy<<" instead of "<<get_e() 
             <<"). Adjusting.\n";

    double energy_scale_factor = get_e() / total_deposited_energy;
    for(auto& deposit : calorimeter_deposits)
    {
      deposit *= energy_scale_factor;
    }
  }
}

void Electron::print() const
{
  Lepton::print(); // Call the base class print function
  std::cout<<"  Calorimeter Deposits: ";
  for(const auto& deposit : calorimeter_deposits)
  {
    std::cout<<std::fixed<<std::setprecision(2)<<deposit<<" ";
  }
  std::cout<<"\n";
}


Muon::Muon(double px, double py, double pz, bool isolated, bool is_anti)
  : Lepton(muon_mass, is_anti ? 1 : -1, px, py, pz, is_anti ? "AntiMuon" : "Muon", is_anti,
           0, is_anti ? -1 : 1, 0), is_isolated(isolated)
{
  if(muon_mass <= 0)
  {
    throw std::invalid_argument("Muon mass must be positive.");
  }
}

Muon::Muon(const Muon& other) : Lepton(other), 
  is_isolated(other.is_isolated) {}

Muon& Muon::operator=(const Muon& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    is_isolated = other.is_isolated;
  }
  return *this;
}

Muon::Muon(Muon&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    is_isolated(other.is_isolated) {} // Move the is_isolated flag

Muon& Muon::operator=(Muon&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    is_isolated = other.is_isolated; // Move the is_isolated flag
  }
  return *this;
}


void Muon::print() const
{
  Lepton::print(); // Call base class print function first
  std::cout<<"  Is Isolated: "<<(is_isolated ? "Yes" : "No")<<"\n";
}

Tau::Tau(double px, double py, double pz, bool is_anti)
  : Lepton(tau_mass, is_anti ? 1 : -1, px, py, pz, is_anti ? "AntiTau" : "Tau", is_anti, 0, 0, is_anti ? -1 : 1)
{
  if(tau_mass <= 0)
  {
    throw std::invalid_argument("Tau mass must be positive.");
  }
}

Tau::Tau(const Tau& other, bool copy_decay_products) : Lepton(other, copy_decay_products),
  decay_type(other.decay_type) {}

Tau& Tau::operator=(const Tau& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    decay_type = other.decay_type;
  }
  return *this;
}

Tau::Tau(Tau&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    decay_type(std::move(other.decay_type)) {} // Move the decay_type

Tau& Tau::operator=(Tau&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    decay_type = std::move(other.decay_type); // Move the decay_type
  }
  return *this;
}

void Tau::print() const
{
  Lepton::print(); // Call base class print function first
  std::cout<<"Decay Type: "<<(decay_type)<<"\n";
  std::cout<<"Decay Products:\n";
  for(const auto& product : decay_products)
  {
    product->print();
  }
}

std::shared_ptr<Particle> Tau::clone() const
{
  auto new_tau = std::make_shared<Tau>(*this);  // Uses copy constructor
  new_tau->clear_decay_products();  // Clear any existing decay products
  for(const auto& dp : decay_products)
  {
    new_tau->add_decay_product(dp->clone());  // Recursively clone decay products
  }
  return new_tau;
}

ElectronNeutrino::ElectronNeutrino(double px, double py, double pz, bool interacted, bool is_anti)
  : Lepton(electron_neutrino_mass, 0, px, py, pz, is_anti ? "AntiElectronNeutrino" : "ElectronNeutrino", is_anti, is_anti ? -1 : 1, 0, 0),
    has_interacted(interacted) {}

ElectronNeutrino::ElectronNeutrino(const ElectronNeutrino& other) : Lepton(other),
  has_interacted(other.has_interacted) {}

ElectronNeutrino& ElectronNeutrino::operator=(const ElectronNeutrino& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    has_interacted = other.has_interacted;
  }
  return *this;
}

ElectronNeutrino::ElectronNeutrino(ElectronNeutrino&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    has_interacted(other.has_interacted) {} // Move the has_interacted flag

ElectronNeutrino& ElectronNeutrino::operator=(ElectronNeutrino&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    has_interacted = other.has_interacted; // Move the has_interacted flag
  }
  return *this;
}

MuonNeutrino::MuonNeutrino(double px, double py, double pz, bool interacted, bool is_anti)
  : Lepton(muon_neutrino_mass, 0, px, py, pz, is_anti ? "AntiMuonNeutrino" : "MuonNeutrino", is_anti,
           0, is_anti ? -1 : 1, 0),
    has_interacted(interacted) {}

MuonNeutrino::MuonNeutrino(const MuonNeutrino& other) : Lepton(other),
  has_interacted(other.has_interacted) {}

MuonNeutrino& MuonNeutrino::operator=(const MuonNeutrino& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    has_interacted = other.has_interacted;
  }
  return *this;
}

MuonNeutrino::MuonNeutrino(MuonNeutrino&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    has_interacted(other.has_interacted) {} // Move the has_interacted flag

MuonNeutrino& MuonNeutrino::operator=(MuonNeutrino&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    has_interacted = other.has_interacted; // Move the has_interacted flag
  }
  return *this;
}

TauNeutrino::TauNeutrino(double px, double py, double pz, bool interacted, bool is_anti)
  : Lepton(tau_neutrino_mass, 0, px, py, pz, is_anti ? "AntiTauNeutrino" : "TauNeutrino", is_anti,
           0, 0, is_anti ? -1 : 1),
    has_interacted(interacted) {}

TauNeutrino::TauNeutrino(const TauNeutrino& other) : Lepton(other),
  has_interacted(other.has_interacted) {}

TauNeutrino& TauNeutrino::operator=(const TauNeutrino& other)
{
  Lepton::operator=(other);
  if(this != &other)
  {
    has_interacted = other.has_interacted;
  }
  return *this;
}

TauNeutrino::TauNeutrino(TauNeutrino&& other) noexcept
  : Lepton(std::move(other)), // Invoke the base class move constructor
    has_interacted(other.has_interacted) {} // Move the has_interacted flag

TauNeutrino& TauNeutrino::operator=(TauNeutrino&& other) noexcept
{
  if(this != &other)
  {
    Lepton::operator=(std::move(other)); // Invoke the base class move assignment operator
    has_interacted = other.has_interacted; // Move the has_interacted flag
  }
  return *this;
}

void ElectronNeutrino::print() const
{
  Lepton::print();
  std::cout<<"  Has Interacted: "<<(has_interacted ? "Yes" : "No")<<"\n";
}

void MuonNeutrino::print() const
{
  Lepton::print();
  std::cout<<"  Has Interacted: "<<(has_interacted ? "Yes" : "No")<<"\n";
}

void TauNeutrino::print() const
{
  Lepton::print();
  std::cout<<"  Has Interacted: "<<(has_interacted ? "Yes" : "No")<<"\n";
}

// Lepton numbers

int Lepton::get_electron_lepton_number() const
{
  return 0;  
}

int Lepton::get_muon_lepton_number() const
{
  return 0;  
}

int Lepton::get_tau_lepton_number() const
{
  return 0;  
}


int Electron::get_electron_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

int ElectronNeutrino::get_electron_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

int Muon::get_muon_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

int MuonNeutrino::get_muon_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

int Tau::get_tau_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

int TauNeutrino::get_tau_lepton_number() const
{
  return is_antiparticle ? -1 : 1;
}

// Decays
void Lepton::decay() {}
void Electron::decay() {}
void Muon::decay() {}
void ElectronNeutrino::decay() {}
void MuonNeutrino::decay() {}
void TauNeutrino::decay() {}

void Tau::decay()
{
  // Assume decay_products is a member variable of Tau
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0); 
  double borrowed_energy = 0;

  if(dis(gen) < 0.33)
  { // Leptonic decay
    decay_type = "Leptonic";
    if(dis(gen) < 0.5)
    { // Muonic decay, equal possibility
      
      auto muon_tau = std::make_shared<Muon>(0, 0, 0, false, is_antiparticle ? true : false); // Initial momenta set to 0
      auto muonAntiNeutrino_tau = std::make_shared<MuonNeutrino>(0, 0, 0, false, is_antiparticle ? false : true);
      auto tauNeutrino_tau = std::make_shared<TauNeutrino>(0, 0, 0, false, is_antiparticle ? true : false);

      // Instead of using a local vector, directly add to this->decay_products
      this->add_decay_product(muon_tau);
      this->add_decay_product(muonAntiNeutrino_tau);
      this->add_decay_product(tauNeutrino_tau);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
    }
    else
    { // Electronic decay, equal possibility

      auto electron_tau = std::make_shared<Electron>(0, 0, 0, std::vector<double>{0.511,0,0,0}, is_antiparticle ? true : false); // Initial momenta set to 0
      auto electronAntiNeutrino_tau = std::make_shared<ElectronNeutrino>(0, 0, 0, false, is_antiparticle ? false : true);
      auto tauNeutrino_tau = std::make_shared<TauNeutrino>(0, 0, 0, false, is_antiparticle ? true : false);

      // Instead of using a local vector, directly add to this->decay_products
      this->add_decay_product(electron_tau);
      this->add_decay_product(electronAntiNeutrino_tau);
      this->add_decay_product(tauNeutrino_tau); 
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);
      electron_tau->adjust_calorimeter_deposits();
    }

  }

  else
  { // Hadronic decay
    decay_type = "Hadronic";
    if(dis(gen) < 0.5)
    { // Up-Down decay

      auto up_tau = std::make_shared<UpQuark>(0, 0, 0, is_antiparticle ? ColourCharge::Red : ColourCharge::AntiRed, is_antiparticle ? false : true); // Initial momenta set to 0
      auto down_tau = std::make_shared<DownQuark>(0, 0, 0, is_antiparticle ? ColourCharge::AntiRed : ColourCharge::Red, is_antiparticle ? true : false);
      auto tauNeutrino_tau = std::make_shared<TauNeutrino>(0, 0, 0, false, is_antiparticle ? true : false);

      this->add_decay_product(up_tau);
      this->add_decay_product(down_tau);
      this->add_decay_product(tauNeutrino_tau);
      distribute_energy_momentum(this->decay_products, this->get_e(), this->get_px(), this->get_py(), this->get_pz(), 0.0);

    }
    else
    { // Up-Strange decay

      auto up_tau = std::make_shared<UpQuark>(0, 0, 0, is_antiparticle ? ColourCharge::Blue : ColourCharge::AntiBlue, is_antiparticle ? false : true); // Initial momenta set to 0
      auto strange_tau = std::make_shared<StrangeQuark>(0, 0, 0, is_antiparticle ? ColourCharge::AntiRed : ColourCharge::Red, is_antiparticle ? true : false);
      auto tauNeutrino_tau = std::make_shared<TauNeutrino>(0, 0, 0, false, is_antiparticle ? true : false);

      this->add_decay_product(up_tau);
      this->add_decay_product(strange_tau);
      this->add_decay_product(tauNeutrino_tau);
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
  if(!(check_charge_conservation(this->decay_products)))
  {
    std::cerr<<"Invalid particle decay: charge conservation violated."<<std::endl;
  }
  if(!(check_invariant_mass(this->decay_products, borrowed_energy)))
  {
    std::cerr<<"Invalid particle decay: invariant mass violated."<<std::endl;
  }   
}

// Clones
std::shared_ptr<Particle> Electron::clone() const
{
  return std::make_shared<Electron>(*this);
}
std::shared_ptr<Particle> Muon::clone() const
{
  return std::make_shared<Muon>(*this);
}
std::shared_ptr<Particle> ElectronNeutrino::clone() const
{
  return std::make_shared<ElectronNeutrino>(*this);
}
std::shared_ptr<Particle> MuonNeutrino::clone() const
{
  return std::make_shared<MuonNeutrino>(*this);
}
std::shared_ptr<Particle> TauNeutrino::clone() const
{
  return std::make_shared<TauNeutrino>(*this);
}
