#include "quark.h"
#include "particle.h"
#include "fourmom.h"
#include <iostream>

Quark::Quark(double mass, double charge, double px, double py, double pz, const std::string& type, ColourCharge colour, bool is_anti, double baryon_number)
  : Particle(mass, charge, 0.5, sqrt(px*px + py*py + pz*pz + mass*mass), px, py, pz, type, is_anti), colour(colour), baryon_number(baryon_number)
{
  check_colour_consistency();
}

UpQuark::UpQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(up_mass, is_anti ? -2.0/3.0 : 2.0/3.0, px, py, pz, is_anti ? "AntiUpQuark" : "UpQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}

DownQuark::DownQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(down_mass, is_anti ? 1.0/3.0 : -1.0/3.0, px, py, pz, is_anti ? "AntiDownQuark" : "DownQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}

CharmQuark::CharmQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(charm_mass, is_anti ? -2.0/3.0 : 2.0/3.0, px, py, pz, is_anti ? "AntiCharmQuark" : "CharmQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}

StrangeQuark::StrangeQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(strange_mass, is_anti ? 1.0/3.0 : -1.0/3.0, px, py, pz, is_anti ? "AntiStrangeQuark" : "StrangeQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}

TopQuark::TopQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(top_mass, is_anti ? -2.0/3.0 : 2.0/3.0, px, py, pz, is_anti ? "AntiTopQuark" : "TopQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}

BottomQuark::BottomQuark(double px, double py, double pz, ColourCharge colour, bool is_anti)
  : Quark(bottom_mass, is_anti ? 1.0/3.0 : -1.0/3.0, px, py, pz, is_anti ? "AntiBottomQuark" : "BottomQuark", colour, is_anti, is_anti ? -1.0/3.0 : 1.0/3.0) {}


// Deep copy functionality

Quark::Quark(const Quark& other)
  : Particle(other), colour(other.colour), baryon_number(other.baryon_number) {}

Quark& Quark::operator=(const Quark& other)
{
  if(this != &other)
  {
    Particle::operator=(other);
    colour = other.colour;
    baryon_number = other.baryon_number;
  }
  return *this;
}

Quark::Quark(Quark&& other) noexcept
  : Particle(std::move(other)), // Invoke the base class move constructor
    colour(other.colour), // Move the colour
    baryon_number(other.baryon_number) {} // Move the baryon number

Quark& Quark::operator=(Quark&& other) noexcept
{
  if(this != &other)
  {
    Particle::operator=(std::move(other)); // Invoke the base class move assignment operator
    colour = other.colour; // Move the colour
    baryon_number = other.baryon_number; // Move the baryon number
  }
  return *this;
}

UpQuark::UpQuark(const UpQuark& other)
  : Quark(other) {}

UpQuark& UpQuark::operator=(const UpQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

UpQuark::UpQuark(UpQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

UpQuark& UpQuark::operator=(UpQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

DownQuark::DownQuark(const DownQuark& other)
  : Quark(other) {}

DownQuark& DownQuark::operator=(const DownQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

DownQuark::DownQuark(DownQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

DownQuark& DownQuark::operator=(DownQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

CharmQuark::CharmQuark(const CharmQuark& other)
  : Quark(other) {}

CharmQuark& CharmQuark::operator=(const CharmQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

CharmQuark::CharmQuark(CharmQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

CharmQuark& CharmQuark::operator=(CharmQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

StrangeQuark::StrangeQuark(const StrangeQuark& other)
  : Quark(other) {}

StrangeQuark& StrangeQuark::operator=(const StrangeQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

StrangeQuark::StrangeQuark(StrangeQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

StrangeQuark& StrangeQuark::operator=(StrangeQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

TopQuark::TopQuark(const TopQuark& other)
  : Quark(other) {}

TopQuark& TopQuark::operator=(const TopQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

TopQuark::TopQuark(TopQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

TopQuark& TopQuark::operator=(TopQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

BottomQuark::BottomQuark(const BottomQuark& other)
  : Quark(other) {}

BottomQuark& BottomQuark::operator=(const BottomQuark& other)
{
  if(this != &other)
  {
    Quark::operator=(other);
  }
  return *this;
}

BottomQuark::BottomQuark(BottomQuark&& other) noexcept
  : Quark(std::move(other)) {} // Invoke the base class move constructor

BottomQuark& BottomQuark::operator=(BottomQuark&& other) noexcept
{
  if(this != &other)
  {
    Quark::operator=(std::move(other)); // Invoke the base class move assignment operator
  }
  return *this;
}

void Quark::check_colour_consistency()
{
  // Lambda for swapping colours to their anticolours and vice versa
  auto swap_colour = [this]()->ColourCharge
  {
    switch (this->colour)
    {
      case ColourCharge::Red: return ColourCharge::AntiRed;
      case ColourCharge::Green: return ColourCharge::AntiGreen;
      case ColourCharge::Blue: return ColourCharge::AntiBlue;
      case ColourCharge::AntiRed: return ColourCharge::Red;
      case ColourCharge::AntiGreen: return ColourCharge::Green;
      case ColourCharge::AntiBlue: return ColourCharge::Blue;
      default: 
        throw std::invalid_argument("Invalid colour charge.");
    }
  };

  // Determine if a swap is needed based on is_anti and the current colour charge
  bool needs_swap = (is_antiparticle && (colour == ColourCharge::Red || colour == ColourCharge::Green || colour == ColourCharge::Blue)) ||
                    (!is_antiparticle && (colour == ColourCharge::AntiRed || colour == ColourCharge::AntiGreen || colour == ColourCharge::AntiBlue));

  if(needs_swap)
  {
    colour = swap_colour(); // Perform the swap
    std::cout<<"Invalid colour assignment for "<<particle_type<<". Swapping "<<(is_antiparticle ? "Colour" : "AntiColour") <<" of "
    <<(is_antiparticle ? "AntiQuark" : "Quark")<<" to its respective "<<(is_antiparticle ? "AntiColour" : "Colour")<<".\n";
  }
}

double Quark::get_baryon_number() const
{
  return is_antiparticle ? -1.0/3.0 : 1.0/3.0;
}


void Quark::print() const
{
  Particle::print(); // Call the base class print function
  std::cout<<"  Colour Charge: "<<colour_charge_to_string(colour)<<"\n"
           <<"  Baron number: "<<(baryon_number)<<"\n"
           <<"  Antiparticle: "<<(is_antiparticle ? "Yes" : "No")<<"\n";
}

// Quark decays
void Quark::decay() {}
void UpQuark::decay() {}
void DownQuark::decay() {}
void CharmQuark::decay() {}
void StrangeQuark::decay() {}
void TopQuark::decay() {}
void BottomQuark::decay() {}


// Changing colour variable to a string name.
std::string colour_charge_to_string(ColourCharge colour)
{
  switch (colour)
  {
    case ColourCharge::Red: return "Red";
    case ColourCharge::Green: return "Green";
    case ColourCharge::Blue: return "Blue";
    case ColourCharge::AntiRed: return "AntiRed";
    case ColourCharge::AntiGreen: return "AntiGreen";
    case ColourCharge::AntiBlue: return "AntiBlue";
    default: return "Unknown Colour";
  }
}

// Clones
std::shared_ptr<Particle> UpQuark::clone() const
{
  return std::make_shared<UpQuark>(*this);
}
std::shared_ptr<Particle> DownQuark::clone() const
{
  return std::make_shared<DownQuark>(*this);
}
std::shared_ptr<Particle> CharmQuark::clone() const
{
  return std::make_shared<CharmQuark>(*this);
}
std::shared_ptr<Particle> StrangeQuark::clone() const
{
  return std::make_shared<StrangeQuark>(*this);
}
std::shared_ptr<Particle> TopQuark::clone() const
{
  return std::make_shared<TopQuark>(*this);
}
std::shared_ptr<Particle> BottomQuark::clone() const
{
  return std::make_shared<BottomQuark>(*this);
}

