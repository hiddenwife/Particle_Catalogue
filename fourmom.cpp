#include "fourmom.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

FourMomentum::FourMomentum(double E, double px, double py, double pz)
{
  const double max_mom = 1e10;
  if(std::abs(px) > max_mom || std::abs(py) > max_mom || std::abs(pz) > max_mom)
  {
    throw std::invalid_argument("Momentum component is out of the allowed range."); // Preventing the addition of the particle to the catalogue.
  }
  this->E = E;
  this->px = px;
  this->py = py;
  this->pz = pz;
}

// Copy constructor
FourMomentum::FourMomentum(const FourMomentum& other)
  : E(other.E), px(other.px), py(other.py), pz(other.pz) {}

// Move constructor
FourMomentum::FourMomentum(FourMomentum&& other) noexcept
  : E(other.E), px(other.px), py(other.py), pz(other.pz)
{
  other.E = 0;
  other.px = 0;
  other.py = 0;
  other.pz = 0;
}

// Move assignment operator
FourMomentum& FourMomentum::operator=(FourMomentum&& other) noexcept
{
  if(this != &other)
  {
    E = other.E;
    px = other.px;
    py = other.py;
    pz = other.pz;

    other.E = 0;
    other.px = 0;
    other.py = 0;
    other.pz = 0;
  }
  return *this;
}

// Destructor
FourMomentum::~FourMomentum() {}

void FourMomentum::set_e(double E)
{
  this->E = E > 0 ? E : 0;
}

void FourMomentum::set_px(double px)
{
  this->px = px;
}

void FourMomentum::set_py(double py)
{
  this->py = py;
}

void FourMomentum::set_pz(double pz)
{
  this->pz = pz;
}

double FourMomentum::get_e() const
{
  return E;
}

double FourMomentum::get_px() const
{
  return px;
}

double FourMomentum::get_py() const
{
  return py;
}

double FourMomentum::get_pz() const
{
  return pz;
}

FourMomentum operator+(const FourMomentum& lhs, const FourMomentum& rhs)
{
  return FourMomentum(lhs.get_e() + rhs.get_e(), lhs.get_px() + rhs.get_px(), lhs.get_py() + rhs.get_py(), lhs.get_pz() + rhs.get_pz());
}

FourMomentum operator-(const FourMomentum& lhs, const FourMomentum& rhs)
{
  return FourMomentum(lhs.get_e() - rhs.get_e(), lhs.get_px() - rhs.get_px(), lhs.get_py() - rhs.get_py(), lhs.get_pz() - rhs.get_pz());
}

double dot_product(const FourMomentum& lhs, const FourMomentum& rhs)
{
  // Metric signature (+, -, -, -) 
  return lhs.get_e() * rhs.get_e() - (lhs.get_px() * rhs.get_px() + lhs.get_py() * rhs.get_py() + lhs.get_pz() * rhs.get_pz());
}

double FourMomentum::invariant_mass() const
{
  double mass_squared = dot_product(*this, *this);
  return mass_squared >= 0 ? std::sqrt(mass_squared) : 0;
}
