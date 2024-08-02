#ifndef BOSON_H
#define BOSON_H

#include "particle.h"
#include "fourmom.h"
#include "quark.h"
#include <string>

class Boson : public Particle
{
public:
  Boson(double mass, double charge, double spin, double px, double py, double pz, const std::string& type);
  Boson(const Boson& other, bool copy_decay_products = true); // Copy constructor
  Boson(Boson&& other) noexcept; // Move constructor
  Boson& operator=(const Boson& other); // Copy assignment operator
  Boson& operator=(Boson&& other) noexcept; // Move assignment operator
  virtual ~Boson() = default;

  virtual void decay() override = 0;
  virtual void print() const override;

};

class Photon : public Boson
{
public:
  Photon(double px, double py, double pz);
  Photon(const Photon& other); // Copy constructor declaration
  Photon& operator=(const Photon& other); // Copy assignment operator declaration
  Photon(Photon&& other) noexcept;
  Photon& operator=(Photon&& other) noexcept;
  virtual ~Photon() = default;

  void decay() override;
  void print() const override;

  std::shared_ptr<Particle> clone() const override;
};

class WBoson : public Boson
{
private:
  static constexpr double W_mass = 80377;
  double borrowed_energy;
  std::string decay_type;

public:
  WBoson(int charge, double px, double py, double pz, double borrowed_energy = 0);
  WBoson(const WBoson& other, bool copy_decay_products = false); // Copy constructor declaration
  WBoson& operator=(const WBoson& other); // Copy assignment operator declaration
  WBoson(WBoson&& other) noexcept;
  WBoson& operator=(WBoson&& other) noexcept;
  virtual ~WBoson() = default;

  void decay() override;
  void print() const override;
  static constexpr double get_W_mass();

  std::shared_ptr<Particle> clone() const override;
};

class ZBoson : public Boson
{
private:
  static constexpr double Z_mass = 91187.6;
  double borrowed_energy;
  std::string decay_type;

public:
  ZBoson(double px, double py, double pz, double borrowed_energy = 0);
  ZBoson(const ZBoson& other, bool copy_decay_products = false); // Copy constructor declaration
  ZBoson& operator=(const ZBoson& other); // Copy assignment operator declaration
  ZBoson(ZBoson&& other) noexcept;
  ZBoson& operator=(ZBoson&& other) noexcept;
  virtual ~ZBoson() = default;

  void decay() override;
  void print() const override;
  static constexpr double get_Z_mass();

  std::shared_ptr<Particle> clone() const override;
};

class HiggsBoson : public Boson
{
private:
  static constexpr double higgs_mass = 125110;
  std::string decay_type;

public:
  HiggsBoson(double px, double py, double pz);
  HiggsBoson(const HiggsBoson& other, bool copy_decay_products = false); // Copy constructor declaration
  HiggsBoson& operator=(const HiggsBoson& other); // Copy assignment operator declaration
  HiggsBoson(HiggsBoson&& other) noexcept;
  HiggsBoson& operator=(HiggsBoson&& other) noexcept;
  virtual ~HiggsBoson() = default;

  void decay() override;
  void print() const override;

  std::shared_ptr<Particle> clone() const override;
};

class Gluon : public Boson
{
private:
  ColourCharge colour1, colour2;
  void check_colour_consistency();

public:
  Gluon(ColourCharge colour1, ColourCharge colour2, double px, double py, double pz);
  Gluon(const Gluon& other);  // Copy constructor
  Gluon(Gluon&& other) noexcept;  // Move constructor
  Gluon& operator=(const Gluon& other);  // Copy assignment operator
  Gluon& operator=(Gluon&& other) noexcept;  // Move assignment operator
  virtual ~Gluon() override = default;  // Destructor

  void decay() override;
  void print() const override;

  std::shared_ptr<Particle> clone() const override;

};

#endif // BOSON_H
