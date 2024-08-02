#ifndef QUARK_H
#define QUARK_H

#include "particle.h"
#include "fourmom.h"
#include <string>
#include <cmath>

enum class ColourCharge
{
  Red, Green, Blue, // Colours
  AntiRed, AntiGreen, AntiBlue, // Anticolours
  Neutral // For default
};

std::string colour_charge_to_string(ColourCharge colour);

class Quark : public Particle
{
protected:
  ColourCharge colour;
  double baryon_number;

public:
  Quark(double mass, double charge, double px, double py, double pz, const std::string& type, ColourCharge colour, bool is_anti,
    double baryon_number);
  Quark(const Quark& other); // Copy constructor
  Quark(Quark&& other) noexcept; // Move constructor
  Quark& operator=(const Quark& other); // Copy assignment operator
  Quark& operator=(Quark&& other) noexcept; // Move assignment operator
  virtual ~Quark() = default;
  virtual void print() const override;
  virtual void decay() override; // Making Quark an abstract class since all Quarks must implement decay()
  void check_colour_consistency();

  double get_baryon_number() const override;
};


class UpQuark : public Quark
{
private:
  static constexpr double up_mass = 2.2;

public:
  UpQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  UpQuark(const UpQuark& other); // Copy constructor
  UpQuark(UpQuark&& other) noexcept; // Move constructor
  UpQuark& operator=(const UpQuark& other); // Copy assignment operator
  UpQuark& operator=(UpQuark&& other) noexcept; // Move assignment operator
  virtual ~UpQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

class DownQuark : public Quark
{
private:
  static constexpr double down_mass = 4.7;

public:
  DownQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  DownQuark(const DownQuark& other); // Copy constructor
  DownQuark(DownQuark&& other) noexcept; // Move constructor
  DownQuark& operator=(const DownQuark& other); // Copy assignment operator
  DownQuark& operator=(DownQuark&& other) noexcept; // Move assignment operator
  virtual ~DownQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

class CharmQuark : public Quark
{
private:
  static constexpr double charm_mass = 1280;

public:
  CharmQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  CharmQuark(const CharmQuark& other); // Copy constructor
  CharmQuark(CharmQuark&& other) noexcept; // Move constructor
  CharmQuark& operator=(const CharmQuark& other); // Copy assignment operator
  CharmQuark& operator=(CharmQuark&& other) noexcept; // Move assignment operator
  virtual ~CharmQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

class StrangeQuark : public Quark
{
private:
  static constexpr double strange_mass = 95;

public:
  StrangeQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  StrangeQuark(const StrangeQuark& other); // Copy constructor
  StrangeQuark(StrangeQuark&& other) noexcept; // Move constructor
  StrangeQuark& operator=(const StrangeQuark& other); // Copy assignment operator
  StrangeQuark& operator=(StrangeQuark&& other) noexcept; // Move assignment operator
  virtual ~StrangeQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

class TopQuark : public Quark
{
private:
  static constexpr double top_mass = 173100;

public:
  TopQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  TopQuark(const TopQuark& other); // Copy constructor
  TopQuark(TopQuark&& other) noexcept; // Move constructor
  TopQuark& operator=(const TopQuark& other); // Copy assignment operator
  TopQuark& operator=(TopQuark&& other) noexcept; // Move assignment operator
  virtual ~TopQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

class BottomQuark : public Quark
{
private:
  static constexpr double bottom_mass = 4180;

public:
  BottomQuark(double px, double py, double pz, ColourCharge colour, bool is_anti = false);
  BottomQuark(const BottomQuark& other); // Copy constructor
  BottomQuark(BottomQuark&& other) noexcept; // Move constructor
  BottomQuark& operator=(const BottomQuark& other); // Copy assignment operator
  BottomQuark& operator=(BottomQuark&& other) noexcept; // Move assignment operator
  virtual ~BottomQuark() = default;
  void decay() override;

  std::shared_ptr<Particle> clone() const override;
};

#endif // QUARK_H
