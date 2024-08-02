#ifndef LEPTON_H
#define LEPTON_H

#include "fourmom.h"
#include "particle.h"
#include <memory>
#include <string>
#include <vector>

class Lepton : public Particle
{
protected:
  int electron_lepton_number;
  int muon_lepton_number;
  int tau_lepton_number;

public:
  Lepton(double mass, double charge, double px, double py, double pz, const std::string& type, bool is_anti,
         int electron_number = 0, int muon_number = 0, int tau_number = 0);
  Lepton(const Lepton& other, bool copy_decay_products = true); // Copy constructor
  Lepton(Lepton&& other) noexcept; // Move constructor
  Lepton& operator=(const Lepton& other); // Copy assignment operator
  Lepton& operator=(Lepton&& other) noexcept; // Move assignment operator
  virtual ~Lepton() = default; // Destructor

  virtual void print() const override;
  virtual void decay() override;

  // Override getters to return the actual numbers
  virtual int get_electron_lepton_number() const override;
  virtual int get_muon_lepton_number() const override;
  virtual int get_tau_lepton_number() const override;
};

class Electron : public Lepton
{
private:
  std::vector<double> calorimeter_deposits;
  static constexpr double electron_mass = 0.511;

public:
  Electron(double px = 0, double py = 0, double pz = 0, std::vector<double> deposits = {}, bool is_anti = false);
  Electron(const Electron& other); // Copy constructor
  Electron(Electron&& other) noexcept; // Move constructor
  Electron& operator=(const Electron& other); // Copy assignment operator
  Electron& operator=(Electron&& other) noexcept; // Move assignment operator
  virtual ~Electron() = default;

  void print() const override;
  void decay() override;
  int get_electron_lepton_number() const override;
  void adjust_calorimeter_deposits();

  std::shared_ptr<Particle> clone() const override;
};

class Muon : public Lepton
{
private:
  bool is_isolated;
  static constexpr double muon_mass = 105.66;

public:
  Muon(double px=0, double py=0, double pz=0, bool isolated = false, bool isAnti = false);
  Muon(const Muon& other); // Copy constructor
  Muon(Muon&& other) noexcept; // Move constructor
  Muon& operator=(const Muon& other); // Copy assignment operator
  Muon& operator=(Muon&& other) noexcept; // Move assignment operator
  virtual ~Muon() = default;

  void print() const override;
  void decay() override;
  int get_muon_lepton_number() const override;

  std::shared_ptr<Particle> clone() const override;
};

class Tau : public Lepton
{
private:
  std::string decay_type;
  static constexpr double tau_mass = 1776.86;

public:
  Tau(double px=0, double py=0, double pz=0, bool isAnti = false);
  Tau(const Tau& other, bool copy_decay_products = false); // Copy constructor
  Tau(Tau&& other) noexcept; // Move constructor
  Tau& operator=(const Tau& other); // Copy assignment operator
  Tau& operator=(Tau&& other) noexcept; // Move assignment operator
  virtual ~Tau() = default;


  void print() const override;
  void decay() override;
  int get_tau_lepton_number() const override;

  std::shared_ptr<Particle> clone() const override;
};

class ElectronNeutrino : public Lepton
{
public:
  // Assuming neutrinos are not anti-particles by default and do not carry lepton numbers other than their own
  ElectronNeutrino(double px=0, double py=0, double pz=0, bool interacted = false, bool isAnti = false);
  ElectronNeutrino(const ElectronNeutrino& other); // Copy constructor
  ElectronNeutrino(ElectronNeutrino&& other) noexcept; // Move constructor
  ElectronNeutrino& operator=(const ElectronNeutrino& other); // Copy assignment operator
  ElectronNeutrino& operator=(ElectronNeutrino&& other) noexcept; // Move assignment operator
  virtual ~ElectronNeutrino() = default;

  void print() const override;
  void decay() override;
  int get_electron_lepton_number() const override;

  std::shared_ptr<Particle> clone() const override;

private:
  bool has_interacted;
  static constexpr double electron_neutrino_mass = 0; // Assuming mass is very small for electron neutrinos
};

class MuonNeutrino : public Lepton
{
public:
  MuonNeutrino(double px=0, double py=0, double pz=0, bool interacted = false, bool isAnti = false);
  MuonNeutrino(const MuonNeutrino& other); // Copy constructor
  MuonNeutrino(MuonNeutrino&& other) noexcept; // Move constructor
  MuonNeutrino& operator=(const MuonNeutrino& other); // Copy assignment operator
  MuonNeutrino& operator=(MuonNeutrino&& other) noexcept; // Move assignment operator
  virtual ~MuonNeutrino() = default;

  void print() const override;
  void decay() override;
  int get_muon_lepton_number() const override;
  std::shared_ptr<Particle> clone() const override;
  
private:
  bool has_interacted;
  static constexpr double muon_neutrino_mass = 0; // Assuming mass is very small for muon neutrinos
};

class TauNeutrino : public Lepton
{
public:
  TauNeutrino(double px=0, double py=0, double pz=0, bool interacted = false, bool isAnti = false);
  TauNeutrino(const TauNeutrino& other); // Copy constructor
  TauNeutrino(TauNeutrino&& other) noexcept; // Move constructor
  TauNeutrino& operator=(const TauNeutrino& other); // Copy assignment operator
  TauNeutrino& operator=(TauNeutrino&& other) noexcept; // Move assignment operator
  virtual ~TauNeutrino() = default;


  void print() const override;
  void decay() override;
  int get_tau_lepton_number() const override;

  std::shared_ptr<Particle> clone() const override;
  
private:
  bool has_interacted;
  static constexpr double tau_neutrino_mass = 0; // Assumning mass is very small for tau neutrinos
};

#endif // LEPTON_H
