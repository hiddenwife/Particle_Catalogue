#ifndef FOURMOM_H
#define FOURMOM_H

#include <array>

class FourMomentum
{
private:
    double E, px, py, pz;

public:
  // Constructors
  FourMomentum(double E, double px, double py, double pz);
  FourMomentum(const FourMomentum& other);
  FourMomentum(FourMomentum&& other) noexcept;
  FourMomentum& operator=(FourMomentum&& other) noexcept;
  ~FourMomentum();

  void set_e(double E);
  void set_px(double px);
  void set_py(double py);
  void set_pz(double pz);

  double get_e() const;
  double get_px() const;
  double get_py() const;
  double get_pz() const;

  double invariant_mass() const;

  friend FourMomentum operator+(const FourMomentum& lhs, const FourMomentum& rhs);
  friend FourMomentum operator-(const FourMomentum& lhs, const FourMomentum& rhs);
  friend double dot_product(const FourMomentum& lhs, const FourMomentum& rhs);
};

#endif // FOURMOM_H
