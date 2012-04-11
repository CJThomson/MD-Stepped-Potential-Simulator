struct Results
{
  Results() {}
Results(double T, double P,double U) :
  temperature(T), pressure(P), potential(U) {}
  const Results& operator += (const Results &r)
  {
    temperature += r.temperature;
    pressure += r.pressure;
    potential += r.potential;
    return *this;
  }
  const Results& operator *= (const double &a)
  {
    temperature *= a;
    pressure *= a;
    potential *= a;
    return *this;
  }
  double temperature;
  double pressure;
  double potential;
};
