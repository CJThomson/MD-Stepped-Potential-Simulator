#pragma once
//include boost libraries for random numbers
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <time.h>
class CRandom
{

 public:
  CRandom()
    {}
  boost::mt19937 randomNumberGenerator; //random number generator
  double var_normal()
  {
    //initialise normally distributed random numbers
    boost::normal_distribution <double> nd(0.0, 1.0);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > var_nor(randomNumberGenerator, nd);
    return var_nor();
  }
  double var_01()
  {
    //initialise uniformly distributed random numbers between 0 and 1
    boost::uniform_real<double> ud_real(0.0, 1.0);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > var_01(randomNumberGenerator, ud_real);

    return var_01();
  }
  double var_uniformInt(int lower, int upper)
  {
    boost::uniform_int<int> ud_int(lower, upper);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<int> > uni_int(randomNumberGenerator, ud_int);
    return uni_int();
  }
  void seed()
  {
    randomNumberGenerator.seed(time(NULL));
  }
};
