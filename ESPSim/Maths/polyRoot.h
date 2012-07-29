#pragma once

/* = = = = = = = = = = = = = = = = = = = = = = = 
   = = = = = = = Polynomial  Class = = = = = = = 
   = = = = = = = = = = = = = = = = = = = = = = = */ 

//This class performs functions on polynomials, including root findings, solving, and differentiating
#include<math.h>
#include<iostream>
class Polynomial{
 private:
  void synthDiv(std::vector<double>& poly, double root)
  {
    std::vector<double> intermediate(poly.size());
    std::vector<double> answer(poly.size() - 1);
    for(int i = 0; i < poly.size() - 1; ++i)
      {
        if(i == 0)
	  intermediate[i] = 0;
        else
	  intermediate[i] = answer[i - 1] * root;

        answer[i] = poly[i] + intermediate[i];
      }
    poly[0] = 0;
    for(int i = 1; i < poly.size(); ++i)
      poly[i] = answer[i - 1];
  }
 public:
  std::vector<double> polynomial; //create a new vector that is going to act as a storage structure of the polynomial
 Polynomial(std::vector<double> poly) : polynomial(poly) {} //Constructor
  Polynomial(){}
  ~Polynomial() { polynomial.clear(); } //Deconstructor

  //Operator Overloads
  inline double& operator[](const size_t idx)
  {
    return polynomial[idx];
  }

  inline const double& operator[](const size_t idx) const
  {
    return polynomial[idx];
  }
  void rootFind(std::vector<double> poly, std::vector<double>& roots, 
		double tolerance, int iterations)
  {
    double root;
    double root0 = 0;
    bool noRoot = false;
    std::vector<double> Pdot;
    std::vector<double> P2dot;
    polyDiff(poly, Pdot);
    polyDiff(Pdot, P2dot);
    for(int n = 0; n < poly.size() - 1; ++n)
      {
	//writePoly(poly);
        int count = 0;
        root = poly[poly.size() - 1];
	root0 = 0;
        do{
	  count++;
	  root0 = root;
	  root = root0 - ( 2 * polySolve(root0, poly) * polySolve(root0,Pdot) )
	    /( 2 * pow(polySolve(root0,Pdot), 2) - polySolve(root0,poly) * 
	       polySolve(root0, P2dot));
	  if(count == iterations)
            {
	      noRoot=true;
	      break;
            }
        }while(fabs(root - root0) > tolerance);

        if(noRoot==true)
	  break;
	//std::cout << root<< std::endl;
	roots.push_back(root);
        synthDiv(poly, root);
        polyDiff(poly, Pdot);
        polyDiff(Pdot, P2dot);
      }

  }

  void polyDiff(std::vector<double>& poly, std::vector<double>& diffPoly)
  {
    //this function returns the derivative of a polynomial
    diffPoly.resize(poly.size() - 1);
    for(int i = 0; i < poly.size()-1; ++i)
      diffPoly[i] = poly[i] * (poly.size() - (1 + i));
  }
  double polySolve(double x, std::vector<double>& poly)
  {
    //this solves a polynomial at a specific x
    double sum=0;
    for(int i = 0; i < poly.size(); ++i)
      sum += poly[i] * pow(x, poly.size() - (1 + i));
    return sum;
  }

  void writePoly(std::vector<double>& poly) 
  {
    for(int i = 0; i < poly.size(); ++i)
      std::cout <<" + "<< poly[i]<< "x^"<< poly.size()- (i + 1);
    std::cout << std::endl;
  }
};
