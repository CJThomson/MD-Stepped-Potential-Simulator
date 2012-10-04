//Class to allow 3D vector calculations
#pragma once
#include<math.h>
#include<cstdlib>
//declarations that classes and operator overloads are templates...this is so the compiler knows what to expect
template <class T> class Vector3;
template <class T> Vector3<T> operator* (const T c,const Vector3<T> &v);
template <class T> Vector3<T> operator/ (const T c,const Vector3<T> &v);

template<class T> //make the vector3 class a template class to allow doubles/ints/ etc to be used as parameters.
class Vector3
{
 public:
  //Constructors
 Vector3(T a1 = 0, T a2 = 0, T a3 = 0) : x(a1), y(a2), z(a3) {}
  //Deconstructor
  ~Vector3(){};
  //Operator Overloads
  Vector3 operator+ (const Vector3 &v) const
  {
    Vector3 v2;
    v2.x = x + v.x;
    v2.y = y + v.y;
    v2.z = z + v.z;
    return v2;
  };
  Vector3 operator- (const Vector3 &v) const
  {
    Vector3 v2;
    v2.x = x - v.x;
    v2.y = y - v.y;
    v2.z = z - v.z;
    return v2;
  };
  Vector3 operator* (const T c) const
    {
      Vector3 v2;
      v2.x = x * c;
      v2.y = y * c;
      v2.z = z * c;
      return v2;
    };
  Vector3 operator/ (const T c) const
  {
    Vector3 v2;
    v2.x = x / c;
    v2.y = y / c;
    v2.z = z / c;
    return v2;
  };
  const Vector3& operator= (const Vector3 &v)
    {
      x = v.x; y = v.y; z = v.z;
      return *this;
    }
  const Vector3& operator+= (const Vector3 &v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  };
  const Vector3& operator-= (const Vector3 &v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  };
  const Vector3& operator*= (const T c) 
  {
    x *= c;
    y *= c;
    z *= c;
    return *this;
  };
  friend Vector3<T> (::operator* <>) (const T c,const Vector3<T> &v);
  friend Vector3<T> (::operator/ <>) (const T c,const Vector3<T> &v);
  //Functions
  void zero()
  {
    x = 0; y = 0; z = 0;
  }
  double length() const
  {
    return sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2) );
  }
  Vector3 normalise() const
  {
    double abs = sqrt(x*x+y*y+z*z);
    double x1=x/abs;
    double y1=y/abs;
    double z1=z/abs;
    return Vector3(x1,y1,z1);
  }
  T dotProd(const Vector3 &v) const
  {
    return x * v.x + y * v.y + z * v.z;
  }

  T lengthSqr() const
  {
    return x * x + y * y + z * z;
  }
  Vector3 crossProd(const Vector3 &v) const
  {
    double x1 = y * v.z - z * v.y;
    double y1 = x * v.z - z * v.y;
    double z1 = x * v.y - y * v.x;
    return Vector3(x1, y1, z1);
  }
  double angle(const Vector3 &v) const
  {
    return acos((x*v.x + y*v.y+z*v.z)/(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*sqrt(pow(v.x,2)+pow(v.y,2)+pow(v.z,2))));
  }

  inline T& operator[](const size_t idx)
  {
    switch (idx)
      {
      case 0:
	return x;
      case 1:
	return y;
      case 2:
	return z;
      default:
	exit(1);
      }
  }

  inline const T& operator[](const size_t idx) const
  {
    switch (idx) { case 0: return x; case 1: return y; case 2: return z; default: exit(1); }
  }
 protected:
  T x, y, z;
};
template <class T>
Vector3<T> operator* (const T c,const Vector3<T> &v)
{
  return v*c;
};
template <class T>
Vector3<T> operator/ (const T c,const Vector3<T> &v)
{
  return v/c;
};

template <class T>
class PBCVector:public Vector3<T>
{
 private:
  T sysLen;
  bool centreOrigin;
  void applyPBC()
  {
    if(centreOrigin)
      {
	this->x -= lrint(this->x / sysLen) * sysLen;
	this->y -= lrint(this->y / sysLen) * sysLen;
	this->z -= lrint(this->z / sysLen) * sysLen;		
      }
    else
      {
	this->x -= floor(this->x / sysLen) * sysLen;
	this->y -= floor(this->y / sysLen) * sysLen;
	this->z -= floor(this->z / sysLen) * sysLen;
      }
  }
 public:
 PBCVector(T length, bool cenOri, Vector3<T> vec) : Vector3<T>(vec), sysLen(length), centreOrigin(cenOri) {applyPBC();}
  const PBCVector operator= (const Vector3<T> &v)
    {
      this->x = v[0]; this->y = v[1]; this->z = v[2];
      applyPBC();
    }
  
};
