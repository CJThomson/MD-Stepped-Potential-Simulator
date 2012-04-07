//Class to allow 3D vector calculations
#pragma once
#include<math.h>
class CVector3
{
    public:
    //Constructors
    CVector3(){x=0;y=0;z=0;};
    CVector3(double a1, double a2, double a3)
    {
        x=a1;y=a2;z=a3;
    };
    //Deconstructor
    ~CVector3(){};
    //Operator Overloads
    CVector3 operator+ (const CVector3 &v) const
    {
        CVector3 v2;
        v2.x=x+v.x;
        v2.y=y+v.y;
        v2.z=z+v.z;
        return v2;
    };
    CVector3 operator- (const CVector3 &v) const
    {
        CVector3 v2;
        v2.x=x-v.x;
        v2.y=y-v.y;
        v2.z=z-v.z;
        return v2;
    };
    CVector3 operator* (const double c) const
    {
        CVector3 v2;
        v2.x=x*c;
        v2.y=y*c;
        v2.z=z*c;
        return v2;
    };
    CVector3 operator/ (const double c) const
    {
        CVector3 v2;
        v2.x=x/c;
        v2.y=y/c;
        v2.z=z/c;
        return v2;
    };
    const CVector3& operator+= (const CVector3 &v)
    {
		x+=v.x;
		y+=v.y;
		z+=v.z;
		return *this;
    };
    const CVector3& operator-= (const CVector3 &v)
    {
		x-=v.x;
		y-=v.y;
		z-=v.z;
		return *this;
    };
    const CVector3& operator*= (const CVector3 &v)
    {
		x*=v.x;
		y*=v.y;
		z*=v.z;
		return *this;
    };
    friend CVector3 operator* (const double c,const CVector3 &v);
    friend CVector3 operator/ (const double c,const CVector3 &v);
    //Functions
	void zero()
	{
		x=0;y=0;z=0;
	}
    double length()
    {
        return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    };
    CVector3 normalise()
    {
        double abs = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        double x1=x/abs;
        double y1=y/abs;
        double z1=z/abs;
        return CVector3(x1,y1,z1);
    };
    double dotProd(const CVector3 &v)
    {
        return x*v.x + y*v.y+z*v.z;
    }
    CVector3 crossProd(const CVector3 &v)
    {
        double x1 = y*v.z-z*v.y;
        double y1 = x*v.z-z*v.y;
        double z1 = x*v.y-y*v.x;
        return CVector3(x1,y1,z1);
    }
    double angle(const CVector3 &v)
    {
        return acos((x*v.x + y*v.y+z*v.z)/(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*sqrt(pow(v.x,2)+pow(v.y,2)+pow(v.z,2))));
    }
    double x,y,z;

    inline double& operator[](const size_t idx)
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

    inline const double& operator[](const size_t idx) const
    {
      switch (idx) { case 0: return x; case 1: return y; case 2: return z; default: exit(1); }
    }
};
CVector3 operator* (const double c,const CVector3 &v)
{
    return v*c;
};
CVector3 operator/ (const double c,const CVector3 &v)
{
    return v/c;
};

