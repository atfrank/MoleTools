//Sean M. Law

#include <cmath>

class Vector {
  private:
    double x;
    double y;
    double z;

  public:
    Vector();
    Vector(double xin, double yin, double zin); //Constructor
    Vector(const Vector& vec); //Overload Constructor 

    Vector& operator= (const Vector& vec);
    Vector& operator= (const double val);
    //Addition
    Vector& operator+ (const Vector& vec);
    Vector& operator+= (const Vector& vec);
    Vector& operator+ (const double val);
    Vector& operator+= (const double val);
    //Subtraction
    Vector& operator- (const Vector& vec);
    Vector& operator-= (const Vector& vec);
    Vector& operator- (const double val);
    Vector& operator-= (const double val);
    //Multiplication
    Vector& operator* (const Vector& vec);
    Vector& operator*= (const Vector& vec);
    Vector& operator* (const double val);
    Vector& operator*= (const double val);
    //Division
    Vector& operator/ (const Vector& vec);
    Vector& operator/= (const Vector& vec);
    Vector& operator/ (const double val);
    Vector& operator/= (const double val);

    Vector operator- () const;
    Vector operator- (const Vector& vec) const;
    Vector operator+ (const Vector& vec) const;
    double operator* (const Vector& vec) const; //Dot Product
    Vector cross (const Vector& vec) const; //Cross Product
    double norm () const;
};


