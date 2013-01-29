//Sean M. Law

class Vector {
  private:
    double x;
    double y;
    double z;

  public:
    Vector();
    Vector(double xin, double yin, double zin); //Constructor
    Vector(const Vector& vec); //Overload constructor 

    Vector& operator= (const Vector& vec);
    Vector& operator= (const double val);
    Vector& operator+ (const Vector& vec);
    Vector& operator+= (const Vector& vec);
    Vector& operator+ (const double val);
    Vector& operator+= (const double val);
    Vector& operator- (const Vector& vec);
    Vector& operator-= (const Vector& vec);
    Vector& operator- (const double val);
    Vector& operator-= (const double val);
    Vector& operator* (const Vector& vec);
    Vector& operator*= (const Vector& vec);
    Vector& operator* (const double val);
    Vector& operator*= (const double val);
    Vector& operator/ (const Vector& vec);
    Vector& operator/= (const Vector& vec);
    Vector& operator/ (const double val);
    Vector& operator/= (const double val);

    Vector operator- () const;
    Vector operator- (const Vector& vec) const;
    Vector operator+ (const Vector& vec) const;
    Vector operator* (const Vector& vec) const;
};


