//Sean M. Law

class Vector {
  private:
    double x;
    double y;
    double z;

  public:
    Vector(double xin=0.0, double yin=0.0, double zin=0.0) :
      x(xin), y(yin), z(zin) {} //constructor initialization

    Vector(const Vector& vec){
      x=vec.x;
      y=vec.y;
      z=vec.z;
    }

    Vector& operator=(const Vector& vec){
      x=vec.x;
      y=vec.y;
      z=vec.z;
      return(*this);
    }

    Vector& operator=(const double val){
      x=val;
      y=val;
      z=val;
      return(*this);
    }

};
