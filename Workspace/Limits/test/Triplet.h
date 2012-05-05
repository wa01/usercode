#ifndef Triplet_h_
#define Triplet_h_
struct Triplet {
  Triplet () : x_(0), y_(0), z_(0.) {}
  Triplet (int x, int y, double z) : x_(x), y_(y), z_(z) {}
  inline int x () const {return x_;}
  inline int y () const {return y_;}
  inline double z () const {return z_;}
  int x_,y_;
  double z_;
};
#endif
