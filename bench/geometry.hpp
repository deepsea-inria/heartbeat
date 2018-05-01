
#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>

#ifndef _PBBS_SPTL_GEOMETRY_
#define _PBBS_SPTL_GEOMETRY_

namespace sptl {

#define EPS 1e-20

// *************************************************************
//    POINTS AND VECTORS (3d),  2d is below
// *************************************************************

template <class _floatT> class _point3d;

template <class _floatT> class _vect3d {
public:
  typedef _floatT floatT;
  typedef _vect3d vectT;
  typedef _point3d<floatT> pointT;
  floatT x; floatT y; floatT z;

  _vect3d(floatT xx,floatT yy, floatT zz) : x(xx),y(yy),z(zz) {}

  _vect3d(): x(0), y(0), z(0) {}

  _vect3d(pointT p);

  _vect3d(floatT* p) : x(p[0]), y(p[1]), z(p[2]) {};

  vectT operator+(vectT op2) {
    return vectT(x + op2.x, y + op2.y, z + op2.z);
  }

  vectT operator-(vectT op2) {
    return vectT(x - op2.x, y - op2.y, z - op2.z);
  }

  pointT operator+(pointT op2);

  vectT operator*(floatT s) {
    return vectT(x * s, y * s, z * s);
  }

  vectT operator/(floatT s) {
    return vectT(x / s, y / s, z / s);
  }

  floatT& operator[] (int i) {
    return (i == 0) ? x : (i == 1) ? y : z;
  }
  floatT dot(vectT v) {return x * v.x + y * v.y + z * v.z;}

  vectT cross(vectT v) {
    return vectT(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
  }

  floatT max_dim() {
    return std::max(x, std::max(y, z));
  }

  void print() {
    std::cout << std::setprecision(10) << ":(" << x << "," << y << "," << z << "):";
  }

  floatT length(void) {
    return sqrt(x * x + y * y + z * z);
  }

  static const int dim = 3;
};

template <class _floatT> class _point3d {
public:
  typedef _floatT floatT;
  typedef _vect3d<floatT> vectT;
  typedef _point3d pointT;
  floatT x; floatT y; floatT z;
  int dimension() {return 3;}
  _point3d(floatT xx,floatT yy, floatT zz) : x(xx),y(yy),z(zz) {}
  _point3d() {x=0;y=0;z=0;}
  _point3d(vectT v) : x(v.x),y(v.y),z(v.z) {};
  _point3d(floatT* p) : x(p[0]), y(p[1]), z(p[2]) {};
  void print() {std::cout << ":(" << x << "," << y << "," << z << "):";}
  vectT operator-(pointT op2) {
    return vectT(x - op2.x, y - op2.y, z - op2.z);}
  pointT operator+(vectT op2) {
    return pointT(x + op2.x, y + op2.y, z + op2.z);}

  pointT min_coord(pointT b) {
    return pointT(std::min(x, b.x), std::min(y, b.y), std::min(z, b.z));
  }

  pointT max_coord(pointT b) {
    return pointT(std::max(x, b.x), std::max(y, b.y), std::max(z, b.z)); 
  }

  floatT& operator[] (int i) {return (i==0) ? x : (i==1) ? y : z;}
  int quadrant(pointT center) {
    int index = 0;
    if (x > center.x) index += 1;
    if (y > center.y) index += 2;
    if (z > center.z) index += 4;
    return index;
  }
  // returns a pointT offset by offset in one of 8 directions
  // depending on dir (an integer from [0..7])
  pointT offset_point(int dir, floatT offset) {
    floatT xx = x + ((dir & 1) ? offset : -offset);
    floatT yy = y + ((dir & 2) ? offset : -offset);
    floatT zz = z + ((dir & 4) ? offset : -offset);
    return pointT(xx, yy, zz);
  }
  // checks if pt is outside of a box centered at this point with
  // radius hsize
  bool out_of_box(pointT pt, floatT hsize) {
    return ((x - hsize > pt.x) || (x + hsize < pt.x) ||
            (y - hsize > pt.y) || (y + hsize < pt.y) ||
            (z - hsize > pt.z) || (z + hsize < pt.z));
  }
  static const int dim = 3;
};

template <class floatT>
inline _point3d<floatT> _vect3d<floatT>::operator+(_point3d<floatT> op2) {
  return _point3d<floatT>(x + op2.x, y + op2.y, z + op2.z);}

template <class floatT>
inline _vect3d<floatT>::_vect3d(_point3d<floatT> p) { x = p.x; y = p.y; z = p.z;}

// *************************************************************
//    POINTS AND VECTORS (2d)
// *************************************************************

template <class floatT> class _point2d;

template <class _floatT> class _vect2d {
public:
  typedef _floatT floatT;
  typedef _point2d<floatT> pointT;
  typedef _vect2d vectT;
  floatT x;
  floatT y;

  _vect2d(floatT xx,floatT yy) : x(xx),y(yy) {}

  _vect2d(): x(0), y(0) { }

  _vect2d(pointT p);

  _vect2d(floatT* p) : x(p[0]), y(p[1]) {};

  vectT operator+(vectT op2) {
   return vectT(x + op2.x, y + op2.y);
  }

  vectT operator-(vectT op2) {
    return vectT(x - op2.x, y - op2.y);
  }

  pointT operator+(pointT op2);

  vectT operator*(floatT s) {
    return vectT(x * s, y * s);
  }

  vectT operator/(floatT s) {
    return vectT(x / s, y / s);
  }

  floatT operator[] (int i) {
    return (i == 0) ? x : y;
  }

  floatT dot(vectT v) {
    return x * v.x + y * v.y;
  }

  floatT cross(vectT v) {
    return x*v.y - y*v.x; 
  }

  floatT max_dim() {
    return std::max(x, y);
  }

  void print() {
    std::cout << ":(" << x << "," << y << "):";
  }

  floatT length(void) {
    return sqrt(x * x + y * y);
  }

  static const int dim = 3;
};

template <class floatT>
static std::ostream& operator<<(std::ostream& os, const _vect3d<floatT> v) {
  return os << v.x << " " << v.y << " " << v.z; }

template <class floatT>
static std::ostream& operator<<(std::ostream& os, const _point3d<floatT> v) {
  return os << v.x << " " << v.y << " " << v.z;
}

template <class _floatT> class _point2d {
public:
  typedef _floatT floatT;
  typedef _vect2d<floatT> vectT;
  typedef _point2d pointT;
  floatT x; floatT y;
  int dimension() {return 2;}
  _point2d(floatT xx,floatT yy) : x(xx),y(yy) {}
  _point2d() {x=0;y=0;}
  _point2d(vectT v) : x(v.x),y(v.y) {};
  _point2d(floatT* p) : x(p[0]), y(p[1]) {};
  void print() {std::cout << ":(" << x << "," << y << "):";}

  vectT operator-(pointT op2) {
    return vectT(x - op2.x, y - op2.y);
  }

  pointT operator+(vectT op2) {
    return pointT(x + op2.x, y + op2.y);
  }

  floatT operator[] (int i) {
    return (i==0) ? x : y;
  }

  pointT min_coord(pointT b) {
    return pointT(std::min(x, b.x), std::min(y, b.y));
  }

  pointT max_coord(pointT b) {
    return pointT(std::max(x, b.x), std::max(y, b.y));
  }

  int quadrant(pointT center) {
    int index = 0;
    if (x > center.x) index += 1;
    if (y > center.y) index += 2;
    return index;
  }
  // returns a pointT offset by offset in one of 4 directions
  // depending on dir (an integer from [0..3])
  pointT offset_point(int dir, floatT offset) {
    floatT xx = x + ((dir & 1) ? offset : -offset);
    floatT yy = y + ((dir & 2) ? offset : -offset);
    return pointT(xx, yy);
  }
  bool out_of_box(pointT pt, floatT hsize) {
    return ((x - hsize > pt.x) || (x + hsize < pt.x) ||
            (y - hsize > pt.y) || (y + hsize < pt.y));
  }
  static const int dim = 2;
};

template <class floatT>
inline _point2d<floatT> _vect2d<floatT>::operator+(_point2d<floatT> op2) {
  return _point2d<floatT>(x + op2.x, y + op2.y);}

template <class floatT>
inline _vect2d<floatT>::_vect2d(_point2d<floatT> p) { x = p.x; y = p.y;}

template <class floatT>
static std::ostream& operator<<(std::ostream& os, const _vect2d<floatT> v) {
  return os << v.x << " " << v.y;}

template <class floatT>
static std::ostream& operator<<(std::ostream& os, const _point2d<floatT> v) {
  return os << v.x << " " << v.y; }

// *************************************************************
//    SPECIALIZING TO DOUBLE PRECISION
// *************************************************************

#define GCTYPE double
typedef _vect2d<GCTYPE> vect2d;
typedef _point2d<GCTYPE> point2d;
typedef _vect3d<GCTYPE> vect3d;
typedef _point3d<GCTYPE> point3d;

// *************************************************************
//    GEOMETRY
// *************************************************************

// Returns twice the area of the oriented triangle (a, b, c)
inline double triangle_area(point2d a, point2d b, point2d c) {
  return (b - a).cross(c - a);
}

inline double triAreaNormalized(point2d a, point2d b, point2d c) {
  return triangle_area(a, b, c) / ((b - a).length() * (c - a).length());
}

// Returns TRUE if the points a, b, c are in a counterclockise order
inline int counter_clockwise(point2d a, point2d b, point2d c) {
  return (b - a).cross(c - a) > 0.0;
}

template <class floatT>
inline _vect3d<floatT> onParabola(_vect2d<floatT> v) {
  return vect3d(v.x, v.y, v.x*v.x + v.y*v.y);}

// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c.
// Projects a, b, c onto a parabola centered with d at the origin
//   and does a plane side test (tet volume > 0 test)
template <class floatT>
inline bool inCircle(_point2d<floatT> a, _point2d<floatT> b,
                     _point2d<floatT> c, _point2d<floatT> d) {
  _vect3d<floatT> ad = onParabola(a-d);
  _vect3d<floatT> bd = onParabola(b-d);
  _vect3d<floatT> cd = onParabola(c-d);
  return (ad.cross(bd)).dot(cd) > 0.0;
}

// returns a number between -1 and 1, such that -1 is out at infinity,
// positive numbers are on the inside, and 0 is at the boundary
template <class floatT>
inline double inCircleNormalized(_point2d<floatT> a, _point2d<floatT> b,
                                 _point2d<floatT> c, _point2d<floatT> d) {
  _vect3d<floatT> ad = onParabola(a - d);
  _vect3d<floatT> bd = onParabola(b - d);
  _vect3d<floatT> cd = onParabola(c - d);
  return (ad.cross(bd)).dot(cd) / (ad.length() * bd.length() * cd.length());
}

// *************************************************************
//    TRIANGLES
// *************************************************************

struct triangle {
  int vertices[3];

  triangle() {}

  triangle(int p1, int p2, int p3) {
    vertices[0] = p1;
    vertices[1] = p2;
    vertices[2] = p3;
  }
};

static std::ostream& operator<<(std::ostream& os, const triangle t) {
  return os << t.vertices[0] << " " << t.vertices[1] << " " << t.vertices[2];
}

template <class point>
struct triangles {
  long num_points;
  long num_triangles;
  point* p;
  triangle* t;
  triangles() {}
  void del() {
    free(p);
    free(t);
  }
  
  triangles(long np, long nt, point* _p, triangle* _t)
    : num_points(np), num_triangles(nt), p(_p), t(_t) {}
};

template <class pointT>
struct ray {
  typedef typename pointT::vectT vectT;
  pointT o;
  vectT d;
  ray(pointT _o, vectT _d) : o(_o), d(_d) {}
  ray() {}
};

inline double angle(point2d a, point2d b, point2d c) {
  vect2d ba = (b - a);
  vect2d ca = (c - a);
  double lba = ba.length();
  double lca = ca.length();
  double pi = 3.14159;
  return 180. / pi * acos(ba.dot(ca) / (lba * lca));
}

inline double minAngleCheck(point2d a, point2d b, point2d c, double angle) {
  vect2d ba = (b - a);
  vect2d ca = (c - a);
  vect2d cb = (c - b);
  double lba = ba.length();
  double lca = ca.length();
  double lcb = cb.length();
  double pi = 3.14159;
  double co = cos(angle * pi / 180.);
  return (ba.dot(ca) / (lba * lca) > co || ca.dot(cb) / (lca * lcb) > co ||
          -ba.dot(cb) / (lba * lcb) > co);
}

inline point2d triangleCircumcenter(point2d a, point2d b, point2d c) {
  vect2d v1 = b - a;
  vect2d v2 = c - a;
  vect2d v11 = v1 * v2.dot(v2);
  vect2d v22 = v2 * v1.dot(v1);
  return  a + vect2d(v22.y - v11.y, v11.x - v22.x)/(2.0 * v1.cross(v2));
}
  
} // end namespace

#endif
