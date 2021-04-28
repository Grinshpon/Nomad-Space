module;
#include "olcPixelGameEngine.h"
export module linalg;

export struct Vec3 {
  float x,y,z;

  void translate(float dx, float dy, float dz) {
    x += dx;
    y += dy;
    z += dz;
  }

  void normalize() {
    float l = sqrtf(x*x + y*y + z*z);
    x /= l;
    y /= l;
    z /= l;
  }
  void scale(float sx, float sy, float sz) {
    x *= sx;
    y *= sy;
    z *= sz;
  }

  Vec3 operator+ (Vec3 &rhs) {
    return {
      x+rhs.x,
      y+rhs.y,
      z+rhs.z
    };
  }
  Vec3 operator- (Vec3 &rhs) {
    return {
      x-rhs.x,
      y-rhs.y,
      z-rhs.z
    };
  }
  float dot(Vec3 v) {
    return x*v.x + y*v.y + z*v.z;
  }
  Vec3 cross(Vec3 &v) {
    return {
      y * v.z - z * v.y,
      z * v.x - x * v.z,
      x * v.y - y * v.x,
    };
  }
};
