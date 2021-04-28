#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <strstream>
#include "olcPixelGameEngine.h"
//import linalg;

struct Vec3 {
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

struct Triangle {
  std::array<Vec3,3> points;
  olc::Pixel color;

  void translate(float dx, float dy, float dz) {
    points[0].translate(dx,dy,dz);
    points[1].translate(dx,dy,dz);
    points[2].translate(dx,dy,dz);
  }
  void scale(float sx, float sy, float sz) {
    points[0].scale(sx,sy,sz);
    points[1].scale(sx,sy,sz);
    points[2].scale(sx,sy,sz);
  }
};

struct Mesh {
  std::vector<Triangle> tris;

  // Load from .obj file
  bool loadObj(std::string fp) {
    std::ifstream f(fp);
    if (!f.is_open()) {
      return false;
    }

    tris.clear();
    std::vector<Vec3> verts;
    while (!f.eof()) {
      char line[256];
      f.getline(line, 256);
      std::strstream s;
      s << line;
      char junk;
      
      if (line[0] == 'v') {
        Vec3 v;
        s >> junk >> v.x >> v.y >> v.z;
        verts.push_back(v);
      }
      else if (line[0] == 'f') {
        int f[3];
        s >> junk >> f[0] >> f[1] >> f[2];
        tris.push_back({verts[f[0]-1], verts[f[1]-1], verts[f[2]-1]});
      }
    }
    return true;
  }
};

// 4x4 Matrix
struct Mat4 {
  std::array<std::array<float,4>, 4> data = {0};

  float get(int row, int col) {
    return data[row][col];
  }
  void set(int row, int col, float n) {
    data[row][col] = n;
  }
};

// Multiply a 3d Vector by the Projection Matrix by first pretending it's a 4d vector
Vec3 mulProj(Vec3 &u, Mat4 &m) {
  Vec3 res;

  res.x = u.x * m.get(0,0) + u.y * m.get(1,0) + u.z * m.get(2,0) + m.get(3,0);
  res.y = u.x * m.get(0,1) + u.y * m.get(1,1) + u.z * m.get(2,1) + m.get(3,1);
  res.z = u.x * m.get(0,2) + u.y * m.get(1,2) + u.z * m.get(2,2) + m.get(3,2);
  float w = u.x * m.get(0,3) + u.y * m.get(1,3) + u.z * m.get(2,3) + m.get(3,3);

  if (w != 0.0f) {
    res.x /= w;
    res.y /= w;
    res.z /= w;
  }

  return res;
}

Triangle project(Triangle &tri, Mat4 &m) {
  Triangle res;
  res.points[0] = mulProj(tri.points[0], m);
  res.points[1] = mulProj(tri.points[1], m);
  res.points[2] = mulProj(tri.points[2], m);
  return res;
}

olc::Pixel getColor(float lum) {
  if (lum < 0.0f) {
    return olc::BLACK;
  }
  else {
    return olc::WHITE * lum;
  }
}


Mat4 projectionMatrix(float fnear, float ffar, float fov, float aspectRatio) {
  //putting an f in the front makes an "expected qualified-id" error go away,
  // because near and far are defined
  float fovRad = 1.0f / tanf(fov * 0.5f / 180.0f * 3.1415926f);
  Mat4 res;
  res.set(0,0, aspectRatio * fovRad);
  res.set(1,1, fovRad);
  res.set(2,2, ffar / (ffar-fnear));
  res.set(3,2, (-ffar * fnear) / (ffar-fnear));
  res.set(2,3, 1.0f);
  res.set(3,3, 0.0f);
  return res;
}

Mat4 matRotZ(float theta) {
  Mat4 res;
  res.set(0,0, cosf(theta));
  res.set(0,1, sinf(theta));
  res.set(1,0, -sinf(theta));
  res.set(1,1, cosf(theta));
  res.set(2,2, 1.0f);
  res.set(3,3, 1.0f);
  return res;
}

Mat4 matRotX(float theta) {
  Mat4 res;
  res.set(0,0, 1.0f);
  res.set(1,1, cosf(theta));
  res.set(1,2, sinf(theta));
  res.set(2,1, -sinf(theta));
  res.set(2,2, cosf(theta));
  res.set(3,3, 1.0f);
  return res;
}

class Game : public olc::PixelGameEngine {
  Mesh meshCube;
  Mat4 projMat;
  Vec3 camera = {0};

  float theta;

public:
  Game() {
    sAppName = "Test";
  }

  bool OnUserCreate() override {
    if(!meshCube.loadObj("assets/cobra.obj")) {
      std::cout << "Could not load object file\n";
    }

    projMat = projectionMatrix(0.1f, 1000.0f, 90.0f, static_cast<float>(ScreenHeight()) / static_cast<float>(ScreenWidth()));
    return true;
  }

  bool OnUserUpdate(float dt) override {
    Clear(olc::BLACK);

    Mat4 rotZ, rotX;
    theta += 1.0f * dt;
    rotZ = matRotZ(theta);
    rotX = matRotX(theta * 0.5f);

    std::vector<Triangle> trisToDraw;

    // Calculate triangle projection and add to draw queue
    for(auto tri: meshCube.tris) {
      Triangle triProj, triTrans, triRotZ, triRotZX;
      // Rotate in Z Axis
      triRotZ = project(tri, rotZ);
      // Rotate in X Axis
      triRotZX = project(triRotZ, rotX);

      // Translate away from screen
      triTrans = triRotZX;
      triTrans.translate(0,0,8.0f);

      Vec3 normal, u,v;
      u = triTrans.points[1] - triTrans.points[0];
      v = triTrans.points[2] - triTrans.points[0];
      normal = u.cross(v);
      normal.normalize();

      if (normal.dot(triTrans.points[0]-camera) < 0.0f) {
        Vec3 light_dir = {0.0f, 0.0f, -1.0f};
        light_dir.normalize();

        float dp = normal.dot(light_dir);
        triTrans.color = getColor(dp);

        triProj = project(triTrans, projMat);
        triProj.translate(1.0f, 1.0f, 0.0f);
        triProj.scale(0.5f * static_cast<float>(ScreenWidth()), 0.5f * static_cast<float>(ScreenHeight()), 1.0f);
        triProj.color = triTrans.color;

        trisToDraw.push_back(triProj);
      }
    }

    // Sort triangles from back to front
    sort(trisToDraw.begin(), trisToDraw.end(), [](Triangle &t1, Triangle &t2) {
      float avg1 = (t1.points[0].z + t1.points[1].z + t1.points[2].z) / 3.0f;
      float avg2 = (t2.points[0].z + t2.points[1].z + t2.points[2].z) / 3.0f;
      return avg1 > avg2;
    });

    // Draw sorted triangles
    for (auto tri: trisToDraw) {
      FillTriangle(
        tri.points[0].x, tri.points[0].y,
        tri.points[1].x, tri.points[1].y,
        tri.points[2].x, tri.points[2].y,
        tri.color
      );
      // Draw wireframe
      /*
      DrawTriangle(
        tri.points[0].x, tri.points[0].y,
        tri.points[1].x, tri.points[1].y,
        tri.points[2].x, tri.points[2].y,
        olc::BLACK
      );*/
    }
    return true;
  }
};

int main() {
  Game demo;
  if (demo.Construct(320,200,2,2,false, true)) {
    demo.Start();
  }
  return 0;
}
