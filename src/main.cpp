#include <array>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <strstream>
#include "olcPixelGameEngine.h"
//import linalg;

struct Vec2 {
  float u,v;


  void normalize() {
    float l = sqrtf(u*u + v*v);
    u /= l;
    v /= l;
  }
  float sqrMag() {
    return (u*u + v*v);
  }
  float mag() {
    return sqrtf(u*u + v*v);
  }
  void scale(float su, float sv) {
    u *= su;
    v *= sv;
  }
  void scale(float s) {
    scale(s,s);
  }
  void translate(float du, float dv) {
    u += du;
    v += dv;
  }

  Vec2 operator+ (const Vec2 &rhs) const {
    return {
      u+rhs.u,
      v+rhs.v,
    };
  }
  Vec2 operator- (const Vec2 &rhs) const {
    return {
      u-rhs.u,
      v-rhs.v,
    };
  }

  Vec2 operator* (float rhs) {
    return {u*rhs, v*rhs};
  }
};

struct Vec3 {
  float x = 0.0f,y=0.0f,z=0.0f;

  void translate(float dx, float dy, float dz) {
    x += dx;
    y += dy;
    z += dz;
  }

  float sqrMag() {
    return (x*x + y*y + z*z);
  }
  float mag() {
    return sqrtf(x*x + y*y + z*z);
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
  void scale(float s) {
    scale(s,s,s);
  }

  Vec3 operator+ (const Vec3 &rhs) const {
    return {
      x+rhs.x,
      y+rhs.y,
      z+rhs.z
    };
  }
  Vec3 operator- (const Vec3 &rhs) const {
    return {
      x-rhs.x,
      y-rhs.y,
      z-rhs.z
    };
  }
  float dot(const Vec3 &v) const {
    return x*v.x + y*v.y + z*v.z;
  }
  Vec3 cross(const Vec3 &v) const {
    return {
      y * v.z - z * v.y,
      z * v.x - x * v.z,
      x * v.y - y * v.x,
    };
  }
};

//given point on plane and normal of plane, as well as start and end points of line segment
// plane_n expected to be normalized
//return point that intersects with plane, and also the ratio of where the intersection is between line start and line end
std::pair<Vec3, float> intersectPlane(Vec3 &plane_p, Vec3 &plane_n, Vec3 &lineStart, Vec3 &lineEnd) {
  plane_n.normalize();
  float plane_d = -plane_n.dot(plane_p);
  float ad = lineStart.dot(plane_n);
  float bd = lineEnd.dot(plane_n);
  float t = (-plane_d-ad) / (bd-ad);
  Vec3 startToEnd = lineEnd - lineStart;
  Vec3 intersect = startToEnd;
  intersect.scale(t);
  return {lineStart + intersect, t};
}

struct Triangle {
  std::array<Vec3,3> points;
  std::array<Vec2,3> uvs;
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

  Triangle copy() {
    return {
      .points = points,
      .uvs = uvs,
      .color = color,
    };
  }

  // return how many of the two output tris to use
  int clipAgainstPlane(Vec3 plane_p, Vec3 plane_n, Triangle &out1, Triangle &out2) {
    //plane_n.normalize();

    auto dist = [&](Vec3 &p) {
      return (
        plane_n.x*p.x + plane_n.y*p.y + plane_n.z*p.z - plane_n.dot(plane_p)
      );
    };

    std::array<Vec3*, 3> insidePts;
    int insidePtCount = 0;
    std::array<Vec3*, 3> outsidePts;
    int outsidePtCount = 0;

    std::array<Vec2*, 3> insideUvs; //texture count
    int insideUvCount = 0;
    std::array<Vec2*, 3> outsideUvs;
    int outsideUvCount = 0;

    float d0 = dist(points[0]);
    float d1 = dist(points[1]);
    float d2 = dist(points[2]);
    if (d0 >= 0) {
      insidePts[insidePtCount++] = &points[0];
      insideUvs[insideUvCount++] = &uvs[0];
    }
    else {
      outsidePts[outsidePtCount++] = &points[0];
      outsideUvs[outsideUvCount++] = &uvs[0];
    }
    if (d1 >= 0) {
      insidePts[insidePtCount++] = &points[1];
      insideUvs[insideUvCount++] = &uvs[1];
    }
    else {
      outsidePts[outsidePtCount++] = &points[1];
      outsideUvs[outsideUvCount++] = &uvs[1];
    }
    if (d2 >= 0) {
      insidePts[insidePtCount++] = &points[2];
      insideUvs[insideUvCount++] = &uvs[2];
    }
    else {
      outsidePts[outsidePtCount++] = &points[2];
      outsideUvs[outsideUvCount++] = &uvs[2];
    }

    if (insidePtCount == 0) {
      // all points are outside plane, so whole triangle should be clipped
      return 0;
    }
    else if (insidePtCount == 3) {
      // all points are inside plane, so triangle should be rendered whole and not clipped
      out1 = copy();
      return 1;
    }
    else if (insidePtCount == 1 && outsidePtCount == 2) {
      // triangle should be clipped, resulting in a smaller triangle
      out1.color = color;
      out1.points[0] = *(insidePts[0]);
      out1.uvs[0] = *(insideUvs[0]);
      float t;
      std::tie(out1.points[1],t) = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[0]));
      out1.uvs[1] = (*(outsideUvs[0]) - *(insideUvs[0]))*t + *(insideUvs[0]);

      std::tie(out1.points[2],t) = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[1]));
      out1.uvs[2] = (*(outsideUvs[1]) - *(insideUvs[0]))*t + *(insideUvs[0]);

      return 1;
    }
    else if (insidePtCount == 2 && outsidePtCount == 1) {
      // triangle should be clipped, resulting in quad (two tris)
      out1.color = color;
      out2.color = color;

      out1.points[0] = *(insidePts[0]);
      out1.points[1] = *(insidePts[1]);
      out1.uvs[0] = *(insideUvs[0]);
      out1.uvs[1] = *(insideUvs[1]);

      float t;
      std::tie(out1.points[2], t) = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[0]));
      out1.uvs[2] = (*(outsideUvs[0]) - *(insideUvs[0]))*t + *(insideUvs[0]);

      out2.points[0] = *(insidePts[1]);
      out2.points[1] = out1.points[2];
      out2.uvs[0] = *(insideUvs[1]);
      out2.uvs[1] = out1.uvs[2];
      std::tie(out2.points[2], t) = intersectPlane(plane_p, plane_n, *(insidePts[1]), *(outsidePts[0]));
      out2.uvs[2] = (*(outsideUvs[0]) - *(insideUvs[1]))*t + *(insideUvs[1]);

      return 2;
    }
    //should be unreachable
    return -1;
  }

};

struct Mesh {
  std::vector<Triangle> tris;

  // Load from .obj file
  bool loadObj(std::string fp) {
    std::ifstream f(fp);
    if (!f.is_open()) {
      std::cout << "Could not load object file: " << fp << std::endl;
      return false;
    }

    tris.clear();
    std::vector<Vec3> verts;
    std::string line {};
    while (getline(f, line)) {
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

const Mesh CUBE = {.tris = {
  // SOUTH
  { 0.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f,  1.0f, 1.0f, 0.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,}, 
  { 0.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,  1.0f, 0.0f, 0.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},

  // EAST
  { 1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,  1.0f, 1.0f, 1.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,},
  { 1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},

  // NORTH
  { 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,},
  { 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},

  // WEST
  { 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,},
  { 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f,  0.0f, 0.0f, 0.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},

  // TOP
  { 0.0f, 1.0f, 0.0f,  0.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,},
  { 0.0f, 1.0f, 0.0f,  1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 0.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},

  // BOTTOM
  { 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f,    0.0f, 1.0f,  0.0f, 0.0f,  1.0f, 0.0f,},
  { 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f,  1.0f, 0.0f, 0.0f,    0.0f, 1.0f,  1.0f, 0.0f,  1.0f, 1.0f,},
}};

// 4x4 Matrix
struct Mat4 {
  std::array<std::array<float,4>, 4> data = {0};

  float get(int row, int col) const {
    return data[row][col];
  }
  void set(int row, int col, float n) {
    data[row][col] = n;
  }

  Mat4 operator*(const Mat4 &rhs) const {
    Mat4 res;
    for (int col = 0; col < 4; col++) {
      for (int row = 0; row < 4; row++) {
        float x = get(row,0)*rhs.get(0,col) + get(row,1)*rhs.get(1,col) + get(row,2)*rhs.get(2,col) + get(row,3)*rhs.get(3,col);
        res.set(row,col, x);
      }
    }
    return res;
  }
};

const Mat4 identity4 {
  1.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 1.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 1.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 1.0f,
};

// Multiply a 3d Vector by the Projection Matrix by first pretending it's a 4d vector
Vec3 mulProj(const Vec3 &u, const Mat4 &m) {
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

// Pretend a 3d vector is a 4d vector with w=1, multiply by the matrix, then return the result, discarding w
Vec3 mulV4M4(const Vec3 &u, const Mat4 &m) {
  Vec3 res;

  res.x = u.x * m.get(0,0) + u.y * m.get(1,0) + u.z * m.get(2,0) + m.get(3,0);
  res.y = u.x * m.get(0,1) + u.y * m.get(1,1) + u.z * m.get(2,1) + m.get(3,1);
  res.z = u.x * m.get(0,2) + u.y * m.get(1,2) + u.z * m.get(2,2) + m.get(3,2);

  return res;
}

Triangle project(const Triangle &tri, const Mat4 &m) {
  Triangle res;
  res.points[0] = mulProj(tri.points[0], m);
  res.points[1] = mulProj(tri.points[1], m);
  res.points[2] = mulProj(tri.points[2], m);
  res.color = tri.color;
  res.uvs = tri.uvs;
  return res;
}

olc::Pixel getColor(float lum) {
  return olc::WHITE * std::max((lum+1.0f)/2.0f, 0.1f);
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

Mat4 matPointAt(const Vec3 &pos, const Vec3 &target, const Vec3 &up) {
  Vec3 newForward = target-pos;
  newForward.normalize();

  Vec3 a, newUp;
  a = newForward;
  a.scale(up.dot(newForward));
  newUp = up-a;
  newUp.normalize();

  Vec3 newRight = newUp.cross(newForward);

  return {
    newRight.x,   newRight.y,   newRight.z,   0.0f,
    newUp.x,      newUp.y,      newUp.z,      0.0f,
    newForward.x, newForward.y, newForward.z, 0.0f, 
    pos.x,        pos.y,        pos.z,        1.0f,
  };
}

// Only works for rotation/translation/point-at matrices
Mat4 quickInverse(Mat4 &m) {
  Mat4 res;
  res.set(0,0, m.get(0,0));
  res.set(0,1, m.get(1,0));
  res.set(0,2, m.get(2,0));
  res.set(0,3, 0.0f);
  res.set(1,0, m.get(0,1));
  res.set(1,1, m.get(1,1));
  res.set(1,2, m.get(2,1));
  res.set(1,3, 0.0f);
  res.set(2,0, m.get(0,2));
  res.set(2,1, m.get(1,2));
  res.set(2,2, m.get(2,2));
  res.set(2,3, 0.0f);
  res.set(3,0, -(m.get(3,0) * res.get(0,0) + m.get(3,1) * res.get(1,0) + m.get(3,2) * res.get(2,0)));
  res.set(3,1, -(m.get(3,0) * res.get(0,1) + m.get(3,1) * res.get(1,1) + m.get(3,2) * res.get(2,1)));
  res.set(3,2, -(m.get(3,0) * res.get(0,2) + m.get(3,1) * res.get(1,2) + m.get(3,2) * res.get(2,2)));
  res.set(3,3, 1.0f);
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

Mat4 matRotY(float theta) {
  Mat4 res;
  res.set(0,0, cosf(theta));
  res.set(0,2, sinf(theta));
  res.set(2,0, -sinf(theta));
  res.set(1,1, 1.0f);
  res.set(2,2, cosf(theta));
  res.set(3,3, 1.0f);
  return res;
}

Mat4 matTrans(float x, float y, float z) {
  Mat4 res;
  res.set(0,0, 1.0f);
  res.set(1,1, 1.0f);
  res.set(2,2, 1.0f);
  res.set(3,3, 1.0f);
  res.set(3,0, x);
  res.set(3,1, y);
  res.set(3,2, z);
  return res;
}

class Game : public olc::PixelGameEngine {
  Mesh meshCube;
  olc::Sprite spr;
  Mat4 projMat;
  Vec3 camPos = {0};
  Vec3 camDir {0,0,1};
  float yaw;

  float theta;

public:
  Game() {
    sAppName = "Nomad Space";
  }

  bool OnUserCreate() override {
    /*
    if(!meshCube.loadObj("assets/axis.obj")) {
      return false;
    }
    */
    meshCube = CUBE;
    spr = olc::Sprite("assets/tex.png");

    projMat = projectionMatrix(0.1f, 1000.0f, 90.0f, static_cast<float>(ScreenHeight()) / static_cast<float>(ScreenWidth()));
    return true;
  }

  bool OnUserUpdate(float dt) override {
    // Update
    if (GetKey(olc::Key::UP).bHeld) {
      camPos.y += 8.0f * dt;
    }
    if (GetKey(olc::Key::DOWN).bHeld) {
      camPos.y -= 8.0f * dt;
    }
    if (GetKey(olc::Key::LEFT).bHeld) {
      camPos.x += 8.0f * dt;
    }
    if (GetKey(olc::Key::RIGHT).bHeld) {
      camPos.x -= 8.0f * dt;
    }
    if (GetKey(olc::Key::A).bHeld) {
      yaw -= 2.0f * dt;
    }
    if (GetKey(olc::Key::D).bHeld) {
      yaw += 2.0f * dt;
    }
    Vec3 moveForward = camDir;
    moveForward.scale(8.0f * dt);
    if (GetKey(olc::Key::W).bHeld) {
      camPos = camPos + moveForward;
    }
    if (GetKey(olc::Key::S).bHeld) {
      camPos = camPos - moveForward;
    }

    // Draw
    Clear(olc::BLACK);

    Mat4 rotZ, rotX;
    //theta += 1.0f * dt;
    rotZ = matRotZ(theta);
    rotX = matRotX(theta * 0.5f);

    Mat4 transM = matTrans(0,0, 16.0f);

    Mat4 matWorld = identity4;
    matWorld = rotZ * rotX;
    matWorld = matWorld * transM;

    Vec3 camUp {0,1,0};
    Vec3 camTarget {0,0,1};
    Mat4 camRot = matRotY(yaw);
    camDir = mulV4M4(camTarget, camRot);
    camTarget = camPos + camDir;
    Mat4 camMat = matPointAt(camPos, camTarget, camUp);
    Mat4 viewMat = quickInverse(camMat);

    std::vector<Triangle> trisToDraw;

    // Calculate triangle projection and add to draw queue
    for(auto &tri: meshCube.tris) {
      Triangle triProj, triTrans, triViewed; //projected and transformed triangles

      triTrans = project(tri, matWorld);
      //triTrans.uvs = tri.uvs;

      Vec3 normal, u,v;
      u = triTrans.points[1] - triTrans.points[0];
      v = triTrans.points[2] - triTrans.points[0];
      normal = u.cross(v);
      normal.normalize();

      if (normal.dot(triTrans.points[0]-camPos) < 0.0f) {
        Vec3 light_dir {0.0f, 0.0f, -1.0f};
        light_dir.normalize();

        float dp = normal.dot(light_dir);
        triTrans.color = getColor(dp);

        triViewed = project(triTrans, viewMat);
        //triViewed.uvs = tri.uvs;

        std::array<Triangle, 2> outTris;
        int clippedTris = triViewed.clipAgainstPlane({0,0,0.1f}, {0,0,1.0f}, outTris[0], outTris[1]);

        for(int i = 0; i < clippedTris; i++) {
          triProj = project(outTris[i], projMat);
          //triProj.uvs = outTris[i].uvs;

          triProj.scale(-1,-1,1);
          triProj.translate(1.0f, 1.0f, 0.0f);
          triProj.scale(0.5f * static_cast<float>(ScreenWidth()), 0.5f * static_cast<float>(ScreenHeight()), 1.0f);

          trisToDraw.push_back(triProj);
        }
      }
    }

    // Sort triangles from back to front
    sort(trisToDraw.begin(), trisToDraw.end(), [](Triangle &t1, Triangle &t2) {
      float avg1 = (t1.points[0].z + t1.points[1].z + t1.points[2].z) / 3.0f;
      float avg2 = (t2.points[0].z + t2.points[1].z + t2.points[2].z) / 3.0f;
      return avg1 > avg2;
    });

    // Draw sorted triangles
    for (auto &tri: trisToDraw) {
      // clip triangles against screen edges
      std::array<Triangle, 2> clipped;
      std::list<Triangle> triQueue;
      triQueue.push_back(tri);
      int newTris = 1;

      for(int p = 0; p < 4; p++) {
        int trisToAdd = 0;
        while (newTris > 0) {
          Triangle test = triQueue.front();
          triQueue.pop_front();
          newTris--;

          switch (p) {
          case 0: trisToAdd = test.clipAgainstPlane({0,0,0}, {0,1.0f,0}, clipped[0], clipped[1]); break;
          case 1: trisToAdd = test.clipAgainstPlane({0,static_cast<float>(ScreenHeight())-1.0f,0}, {0,-1.0f,0}, clipped[0], clipped[1]); break;
          case 2: trisToAdd = test.clipAgainstPlane({0,0,0}, {1.0f,0,0}, clipped[0], clipped[1]); break;
          case 3: trisToAdd = test.clipAgainstPlane({static_cast<float>(ScreenWidth())-1.0f,0,0}, {-1.0f,0,0}, clipped[0], clipped[1]); break;
          }
          for (int w = 0; w < trisToAdd; w++) {
            triQueue.push_back(clipped[w]);
          }
        }
        newTris = triQueue.size();
      }

      for (auto &t : triQueue) {
        TexturedTriangle(
          t.points[0].x, t.points[0].y,  t.uvs[0].u, t.uvs[0].v,
          t.points[1].x, t.points[1].y,  t.uvs[1].u, t.uvs[1].v,
          t.points[2].x, t.points[2].y,  t.uvs[2].u, t.uvs[2].v,
          spr
        );
        /*
        FillTriangle(
          t.points[0].x, t.points[0].y,
          t.points[1].x, t.points[1].y,
          t.points[2].x, t.points[2].y,
          t.color
        );
        */
        // Draw wireframe
        DrawTriangle(
          t.points[0].x, t.points[0].y,
          t.points[1].x, t.points[1].y,
          t.points[2].x, t.points[2].y,
          olc::WHITE
        );
      }
    }
    DrawSprite(0,0, &spr);
    return true;
  }

  void TexturedTriangle(int x1, int y1, float u1, float v1,
                        int x2, int y2, float u2, float v2,
                        int x3, int y3, float u3, float v3,
                        const olc::Sprite &spr) {
    if (y2 < y1) {
      std::swap(y1,y2);
      std::swap(x1,x2);
      std::swap(u1,u2);
      std::swap(v1,v2);
    }
    if (y3 < y1) {
      std::swap(y1,y3);
      std::swap(x1,x3);
      std::swap(u1,u3);
      std::swap(v1,v3);
    }
    if (y3 < y2) {
      std::swap(y3,y2);
      std::swap(x3,x2);
      std::swap(u3,u2);
      std::swap(v3,v2);
    }
    int dy1 = y2-y1;
    int dx1 = x2-x1;
    float du1 = u2-u1;
    float dv1 = v2-v1;

    int dy2 = y3-y1;
    int dx2 = x3-x1;
    float du2 = u3-u1;
    float dv2 = v3-v1;

    float dax_step = 0, dbx_step = 0,
      du1_step = 0, dv1_step = 0,
      du2_step = 0, dv2_step = 0;

    if (dy1) {
      dax_step = dx1 / static_cast<float>(abs(dy1));
      du1_step = du1 / static_cast<float>(abs(dy1));
      dv1_step = dv1 / static_cast<float>(abs(dy1));
    }
    if (dy2) {
      dbx_step = dx2 / static_cast<float>(abs(dy2));
      du2_step = du2 / static_cast<float>(abs(dy2));
      dv2_step = dv2 / static_cast<float>(abs(dy2));
    }

    if (dy1) {
      float tu, tv; //final uv points
      for (int i=y1; i <= y2; i++) {
        int ax = x1 + static_cast<float>(i-y1) * dax_step;
        int bx = x1 + static_cast<float>(i-y1) * dbx_step;

        //starting uv coord
        float su = u1 + static_cast<float>(i-y1) * du1_step;
        float sv = v1 + static_cast<float>(i-y1) * dv1_step;

        //ending uv coord
        float eu = u1 + static_cast<float>(i-y1) * du2_step;
        float ev = v1 + static_cast<float>(i-y1) * dv2_step;

        if (ax > bx) {
          std::swap(ax,bx);
          std::swap(su,eu);
          std::swap(sv,ev);
        }
        tu = su;
        tv = sv;
        float t_step = 1.0f / static_cast<float>(bx - ax);
        float t = 0.0f;

        for (int j = ax; j < bx; j++) {
          tu = (1.0f - t)*su + t*eu;
          tv = (1.0f - t)*sv + t*ev;
          t += t_step;
          int w = static_cast<float>(spr.width-1) * tu;
          int h = static_cast<float>(spr.height-1) * tv;
          //std::cout << tu << ", " << tv << "  --  " << w << ", " << h << std::endl;
          Draw(j,i,spr.GetPixel(w,h));
        }
      }
      dy1 = y3 - y2;
      dx1 = x3 - x2;
      du1 = u3 - u2;
      dv1 = v3 - v2;
      if (dy1) {
        dax_step = dx1 / static_cast<float>(abs(dy1));
        du1_step = du1 / static_cast<float>(abs(dy1));
        dv1_step = dv1 / static_cast<float>(abs(dy1));
      }
      else {
        du1_step = 0;
        dv1_step = 0;
      }
      if (dy2) {
        dbx_step = dx2 / static_cast<float>(abs(dy2));
      }

      for (int i=y2; i <= y3; i++) {
        int ax = x2 + static_cast<float>(i-y2) * dax_step;
        int bx = x1 + static_cast<float>(i-y1) * dbx_step;

        //starting uv coord
        float su = u2 + static_cast<float>(i-y2) * du1_step;
        float sv = v2 + static_cast<float>(i-y2) * dv1_step;

        //ending uv coord
        float eu = u1 + static_cast<float>(i-y1) * du2_step;
        float ev = v1 + static_cast<float>(i-y1) * dv2_step;

        if (ax > bx) {
          std::swap(ax,bx);
          std::swap(su,eu);
          std::swap(sv,ev);
        }
        tu = su;
        tv = sv;
        float t_step = 1.0f / static_cast<float>(bx - ax);
        float t = 0.0f;

        for (int j = ax; j < bx; j++) {
          tu = (1.0f - t)*su + t*eu;
          tv = (1.0f - t)*sv + t*ev;
          t += t_step;
          int w = static_cast<float>(spr.width-1) * tu;
          int h = static_cast<float>(spr.height-1) * tv;
          //std::cout << tu << ", " << tv << "  --  " << w << ", " << h << std::endl;
          Draw(j,i,spr.GetPixel(w,h));
        }
        //DrawRect(ax-2,i-2,4,4);
        //DrawRect(bx-2,i-2,4,4);
      }
    }
  }
};

int main() {
  Game demo;
  if (demo.Construct(320,200,2,2,false, true)) {
    demo.Start();
  }
  return 0;
}
