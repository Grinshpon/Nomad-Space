#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <strstream>
#include "olcPixelGameEngine.h"
//import linalg;

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

// plane_n expected to be normalized
Vec3 intersectPlane(Vec3 &plane_p, Vec3 &plane_n, Vec3 &lineStart, Vec3 &lineEnd) {
  plane_n.normalize();
  float plane_d = -plane_n.dot(plane_p);
  float ad = lineStart.dot(plane_n);
  float bd = lineEnd.dot(plane_n);
  float t = (-plane_d-ad) / (bd-ad);
  Vec3 startToEnd = lineEnd - lineStart;
  Vec3 intersect = startToEnd;
  intersect.scale(t);
  return lineStart + intersect;
}

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

  Triangle copy() {
    return {
      .points = points,
      .color = color
    };
  }

  // return how many of the two output tris to use
  int clipAgainstPlane(Vec3 plane_p, Vec3 plane_n, Triangle &out1, Triangle &out2) {
    plane_n.normalize();

    auto dist = [&](Vec3 &p) {
      return (
        plane_n.x*p.x + plane_n.y*p.y + plane_n.z*p.z - plane_n.dot(plane_p)
      );
    };

    std::array<Vec3*, 3> insidePts;
    int insidePtCount = 0;
    std::array<Vec3*, 3> outsidePts;
    int outsidePtCount = 0;

    float d0 = dist(points[0]);
    float d1 = dist(points[1]);
    float d2 = dist(points[2]);
    if (d0 >= 0) { insidePts[insidePtCount++] = &points[0];}
    else { outsidePts[outsidePtCount++] = &points[0];}
    if (d1 >= 0) { insidePts[insidePtCount++] = &points[1];}
    else { outsidePts[outsidePtCount++] = &points[1];}
    if (d2 >= 0) { insidePts[insidePtCount++] = &points[2];}
    else { outsidePts[outsidePtCount++] = &points[2];}

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
      out1.points[1] = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[0]));
      out1.points[2] = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[1]));
      return 1;
    }
    else if (insidePtCount == 2 && outsidePtCount == 1) {
      // triangle should be clipped, resulting in quad (two tris)
      out1.color = color;
      out2.color = color;

      out1.points[0] = *(insidePts[0]);
      out1.points[1] = *(insidePts[1]);
      out1.points[2] = intersectPlane(plane_p, plane_n, *(insidePts[0]), *(outsidePts[0]));
      out2.points[0] = *(insidePts[1]);
      out2.points[1] = out1.points[2];
      out2.points[2] = intersectPlane(plane_p, plane_n, *(insidePts[1]), *(outsidePts[0]));
      return 2;
    }
    return -1;
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

  Mat4 operator*(Mat4 &rhs) {
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

const Mat4 identity4 = {
  1.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 1.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 1.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 1.0f,
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

// Pretend a 3d vector is a 4d vector with w=1, multiply by the matrix, then return the result, discarding w
Vec3 mulV4M4(Vec3 &u, Mat4 &m) {
  Vec3 res;

  res.x = u.x * m.get(0,0) + u.y * m.get(1,0) + u.z * m.get(2,0) + m.get(3,0);
  res.y = u.x * m.get(0,1) + u.y * m.get(1,1) + u.z * m.get(2,1) + m.get(3,1);
  res.z = u.x * m.get(0,2) + u.y * m.get(1,2) + u.z * m.get(2,2) + m.get(3,2);

  return res;
}

Triangle project(Triangle &tri, Mat4 &m) {
  Triangle res;
  res.points[0] = mulProj(tri.points[0], m);
  res.points[1] = mulProj(tri.points[1], m);
  res.points[2] = mulProj(tri.points[2], m);
  res.color = tri.color;
  return res;
}

olc::Pixel getColor(float lum) {
  return olc::WHITE * std::max(lum, 0.1f);
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

Mat4 matPointAt(Vec3 &pos, Vec3 &target, Vec3 &up) {
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
  Mat4 projMat;
  Vec3 camPos = {0};
  Vec3 camDir = {0,0,1};
  float yaw;

  float theta;

public:
  Game() {
    sAppName = "Test";
  }

  bool OnUserCreate() override {
    if(!meshCube.loadObj("assets/axis.obj")) {
      std::cout << "Could not load object file\n";
    }

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

    Vec3 camUp = {0,1,0};
    Vec3 camTarget = {0,0,1};
    Mat4 camRot = matRotY(yaw);
    camDir = mulV4M4(camTarget, camRot);
    camTarget = camPos + camDir;
    Mat4 camMat = matPointAt(camPos, camTarget, camUp);
    Mat4 viewMat = quickInverse(camMat);

    std::vector<Triangle> trisToDraw;

    // Calculate triangle projection and add to draw queue
    for(auto tri: meshCube.tris) {
      //tri.translate(10,0,0);
      Triangle triProj, triTrans, triViewed; //projected and transformed triangles

      triTrans = project(tri, matWorld);

      Vec3 normal, u,v;
      u = triTrans.points[1] - triTrans.points[0];
      v = triTrans.points[2] - triTrans.points[0];
      normal = u.cross(v);
      normal.normalize();

      if (normal.dot(triTrans.points[0]-camPos) < 0.0f) {
        Vec3 light_dir = {0.0f, 0.0f, -1.0f};
        light_dir.normalize();

        float dp = normal.dot(light_dir);
        triTrans.color = getColor(dp);

        triViewed = project(triTrans, viewMat);

        std::array<Triangle, 2> outTris;
        int clippedTris = triViewed.clipAgainstPlane({0,0,0.1f}, {0,0,1.0f}, outTris[0], outTris[1]);

        for(int i = 0; i < clippedTris; i++) {
          triProj = project(outTris[i], projMat);
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
        FillTriangle(
          t.points[0].x, t.points[0].y,
          t.points[1].x, t.points[1].y,
          t.points[2].x, t.points[2].y,
          t.color
        );
        // Draw wireframe
        DrawTriangle(
          t.points[0].x, t.points[0].y,
          t.points[1].x, t.points[1].y,
          t.points[2].x, t.points[2].y,
          olc::BLACK
        );
      }
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
