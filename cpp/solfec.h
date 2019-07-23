/*
The MIT License (MIT)

Copyright (c) 2019 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>

#ifndef __solfec__
#define __solfec__

/* spline */
struct spline
{
  std::vector<std::array<REAL,2>> points; /* vector of (x, y) pairs */
  size_t cache; /* partial cache size */
  std::string path; /* file path for partially cached splines */
  std::vector<off_t> offset; /* file offsets */
  std::vector<REAL> xval; /* x values maching file offsets */
  size_t marker; /* index of the last read interval */
  spline(): marker(0), cache(0) {}
};

/* material */
struct material
{
  REAL density; /* material density */
  REAL young; /* Young's modulus */
  REAL poisson; /* Poisson's ratio */
  REAL viscosity; /* viscosisty coefficient */
};

/* mesh */
struct mesh
{
  std::array<std::vector<REAL>,3> nodes; /* list of nodes */
  std::vector<size_t> elements; /* list of elements */
  std::vector<size_t> colors; /* face colors */
  size_t nhex, nwed, npyr, ntet; /* numbers of elements per kind */
  size_t nfaces; /* number of faces */
  size_t matnum; /* material number */
  size_t gcolor; /* global color */
  std::set<size_t> restrains; /* applied restrains */
  std::set<size_t> prescribes; /* applied prescribes */
  mesh () : nhex(0), nwed(0), npyr(0), ntet(0), nfaces(0), matnum(0), gcolor(0) { }
};

/* restrain */
struct restrain
{
  size_t bodnum; /* body number */
  std::vector<std::array<REAL,3>> points; /* list of points */
  size_t color; /* surface color */
  restrain() : bodnum(0), color(0) { }
};

/* prescribe */
struct prescribe
{
  size_t bodnum; /* body number */
  std::vector<std::array<REAL,3>> points; /* list of points */
  size_t color; /* surface color */
  bool linear_applied; /* linear applied flag */
  REAL linear_values [3]; /* constant linear velocity */
  int64_t linear_splines [3]; /* linear splines (when >= 0) */
  void* linear_callback; /* Python callback function */
  bool angular_applied;
  REAL angular_values [3];
  int64_t angular_splines [3];
  void* angular_callback;
  prescribe() : bodnum(0), color(0), linear_applied(false), angular_applied(false) { }
};

/* velocity */
struct velocity
{
  size_t bodnum; /* body number */
  REAL linear_values [3]; /* linear velocity */
  REAL angular_values [3]; /* linear velocity */
};

/* friction */
struct friction
{
  size_t color1; /* first face color */
  size_t color2; /* second face color */
  REAL static_friction; /* static Coulomb's friction */
  REAL dynamic_friction;/* dynamic Coulomb's friction */
};

/* friction comparison */
struct friction_compare {
  bool operator() (const struct friction& lhs, const struct friction& rhs) const
  {return lhs.color1 == rhs.color1 ? lhs.color2 < rhs.color2 : lhs.color1 < lhs.color2;}
};

/* gravity */
struct gravity
{
  REAL gvalue [3]; /* constant values */
  int64_t gspline [3]; /* linear splines (when >= 0) */
  void* gcallback [3]; /* Python callback functions */
  gravity(): gvalue{0,0,0}, gspline{-1,-1,-1}, gcallback{NULL,NULL,NULL} {}
  void clear() {gvalue[0]=gvalue[1],gvalue[2]=0.;gspline[0]=gspline[1]=gspline[1]=-1;gcallback[0]=gcallback[1]=gcallback[2]=NULL;}
};

/* history */
struct history
{
  std::string entity; /* entity name */
  REAL point [3]; /* referential point */
  size_t bodnum; /* body number */
  std::string filepath; /* text file path */
  std::vector<REAL> values; /* stored values */
};

/* output */
struct output
{
  std::vector<std::string> entities; /* output entities */
  std::vector<size_t> subset; /* subset of body numbers */
  std::set<std::string> modes; /* output modes */
  std::set<std::string> formats; /* output formats */
};

/* solfec global variables */
namespace solfec
{
extern int argc;
extern char **argv;
extern std::string outname;
extern std::map<size_t,spline> splines;
extern size_t splines_count;
extern std::map<size_t,material> materials;
extern size_t materials_count;
extern std::map<size_t,mesh> meshes;
extern size_t meshes_count;
extern std::map<size_t,restrain> restrains;
extern size_t restrains_count;
extern std::map<size_t,prescribe> prescribes;
extern size_t prescribes_count;
extern std::vector<velocity> velocities;
extern std::set<friction,friction_compare> frictions;
extern struct gravity gravity;
extern std::vector<history> histories;
extern std::vector<output> outputs;
extern bool notrun;
};

/* read spline from file */
void spline_from_file (const char *path, int cache, struct spline *spline);

/* calculate splane value */
REAL spline_value (struct spline *spline, REAL xval);

/* interpret an input file (return 0 on success) */
int input (const char *path);

#endif
