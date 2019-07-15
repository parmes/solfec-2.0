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
  spline(): marker(0) {}
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
  std::vector<std::array<REAL,3>> nodes; /* list of nodes */
  std::vector<size_t> elements; /* list of elements */
  size_t matnum; /* material number */
  std::vector<size_t> colors; /* face colors */
};

/* friction */
struct friction
{
  size_t color1; /* first face color */
  size_t color2; /* second face color */
  REAL static_friction; /* static Coulomb's friction */
  REAL dynamic_friction;/* dynamic Coulomb's friction */
};

/* restrain */
struct restrain
{
  size_t bodnum; /* body number */
  std::vector<std::array<REAL,3>> points; /* list of points */
  size_t size; /* list size */
  size_t color; /* surface color */
};

/* prescribe */
struct prescribe
{
  size_t bodnum; /* body number */
  std::vector<std::array<REAL,3>> points; /* list of points */
  size_t size; /* list size */
  size_t color; /* surface color */
  REAL linear_values [3]; /* constant linear velocity */
  int64_t linear_splines [3]; /* linear splines (when >= 0) */
  void* linear_callback; /* Python callback function */
  REAL angular_values [3];
  int64_t angular_splines [3];
  void* angular_callback;
};

/* velocity */
struct velocity
{
  size_t bodnum; /* body number */
  REAL linear_values [3]; /* linear velocity */
  REAL angular_values [3]; /* linear velocity */
};

/* gravity */
struct gravity
{
  REAL gvalues [3]; /* constant values */
  int64_t gsplines [3]; /* linear splines (when >= 0) */
  void* gcallback [3]; /* Python callback functions */
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
extern std::string output_path;
extern std::vector<spline> splines;
extern std::vector<material> materials;
extern std::vector<mesh> bodies;
extern std::vector<friction> frictions;
extern std::vector<restrain> restrains;
extern std::vector<prescribe> prescribes;
extern std::vector<velocity> velocities;
extern struct gravity gravity;
extern std::vector<history> histories;
extern std::vector<output> outputs;
};

/* read spline from file */
void spline_from_file (char *path, int cache, struct spline *spline);

/* calculate splane value */
REAL spline_value (struct spline *spline, REAL xval);

/* interpret an input file (return 0 on success) */
int input (const char *path);

#endif
