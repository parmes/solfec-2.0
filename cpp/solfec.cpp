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

/* Contributors: Tomasz Koziara */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "real.h"
#include "version.h"
#include "solfec.hpp"
#include "compute.hpp"

/* solfec global variables */
namespace solfec
{
int argc = 0;
char **argv = NULL;
std::string outname;
std::map<uint64_t,spline> splines;
uint64_t splines_count = 0;
std::map<uint64_t, material> materials;
uint64_t materials_count = 0;
std::map<uint64_t,mesh> meshes;
std::map<uint64_t,ellip> ellips;
uint64_t bodies_count = 0;
std::map<uint64_t,std::set<uint64_t>> body_restrains;
std::map<uint64_t, restrain> restrains;
uint64_t restrains_count = 0;
std::map<uint64_t,std::set<uint64_t>> body_prescribes;
std::map<uint64_t,prescribe> prescribes;
uint64_t prescribes_count = 0;
std::map<uint64_t, std::vector<velocity>> velocities;
std::set<friction,friction_compare> frictions;
struct gravity gravity;
std::vector<history> histories;
std::vector<output> outputs;
bool notrun = true;
struct interval interval;
};

int main (int argc, char *argv[])
{
  if (argc == 1)
  {
    printf ("VERSION: 2.%s (%s)\n", VERSION_HASH, VERSION_DATE);
    printf ("SYNOPSIS: [mpirun -np N] solfec%d path/to/input/file.py\n", REALSIZE);
    return 0;
  }

  solfec::argc = argc;
  solfec::argv = argv;

  int rank;

  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) input_python (argv[1]);
  else compute_main_loop();
 
  MPI_Finalize ();

  return 0;
}
