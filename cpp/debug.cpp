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

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <mpi.h>
#include "real.h"
#include "debug.hpp"
#include "solfec.hpp"
#include "compute.hpp"

/* https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf */
template<typename ... Args>
std::string sformat (const std::string& format, Args ... args)
{
  size_t size = snprintf (nullptr, 0, format.c_str(), args ... ) + 1; /* Extra space for '\0' */
  std::unique_ptr<char[]> buf(new char[size]); 
  snprintf (buf.get(), size, format.c_str(), args ...);
  return std::string (buf.get(), buf.get() + size - 1); /* We don't want the '\0' inside */
}

/* output compute data into per rank files */
void debug_output_compute_data()
{
  using namespace compute;
  std::ofstream out;
  std::string path;
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  /* materials */
  {
    path = sformat ("debug_materials%d.txt", rank);
    out.open(path);
    out << "TIME " << solfec::simulation_time << std::endl;

    uint64_t count;
    ga_counters->get(rank, cn_materials, cn_materials+1, 0, 1, &count);
    out << "COUNT " << count << std::endl;

    REAL *data = new REAL[count*mt_last];
    ga_materials->get(rank, 0, count, 0, mt_last, data);

    for (uint64_t i = 0; i < count; i ++)
    {
      out << "DENSITY " << data[mt_density*count+i] << std::endl;
      out << "YOUNG " << data[mt_young*count+i] << std::endl;
      out << "POISSON " << data[mt_poisson*count+i] << std::endl;
      out << "VISCOSITY " << data[mt_viscosity*count+i] << std::endl;
    }

    delete [] data;
    out.close();
  }

  /* nodes */
  {
    path = sformat ("debug_nodes%d.txt", rank);
    out.open(path);
    out << "TIME " << solfec::simulation_time << std::endl;

    uint64_t count;
    ga_counters->get(rank, cn_nodes, cn_nodes+1, 0, 1, &count);

    REAL *data = new REAL[count*4];
    ga_nodes->get(rank, 0, count, nd_X, nd_unused+1, data);

    REAL *unused = data+3*count;
    uint64_t used = 0;
    for (uint64_t i = 0; i < count; i ++)
      if (unused[i] == 0.) used ++;

    out << "COUNT " << used << std::endl;

    for (uint64_t i = 0; i < count; i ++)
    {
      if (unused[i] == 0.)
      {
        out << i << " "
            << data[0*count+i] << " "
            << data[1*count+i] << " "
            << data[2*count+i] << std::endl;
      }
    }

    delete [] data;
    out.close();
  }

  /* elements */
  {
    path = sformat ("debug_elements%d.txt", rank);
    out.open(path);
    out << "TIME " << solfec::simulation_time << std::endl;

    uint64_t count;
    ga_counters->get(rank, cn_elements, cn_elements+1, 0, 1, &count);
    out << "COUNT " << count << std::endl;

    uint64_t *data = new uint64_t[count*el_last];
    ga_elements->get(rank, 0, count, 0, el_last, data);

    for (uint64_t i = 0; i < count; i ++)
    {
      uint64_t type = data[el_type*count+i];
      out << "BODNUM " << data[el_bodnum*count+i] << std::endl;
      out << "MATNUM " << data[el_matnum*count+i] << std::endl;
      out << "TYPE " << type << std::endl;
      out << data[el_nd0_rnk*count+i] << " " << data[el_nd0_idx*count+i] << std::endl;
      out << data[el_nd1_rnk*count+i] << " " << data[el_nd1_idx*count+i] << std::endl;
      out << data[el_nd2_rnk*count+i] << " " << data[el_nd2_idx*count+i] << std::endl;
      out << data[el_nd3_rnk*count+i] << " " << data[el_nd3_idx*count+i] << std::endl;
      if (type >= 5) {
	out << data[el_nd4_rnk*count+i] << " " << data[el_nd4_idx*count+i] << std::endl; }
      if (type >= 6) {
	out << data[el_nd5_rnk*count+i] << " " << data[el_nd5_idx*count+i] << std::endl; }
      if (type >= 8) {
	out << data[el_nd6_rnk*count+i] << " " << data[el_nd6_idx*count+i] << std::endl;
	out << data[el_nd7_rnk*count+i] << " " << data[el_nd7_idx*count+i] << std::endl; }
    }

    delete [] data;
    out.close();
  }

  /* ellipsoids  */
  {
    path = sformat ("debug_ellips%d.txt", rank);
    out.open(path);
    out << "TIME " << solfec::simulation_time << std::endl;

    uint64_t count;
    ga_counters->get(rank, cn_ellips, cn_ellips+1, 0, 1, &count);
    out << "COUNT " << count << std::endl;

    REAL *data0 = new REAL[count*ll_last0];
    uint64_t *data1 = new uint64_t[count*ll_last1];
    ga_elldata->get(rank, 0, count, 0, ll_last0, data0);
    ga_ellips->get(rank, 0, count, 0, ll_last1, data1);

    for (uint64_t i = 0; i < count; i ++)
    {
      out << "BODNUM " << data1[ll_bodnum*count+i] << std::endl;
      out << "MATNUM " << data1[ll_matnum*count+i] << std::endl;
      out << "COLOR " << data1[ll_color*count+i] << std::endl;
      out << "CENTER " << data0[ll_x*count+i] << " " << data0[ll_y*count+i] << " " << data0[ll_z*count+i] << std::endl;
      out << "RADIUS " << data0[ll_a*count+i] << " " << data0[ll_b*count+i] << " " << data0[ll_c*count+i] << std::endl;
      out << "ROTATION " << data0[ll_r0*count+i] << " " << data0[ll_r1*count+i] << " " << data0[ll_r2*count+i] << " "
                         << data0[ll_r3*count+i] << " " << data0[ll_r4*count+i] << " " << data0[ll_r5*count+i] << " "
                         << data0[ll_r6*count+i] << " " << data0[ll_r7*count+i] << " " << data0[ll_r8*count+i] << std::endl;
    }

    delete [] data0;
    delete [] data1;
    out.close();
  }
}
