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

#define AMGCL_NO_BOOST
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/lgmres.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <set>
#include <map>
#include <vector>
#include "real.h"
#include "alg.h"
#include "compute.hpp"
#include "timestep.hpp"

#define NELPTSK 32 /* number of elements per task */

/* example application with separated preconditioner and solver;
 * this allows to update the preconditioner less frequently then
 * the system matrix */

typedef amgcl::backend::builtin<REAL> Backend;

typedef amgcl::solver::bicgstab<Backend> Solver;

typedef amgcl::amg<Backend,
  amgcl::coarsening::smoothed_aggregation,
  amgcl::relaxation::spai0> Precond;

typedef std::tuple<uint64_t, /* matrix size */
  std::vector<uint64_t>, /* row pointers */
  std::vector<uint64_t>, /* column indices */
  std::vector<REAL>> Matrix; /* values of a matrix in CRS format */

inline Precond* update_precond (Matrix *matrix, int block_size = 3, int coarse_enough = 3000)
{
  Precond::params pprm;

  pprm.coarsening.aggr.block_size = block_size; /* dense block size */
  pprm.coarse_enough = coarse_enough; /* use direct solve below this size */

  return new Precond(*matrix, pprm);
}

inline std::tuple<int,REAL> solve_system (Matrix *matrix, Precond *precond,
        std::vector<REAL> &b, std::vector<REAL> &x, int maxiter, REAL tol)
{
  Solver::params sprm;
  REAL error;
  int iters;

  sprm.maxiter = maxiter;
  sprm.tol = tol;
  Solver solver(std::get<0>(*matrix), sprm);

  return std::tie(iters, error) = solver(*matrix, *precond, b, x);
}

/* returns a time step task that can be preceded or followed by any
 * externally added task */
tf::Task time_step_task (int rank, REAL time, REAL step, tf::Taskflow& taskflow)
{
  return taskflow.emplace([rank, time, step]()
  {
    using namespace compute;
    uint64_t count;
    ga_counters->get (rank, cn_elements, cn_elements+1, 0, 1, &count);
    std::vector<uint64_t> data(count*el_last);
    ga_elements->get(rank, 0, count, 0, el_last, &data[0]);

    uint64_t *matnum = &data[count*el_matnum],
             *type = &data[count*el_type];
    uint64_t *nd_ri[8][2] = {{&data[count*el_nd0_rnk], &data[count*el_nd0_idx]},
                             {&data[count*el_nd1_rnk], &data[count*el_nd1_idx]},
                             {&data[count*el_nd2_rnk], &data[count*el_nd2_idx]},
                             {&data[count*el_nd3_rnk], &data[count*el_nd3_idx]},
                             {&data[count*el_nd4_rnk], &data[count*el_nd4_idx]},
                             {&data[count*el_nd5_rnk], &data[count*el_nd5_idx]},
                             {&data[count*el_nd6_rnk], &data[count*el_nd6_idx]},
                             {&data[count*el_nd7_rnk], &data[count*el_nd7_idx]}};

    uint64_t nc = 0;
    std::set<int> ranks;
    for (uint64_t i = 0; i < count; i ++)
    {
      nc += type[i];
      for (int j = 0; j < type[i]; j ++)
      {
        ranks.insert(nd_ri[j][0][i]);
      }
    }

    std::vector<REAL> X(nc), Y(nc), Z(nc), dx(nc), dy(nc), dz(nc);
    std::vector<uint64_t> index(nc);
    {
      std::map<int,uint64_t> cnt0;
      std::map<int,std::vector<REAL>> dat0;
      std::map<int,std::vector<uint64_t>> idx0;

      for (auto r : ranks)
      {
        uint64_t count;
        ga_counters->get (r, cn_nodes, cn_nodes+1, 0, 1, &count);
        cnt0[r] = count;
        dat0[r].resize(count*nd_last);
        idx0[r].resize(count);
        ga_nodes->get (r, 0, count, 0, nd_last, &dat0[r][0]);
        ga_nodeindex->get (r, 0, count, 0, 1, &idx0[r][0]);
      }

      for (uint64_t i = 0, k = 0; i < count; i ++)
      {
        for (int j = 0; j < type[i]; j ++, k ++)
        {
          int r = nd_ri[j][0][i];
          uint64_t cnt = cnt0[r];
          uint64_t idx = nd_ri[j][1][i];
          X[k] = dat0[r][cnt*nd_X + idx];
          Y[k] = dat0[r][cnt*nd_Y + idx];
          Z[k] = dat0[r][cnt*nd_Z + idx];
          dx[k] = dat0[r][cnt*nd_dx + idx];
          dy[k] = dat0[r][cnt*nd_dy + idx];
          dz[k] = dat0[r][cnt*nd_dz + idx];
          index[k] = idx0[r][idx];
        }
      }
    }

    Matrix A;
    uint64_t n = std::get<0>(A); /* system size */
    std::vector<uint64_t> &ptr = std::get<1>(A); /* row pointers */
    std::vector<uint64_t> &col = std::get<2>(A); /* column indices */
    std::vector<REAL> &val = std::get<3>(A); /* matrix values */
    std::vector<REAL> rhs; /* right hand side */
    std::vector<REAL> x; /* unkonwn solution */

    Precond *precond = update_precond (&A);

    auto [iters, error] = solve_system (&A, precond, rhs, x, 100, 1E-10);

    delete precond;
  });
}
