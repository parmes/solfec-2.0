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

#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <inttypes.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include "real.h"
#include "err.h"
#include "solfec.hpp"
#include "compute.hpp"

/* rank 0 --- */
uint64_t output_frame = 0;

std::map<std::pair<uint64_t,std::vector<uint64_t>>,uint64_t> topo_map;
  /* maps hashed ordered bonum subset vectors to their initial output frames;
     pair of (hash(vector), vector itself) is used as map key for faster key
     comparisons in case of uncolliding hash values */
/* --- rank 0 */

/* append an XMF file */
static void append_xmf_file (const char *xmf_path, std::string mode,
  uint64_t elements, uint64_t nodes, uint64_t topo_size, uint64_t topo_frame,
  const char *label, const char *h5file, std::set<std::string> &ents)
{
  FILE *xmf_file;

  if (solfec::simulation_time == 0.0)
  {
    ASSERT(xmf_file = fopen (xmf_path, "w"), "XMF markup file open failed");

    fprintf (xmf_file, "<Xdmf>\n");
    fprintf (xmf_file, "<Domain>\n");
    fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
  }
  else
  {
    ASSERT (xmf_file = fopen (xmf_path, "r+"), "XMF markup file open failed");
    fseek (xmf_file, -1, SEEK_END);
    int ln = 0;
    for(;;) /* move back three lines */
    {
      if (fgetc(xmf_file) == '\n') ln ++;
      if (ln == 4) break;
      else fseek (xmf_file, -2, SEEK_CUR);
    }
    int pos = ftell(xmf_file);
    char *mem = (char*)malloc(pos);
    ERRMEM (mem);
    fseek (xmf_file, 0, SEEK_SET);
    ASSERT (fread (mem, sizeof(char), pos, xmf_file) == pos, "XMF markup file read failed"); /* read until the last three lines */
    fclose (xmf_file);
    ASSERT (xmf_file = fopen (xmf_path, "w"), "XMF markup file open failed");
    fwrite (mem, sizeof(char), pos, xmf_file); /* effectively truncate the last three lines */
    free (mem);
  }

  ASSERT (xmf_file, "XDMF file open error");

  fprintf (xmf_file, "<Grid Name=\"%s\" Type=\"Uniform\">\n", label);
  fprintf (xmf_file, "<Time Type=\"Single\" Value=\"%f\" />\n", (float)solfec::simulation_time);

  if (mode == "MESH")
  {
    fprintf (xmf_file, "<Topology Type=\"Mixed\" NumberOfElements=\"%d\">\n", elements);
    fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", topo_size);
    fprintf (xmf_file, "%s:/%" PRIu64 "/TOPO\n", h5file, topo_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Topology>\n");
  }
  else
  {
    fprintf (xmf_file, "<Topology Type=\"Polyvertex\" NumberOfElements=\"%d\" NodesPerElement=\"%d\">\n", nodes, 1);
    fprintf (xmf_file, "</Topology>\n");
  }

  fprintf (xmf_file, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
  fprintf (xmf_file, "%s:/%" PRIu64 "/GEOM\n", h5file, output_frame);
  fprintf (xmf_file, "</DataStructure>\n");
  fprintf (xmf_file, "</Geometry>\n");

  std::set<std::string> modes0 = {"MESH", "SURF", "EL"};
  std::set<std::string> modes1 = {"SURF", "EL"};

  if (ents.count("NUMBER") && modes0.count(mode))
  {
    if (mode == "MESH" || mode == "SURF")
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Cell\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 "\" NumberType=\"Int\" Format=\"HDF\">\n", elements);
    }
    else
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 "\" NumberType=\"Int\" Format=\"HDF\">\n", nodes);
    }
    fprintf (xmf_file, "%s:/%" PRIu64 "/NUMBER\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }

  if (ents.count("COLOR") && modes1.count(mode))
  {
    if (mode == "SURF")
    {
      fprintf (xmf_file, "<Attribute Name=\"COLOR\" Center=\"Cell\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 "\" NumberType=\"Int\" Format=\"HDF\">\n", elements);
    }
    else
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 "\" NumberType=\"Int\" Format=\"HDF\">\n", nodes);
    }
    fprintf (xmf_file, "%s:/%" PRIu64 "/COLOR\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }

  if (ents.count("DISPL") && modes0.count(mode))
  {
    fprintf (xmf_file, "<Attribute Name=\"DISPL\" Center=\"Node\" AttributeType=\"Vector\">\n");
    fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 " 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
    fprintf (xmf_file, "%s:/%" PRIu64 "/DISPL\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }

  if (ents.count("LINVEL") && modes0.count(mode))
  {
    fprintf (xmf_file, "<Attribute Name=\"LINVEL\" Center=\"Node\" AttributeType=\"Vector\">\n");
    fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 " 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
    fprintf (xmf_file, "%s:/%" PRIu64 "/LINVEL\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }

  if (ents.count("STRESS") && modes0.count(mode))
  {
    fprintf (xmf_file, "<Attribute Name=\"STRESS\" Center=\"Node\" AttributeType=\"Vector\">\n");
    fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 " 6\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
    fprintf (xmf_file, "%s:/%" PRIu64 "/STRESS\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }

  /* TODO: other entities */ 
  
  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Domain>\n");
  fprintf (xmf_file, "</Xdmf>\n");

  fclose (xmf_file);
}

/* write hdf5 dataset from meshes; outputs elements, nodes, topo_size counters */
static std::tuple<uint64_t, uint64_t, uint64_t> h5_mesh_dataset (std::vector<uint64_t>& subset, bool topology, std::set<std::string> &ents, hid_t h5_step)
{
  using namespace compute;
  uint64_t elements = 0,
           nodes = 0,
           topo_size = 0;

  for (auto bodnum : subset)
  {
    struct mesh &mesh = solfec::meshes[bodnum];

    nodes += mesh.nodes[0].size();
    elements += mesh.nhex+mesh.nwed+mesh.npyr+mesh.ntet;
    topo_size += 9*mesh.nhex+7*mesh.nwed+6*mesh.npyr+5*mesh.ntet;
  }

  if (topology) /* TOPO */
  {
    std::vector<uint64_t> data;
    uint64_t nc = 0;

    for (auto bodnum : subset)
    {
      struct mesh &mesh = solfec::meshes[bodnum];

      auto it = mesh.elements.begin();
      for (; it != mesh.elements.end(); )
      {
        switch (*it)
        {
         case 8:
           data.push_back(9); /* hex code */
           data.push_back(nc+it[1]); /* nodes */
           data.push_back(nc+it[2]);
           data.push_back(nc+it[3]);
           data.push_back(nc+it[4]);
           data.push_back(nc+it[5]);
           data.push_back(nc+it[6]);
           data.push_back(nc+it[7]);
           data.push_back(nc+it[8]);
           it += 10;
         break;
         case 6:
           data.push_back(8); /* wed code */
           data.push_back(nc+it[1]); /* nodes */
           data.push_back(nc+it[2]);
           data.push_back(nc+it[3]);
           data.push_back(nc+it[4]);
           data.push_back(nc+it[5]);
           data.push_back(nc+it[6]);
           it += 8;
         break;
         case 5:
           data.push_back(7); /* pyr code */
           data.push_back(nc+it[1]); /* nodes */
           data.push_back(nc+it[2]);
           data.push_back(nc+it[3]);
           data.push_back(nc+it[4]);
           data.push_back(nc+it[5]);
           it += 7;
         break;
         case 4: 
           data.push_back(6); /* tet code */
           data.push_back(nc+it[1]); /* nodes */
           data.push_back(nc+it[2]);
           data.push_back(nc+it[3]);
           data.push_back(nc+it[4]);
           it += 6;
         break;
        }
      }

      nc += mesh.nodes[0].size();
    }

    hsize_t length = data.size();
    ASSERT (H5LTmake_dataset (h5_step, "TOPO", 1, &length, H5T_NATIVE_UINT64, &data[0]) >= 0, "HDF5 file write error");
  }

  bool displ = (bool)ents.count("DISPL");
  std::vector<REAL> displ_data;
  bool linvel = (bool)ents.count("LINVEL");
  std::vector<REAL> linvel_data;

  /* GEOM */
  {
    std::vector<REAL> data;
    uint64_t size = 0;

    for (auto bodnum : subset)
    {
      struct mapping &mapping = mesh_mapping[bodnum];

      for (auto rng : mapping.ga_nranges)
      {
        auto nrng = rng[2]-rng[1];
        std::vector<REAL> vals(nrng*9);

        compute::ga_nodes->get(rng[0], rng[1], rng[2], nd_vx, nd_Z+1, &vals[0]);

        for (uint64_t i = 0; i < nrng; i ++)
        {
          REAL x = vals[nrng*nd_x+i],
               y = vals[nrng*nd_y+i],
               z = vals[nrng*nd_z+i];

          data.push_back (x);
          data.push_back (y);
          data.push_back (z);

          if (displ)
          {
            displ_data.push_back (x-vals[nrng*nd_X+i]);
            displ_data.push_back (y-vals[nrng*nd_Y+i]);
            displ_data.push_back (z-vals[nrng*nd_Z+i]);
          }

          if (linvel)
          {
            linvel_data.push_back (vals[nrng*nd_vx+i]);
            linvel_data.push_back (vals[nrng*nd_vy+i]);
            linvel_data.push_back (vals[nrng*nd_vz+i]);
          }
        }

        size += nrng;
      }
    }

    ASSERT (size == nodes, "Inconsistent nodes count during file output");

    hsize_t dims[2] = {size, 3};
#if REALSIZE==4
    ASSERT (H5LTmake_dataset_float (h5_step, "GEOM", 2, dims, &data[0]) >= 0, "HDF5 file write error");
#else
    ASSERT (H5LTmake_dataset_double (h5_step, "GEOM", 2, dims, &data[0]) >= 0, "HDF5 file write error");
#endif
  }

  if (ents.count("NUMBER"))
  {
    std::vector<uint64_t> data;

    for (auto bodnum : subset)
    {
      struct mesh &mesh = solfec::meshes[bodnum];
      auto ne = mesh.nhex+mesh.nwed+mesh.npyr+mesh.ntet;
      std::vector<uint64_t> chunk(ne,bodnum);
      data.insert(data.end(), chunk.begin(), chunk.end());
    }

    hsize_t length = data.size();
    ASSERT (H5LTmake_dataset (h5_step, "NUMBER", 1, &length, H5T_NATIVE_UINT64, &data[0]) >= 0, "HDF5 file write error");
  }

  if (ents.count("DISPL"))
  {
    hsize_t dims[2] = {nodes, 3};
#if REALSIZE==4
    ASSERT (H5LTmake_dataset_float (h5_step, "DISPL", 2, dims, &displ_data[0]) >= 0, "HDF5 file write error");
#else
    ASSERT (H5LTmake_dataset_double (h5_step, "DISPL", 2, dims, &displ_data[0]) >= 0, "HDF5 file write error");
#endif
  }

  if (ents.count("LINVEL"))
  {
    hsize_t dims[2] = {nodes, 3};
#if REALSIZE==4
    ASSERT (H5LTmake_dataset_float (h5_step, "LINVEL", 2, dims, &linvel_data[0]) >= 0, "HDF5 file write error");
#else
    ASSERT (H5LTmake_dataset_double (h5_step, "LINVEL", 2, dims, &linvel_data[0]) >= 0, "HDF5 file write error");
#endif
  }

  if (ents.count("STRESS"))
  {
    std::vector<REAL> data(nodes*6, 0.); /* FIXME: TODO (add extra nodal valus) */
    hsize_t dims[2] = {nodes, 6};
#if REALSIZE==4
    ASSERT (H5LTmake_dataset_float (h5_step, "STRESS", 2, dims, &data[0]) >= 0, "HDF5 file write error");
#else
    ASSERT (H5LTmake_dataset_double (h5_step, "STRESS", 2, dims, &data[0]) >= 0, "HDF5 file write error");
#endif
  }

  return  {elements, nodes, topo_size};
}

/* https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector */
uint64_t hashfunc(std::vector<uint64_t> const& vec)
{
  uint64_t seed = vec.size();
  for(auto& i : vec) {
    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

/* output current results */
void output_current_results()
{
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank > 0) return; /* XXX: start with a simple rank 0 write */

  std::ostringstream h5_path, h5_text, xmf_path;
  FILE *xmf_file;
  hid_t h5_file;
  hid_t h5_step;

  for (auto iout = solfec::outputs.begin(); iout != solfec::outputs.end(); iout ++)
  {
    if (iout->modes.count("MESH"))
    {
      std::vector<uint64_t> subset;

      for (auto bodnum : iout->subset)
      {
        if (solfec::meshes.count(bodnum)) subset.push_back(bodnum); /* ordered */
      }

      if (subset.empty()) /* use all meshes */
      {
        for (const auto& pair : solfec::meshes)
        {
          subset.push_back(pair.first);
        }
      }

      uint64_t hash = hashfunc(subset);

      auto key = std::make_pair(hash,subset);

      uint64_t topo_frame;

      if (!topo_map.count(key)) /* map current output frame for this subset: this is where mesh topology dataset will be stored */
      {
        topo_map[key] = output_frame;

        topo_frame = output_frame;
      }
      else topo_frame = topo_map[key]; /* use previous output frame where mesh topology was stored */
      
      if (subset.size() > 0)
      {
        auto nout = iout-solfec::outputs.begin();

        h5_path.str("");
        h5_path.clear();
        h5_path << solfec::outname << nout << ".h5";

        if (solfec::simulation_time == 0.0)
        {
          ASSERT ((h5_file = H5Fcreate(h5_path.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file open error");
        }
        else
        {
          ASSERT((h5_file = H5Fopen(h5_path.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT)) >= 0, "HDF5 file open error");
        }

        h5_text.str("");
        h5_text.clear();
        h5_text << output_frame;
        ASSERT ((h5_step = H5Gcreate (h5_file, h5_text.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
        double time = solfec::simulation_time;
        ASSERT (H5LTset_attribute_double (h5_step, ".", "TIME", &time, 1) >= 0, "HDF5 file write error"); 

        xmf_path.str("");
        xmf_path.clear();
        xmf_path << solfec::outname << nout << ".xmf";
        const char *label = "SOLFEC-2.0 meshes";
        std::string h5file = h5_path.str().substr(h5_path.str().find_last_of('/')+1);

        auto [elements, nodes, topo_size] = h5_mesh_dataset (subset, output_frame == topo_frame, iout->entities, h5_step); /* append h5 file */

        append_xmf_file (xmf_path.str().c_str(), "MESH", elements, nodes, topo_size, topo_frame, label, h5file.c_str(), iout->entities); /* append xmf file */

        H5Gclose (h5_step);
        H5Fclose (h5_file);
      }
    }
  }

  /* TODO: other modes */

  output_frame ++;
}
