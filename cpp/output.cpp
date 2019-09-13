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
#include "real.h"
#include "err.h"
#include "solfec.hpp"

static uint64_t output_frame = 0;

/* append an XMF file */
static void append_xmf_file (const char *xmf_path, std::string mode,
  uint64_t elements, uint64_t nodes, uint64_t topo_size, const char *label,
  const char *h5file, std::string ent)
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
    fprintf (xmf_file, "%s:/%" PRIu64 "/TOPO\n", h5file, output_frame);
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

  if (ent == "NUMBER")
  {
    if (mode == "MESH")
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
  else if (ent == "COLOR")
  {
    if (mode == "MESH")
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
  else if (ent == "DISPL")
  {
    fprintf (xmf_file, "<Attribute Name=\"DISPL\" Center=\"Node\" AttributeType=\"Vector\">\n");
    fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 " 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
    fprintf (xmf_file, "%s:/%" PRIu64 "/DISPL\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }
  else if (ent == "LINVEL")
  {
    fprintf (xmf_file, "<Attribute Name=\"LINVEL\" Center=\"Node\" AttributeType=\"Vector\">\n");
    fprintf (xmf_file, "<DataStructure Dimensions=\"%" PRIu64 " 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
    fprintf (xmf_file, "%s:/%" PRIu64 "/LINVEL\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Attribute>\n");
  }
  else if (ent == "STRESS")
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

/* output hdf5 dataset from meshes */
static std::tuple<uint64_t, uint64_t, uint64_t> h5_mesh_dataset (std::vector<uint64_t>& subset, std::string ent, hid_t h5_step)
{
  return  {0, 0, 0};
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
        if (solfec::meshes.count(bodnum)) subset.push_back(bodnum);
      }

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

        std::set<std::string> ents = {"NUMBER", "COLOR", "DISPL", "LINVEL", "STRESS"};

        for (auto ent : iout->entities)
        {
          if (ents.count(ent))
          {
            auto [elements, nodes, topo_size] = h5_mesh_dataset (subset, ent, h5_step); /* append h5 file */
            append_xmf_file (xmf_path.str().c_str(), "MESH", elements, nodes, topo_size, label, h5file.c_str(), ent); /* append xmf file */
          }
        }

        H5Gclose (h5_step);
        H5Fclose (h5_file);
      }
    }
  }

  /* TODO: other modes */

  output_frame ++;
}
