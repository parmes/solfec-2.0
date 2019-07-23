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

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "real.h"
#include "err.h"
#include "alg.h"
#include "fmt.h"
#include "solfec.h"
#include "util_ispc.h"

#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsSize_t
#define PyString_FromString PyUnicode_FromString
#define PyString_Check PyUnicode_Check
#define PyString_AsString PyUnicode_AsUTF8
#define Py_RETURN_uint64_t(arg) return PyLong_FromUnsignedLongLong(arg)
#else
#define Py_RETURN_uint64_t(arg) return Py_BuildValue("K", arg)
#endif


#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#endif

#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE return Py_INCREF(Py_True), Py_True
#endif

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

/* string buffer length */
#define BUFLEN 512

/* minimal type initialization */
#define TYPEINIT(typedesc, type, name, flags, dealloc, new, methods, members, getset)\
memset (&(typedesc), 0, sizeof (PyTypeObject));\
(typedesc).tp_basicsize = sizeof (type);\
(typedesc).tp_name = name;\
(typedesc).tp_flags = flags;\
(typedesc).tp_dealloc = (destructor)dealloc;\
(typedesc).tp_new = new;\
(typedesc).tp_methods = methods;\
(typedesc).tp_members = members;\
(typedesc).tp_getset = getset

/* bool test */
static int is_bool (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!PyBool_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be Boolean (True/False)", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* string test */
static int is_string (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!PyString_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a string", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a list, of length >= len divisible by div, or a string */
static int is_list_or_string (PyObject *obj, const char *var, int div, int len)
{
  if (obj)
  {
    if (!(PyList_Check (obj) || PyString_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list or a string object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have N * %d elements, where N >= %d", var, div, len / div);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* list test */
static int is_list (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyList_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d items", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* tuple test */
static int is_tuple (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!PyTuple_Check (obj))
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyTuple_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d elements", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* positive test */
static int is_positive (double num, const char *var)
{
  if (num <= 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be positive", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* non-negative test */
static int is_non_negative (double num, const char *var)
{
  if (num < 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be non-negative", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* in map test */
template <typename map_type>
static int is_in_map (double objnum, const char *var, map_type objmap)
{
  typename map_type::iterator it = objmap.find(objnum);

  if (it == objmap.end())
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' object with number %zu does not exist", var, objnum);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* test whether an object is a list (details as above) or a number */
static int is_list_or_number (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!(PyList_Check (obj) || PyNumber_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list or a number object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (len > 0 && PyList_Size (obj) != len)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have %d items", var, len);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* tuple or list of tuples check */
static int is_tuple_or_list_of_tuples (PyObject *obj, const char *var, size_t *tuple_lengths, size_t n_tuple_lengths)
{
  size_t i, j, k, n;

  if (obj)
  {
    if (PyTuple_Check (obj))
    {
      j = PyTuple_Size (obj);

      for (i = 0; i < n_tuple_lengths; i ++)
      {
        if (j == tuple_lengths[i]) break;
      }

      if (i == n_tuple_lengths)
      {
	char buf [BUFLEN];
	snprintf (buf, BUFLEN, "tuple '%s' length is invalid", var);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
    else if (PyList_Check (obj))
    {
      n = PyList_Size (obj);

      for (i = 0; i < n; i ++)
      {
	PyObject *item = PyList_GetItem (obj, i);

	if (!PyTuple_Check (item))
	{
	  char buf [BUFLEN];
	  sprintf (buf, "'%s' must be a list of tuples: item %d is not a tuple", var, i);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return 0;
	}

	j = PyTuple_Size (item);

	for (k = 0; k < n_tuple_lengths; k ++)
	{
	  if (j == tuple_lengths[k]) break;
	}

	if (k == n_tuple_lengths)
	{
	  char buf [BUFLEN];
	  sprintf (buf, "'%s' list item %d (a tuple) has invalid length", var, i);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return 0;
	}
      }

      return n;
    }
    else
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple or a list of tuples", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* transform mesh nodes */
static void dotransform (REAL *transform, size_t type, std::array<std::vector<REAL>,3> &nodes)
{
  switch (type)
  {
  case 3: /* translate */
    ispc::translate (transform, &nodes[0][0], &nodes[1][0], &nodes[2][0], nodes[0].size());
  break;
  case 7: /* axis rotate */
  {
    REAL rotation[9];
    ROTATION_MATRIX (transform+3, transform[7], rotation);
    ispc::axis_rotate (transform, rotation, &nodes[0][0], &nodes[1][0], &nodes[2][0], nodes[0].size());
  }
  break;
  case 9: /* matrix transform */
    ispc::matrix_transform (transform, &nodes[0][0], &nodes[1][0], &nodes[2][0], nodes[0].size());
  break;
  case 4: /* scale */
    ispc::scale (transform, transform[3], &nodes[0][0], &nodes[1][0], &nodes[2][0], nodes[0].size());
  break;
  }
}

/* define keywords */
#define KEYWORDS(...) const char *kwl [] = {__VA_ARGS__, NULL}

/* parse arguments with keywords */
#define PARSEKEYS(fmt, ...) if (!PyArg_ParseTupleAndKeywords (args, kwds, fmt, (char**)kwl, __VA_ARGS__)) return NULL

/* parse arguments without keywords */
#define PARSE(fmt, ...) if (!PyArg_ParseTuple (args, fmt, __VA_ARGS__)) return NULL

/* object types assertion */
#define TYPETEST(test) if(!(test)) return NULL

/* string argument if block comparison */
#define IFIS(obj, val) if (strcmp (PyString_AsString (obj), val) == 0)
#define ELIF(obj, val) else if (strcmp (PyString_AsString (obj), val) == 0)
#define ELSE else

/* test if a string has a specific ending */
static int endswith (const char *string, const char *ending)
{
  if (strlen (string) >= strlen (ending) &&
      strcmp (string+strlen(string)-strlen(ending), ending) == 0) return 1;
  else return 0;
}

/* command line arguments */
static PyObject* ARGV (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nonsolfec");
  PyObject *nonsolfec, *list;
  using namespace solfec; /* argc, argv */

  nonsolfec = Py_True;

  PARSEKEYS ("|O", &nonsolfec);

  TYPETEST (is_bool (nonsolfec, kwl [0]));

  if (!(list = PyList_New (0))) return NULL;

  if (nonsolfec == Py_True)
  {
    for (int i = 0; i < argc; i ++)
    {
      if (endswith (argv[i], "solfec4")) continue;
      else if (endswith (argv[i], "solfec8")) continue;
      else if (strcmp (argv[i], "-ntasks") == 0)
      {
	i ++;
	continue;
      }
      else if (strlen (argv[i]) > 3 && strcmp(argv[i]+strlen(argv[i])-3, ".py") == 0) continue;
      else PyList_Append (list, PyString_FromString (argv[i]));
    }
  }
  else
  {
    for (int i = 0; i < argc; i ++)
    {
      if (endswith (argv[i], "solfec4")) continue;
      else if (endswith (argv[i], "solfec8")) continue;
      else PyList_Append (list, PyString_FromString (argv[i]));
    }
  }

  return list;
}

/* reset simulation */
static PyObject* RESET (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("outname");
  PyObject *outname;

  outname = NULL;

  PARSEKEYS ("|O", &outname);

  TYPETEST (is_string (outname, kwl [0]));

  if (outname)
  {
    solfec::outname.assign(PyString_AsString(outname));
  }

  solfec::splines.clear();
  solfec::splines_count = 0;
  solfec::materials.clear();
  solfec::materials_count = 0;
  solfec::meshes.clear();
  solfec::meshes_count = 0;
  solfec::restrains.clear();
  solfec::restrains_count = 0;
  solfec::prescribes.clear();
  solfec::prescribes_count = 0;
  solfec::velocities.clear();
  solfec::frictions.clear();
  solfec::gravity.clear();
  solfec::histories.clear();
  solfec::outputs.clear();
  solfec::notrun = true;

  Py_RETURN_NONE;
}

/* Create linear spline */
static PyObject* SPLINE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("points", "cache");
  std::array<REAL,2> temp;
  PyObject *points;
  int cache, i, n;

  cache = 0;

  PARSEKEYS ("O|i", &points, &cache);

  TYPETEST (is_list_or_string (points, kwl [0], 1, 0));

  if (cache < 0)
  {
    PyErr_SetString (PyExc_ValueError, "Negative cache size");
    return NULL;
  }

  struct spline &spline = solfec::splines[solfec::splines_count ++];

  if (PyString_Check (points))
  {
    spline_from_file (PyString_AsString(points), cache, &spline);
  }
  else if (PyList_Check (points))
  {
    if (PyList_Check (PyList_GetItem (points, 0)))
    {
      n = PyList_Size (points);

      for (i = 0; i < n; i ++)
      {
	PyObject *pv = PyList_GetItem (points, i);

	TYPETEST (is_list (pv, "[x,y]", 2));

	temp[0] = PyFloat_AsDouble (PyList_GetItem (pv, 0));
	temp[1] = PyFloat_AsDouble (PyList_GetItem (pv, 1));

	spline.points.push_back (temp);
      }
    }
    else if (PyTuple_Check (PyList_GetItem (points, 0)))
    {
      n = PyList_Size (points);

      for (i = 0; i < n; i ++)
      {
	PyObject *pv = PyList_GetItem (points, i);

	TYPETEST (is_tuple (pv, "(x,y)", 2));

	temp[0] = PyFloat_AsDouble (PyTuple_GetItem (pv, 0));
	temp[1] = PyFloat_AsDouble (PyTuple_GetItem (pv, 1));

	spline.points.push_back (temp);
      }
    }
    else
    {
      TYPETEST (is_list_or_string (points, kwl [0], 2, 4));

      n = PyList_Size (points) / 2;

      for (i = 0; i < n; i ++)
      {
	temp[0] = PyFloat_AsDouble (PyList_GetItem (points, 2*i));
	temp[1] = PyFloat_AsDouble (PyList_GetItem (points, 2*i + 1));

	spline.points.push_back (temp);
      }
    }
  }
  else
  {
    PyErr_SetString (PyExc_TypeError, "Invalid points format");
    return NULL;
  }

  Py_RETURN_uint64_t (solfec::splines_count-1);
}

/* print linear spline */
static PyObject* print_SPLINE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("splnum");
  size_t splnum;

  if (solfec::notrun)
  {
    PARSEKEYS ("K", &splnum);

    TYPETEST (is_in_map (splnum, "SPLINE", solfec::splines));

    struct spline &spline = solfec::splines[splnum];

    std::cout << "SPLINE_" << splnum << "_points = ";
    std::vector<std::array<REAL,2>>::iterator it = spline.points.begin();
    if (it == spline.points.end()) std::cout << std::endl; else std::cout << "[";
    for (; it != spline.points.end(); it ++)
    {
      std::cout << "(" << (*it)[0] << "," << (*it)[1];
      if (it+1 != spline.points.end()) std::cout << "),";
      else std::cout << ")]" << std::endl;
    }
    std::cout << "SPLINE_" << splnum << "_cache = " << spline.cache << std::endl;
    std::cout << "SPLINE_" << splnum << "_path = " << spline.path << std::endl;
  }

  Py_RETURN_NONE;
}

/* Create material */
static PyObject* MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("density", "young", "poisson", "viscosity");
  double density, young, poisson, viscosity;

  PARSEKEYS ("dddd", &density, &young, &poisson, &viscosity);

  TYPETEST (is_positive (density, kwl[0]) && is_positive (young, kwl[1]) &&
            is_positive (poisson, kwl[2]) && is_non_negative (viscosity, kwl[3]));

  struct material &material = solfec::materials[solfec::materials_count++];

  material.density = (REAL)density;
  material.young = (REAL)young;
  material.poisson = (REAL)poisson;
  material.viscosity = (REAL)viscosity;

  Py_RETURN_uint64_t (solfec::materials_count-1);
}

/* print material */
static PyObject* print_MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("matnum");
  size_t matnum;

  if (solfec::notrun)
  {
    PARSEKEYS ("K", &matnum);

    TYPETEST (is_in_map (matnum, "MATERIAL", solfec::materials));

    struct material &material = solfec::materials[matnum];

    std::cout << "MATERIAL_" << matnum << "_density = " << material.density<< std::endl;
    std::cout << "MATERIAL_" << matnum << "_young = " << material.young << std::endl;
    std::cout << "MATERIAL_" << matnum << "_poisson = " << material.poisson << std::endl;
    std::cout << "MATERIAL_" << matnum << "_viscosity = " << material.viscosity << std::endl;
  }

  Py_RETURN_NONE;
}

/* Create meshed body */
static PyObject* MESH (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nodes", "elements", "matnum", "colors", "transform");
  size_t tuple_lengths[4] = {3, 7, 9, 4}, i, j, k, l, n, m, e;
  PyObject *nodes, *elements, *colors, *transform;
  size_t matnum;

  transform = NULL;

  PARSEKEYS ("OOKO|O", &nodes, &elements, &matnum, &colors, &transform);

  TYPETEST (is_list (nodes, kwl[0], 0) && is_list (elements, kwl[1], 0) &&
            is_list_or_number (colors, kwl[3], 0) &&
            is_tuple_or_list_of_tuples (transform, kwl[4], tuple_lengths, 4));

  struct mesh &mesh = solfec::meshes[solfec::meshes_count++];

  mesh.matnum = matnum;

  /* read nodes */
  if (PyList_Size(nodes) % 3)
  {
    PyErr_SetString (PyExc_ValueError, "Nodes list length must be a multiple of 3");
    return NULL;
  }

  n = PyList_Size (nodes) / 3; /* nodes count */

  for (i = 0; i < n; i ++) 
  {
    std::array<REAL,3> temp = {(REAL)PyFloat_AsDouble (PyList_GetItem (nodes, 3*i+0)),
			       (REAL)PyFloat_AsDouble (PyList_GetItem (nodes, 3*i+1)),
			       (REAL)PyFloat_AsDouble (PyList_GetItem (nodes, 3*i+2))};

    mesh.nodes[0].push_back (temp[0]);
    mesh.nodes[1].push_back (temp[1]);
    mesh.nodes[2].push_back (temp[2]);
  }

  /* read elements */
  l = PyList_Size (elements);
  for (e = i = 0; i < l; e ++)
  {
    k = PyInt_AsLong (PyList_GetItem (elements, i));

    switch(k)
    {
    case 4:
      mesh.ntet++;
    break;
    case 5:
      mesh.npyr++;
    break;
    case 6:
      mesh.nwed++;
    break;
    case 8:
      mesh.nhex++;
    break;
    default:
      PyErr_SetString (PyExc_ValueError, "An element must have 4, 5, 6, or 8 nodes");
      return NULL;
    break;
    }

    mesh.elements.push_back (k);

    for (j = i+1; j <= i+k; j ++) /* nodes */
    {
      m = PyInt_AsLong (PyList_GetItem (elements, j));

      if (m < 0 || m >= n)
      {
	char buf [BUFLEN];
	sprintf (buf, "Node %d in element %d is outside of range [0, %d]", j-i, e, n-1);
	PyErr_SetString (PyExc_ValueError, buf);
	return NULL;
      }

      if (std::count(mesh.elements.end()-(j-i-1),mesh.elements.end(),m)>0)
      {
	char buf [BUFLEN];
	sprintf (buf, "Node %d repeats in element %d", m, e);
	PyErr_SetString (PyExc_ValueError, buf);
	return NULL;
      }

      mesh.elements.push_back (m);
    }

    m = PyInt_AsLong (PyList_GetItem (elements, j)); /* material */

    mesh.elements.push_back (m);

    i = ++j;

    if (i > l) /* incomplete */
    {
      PyErr_SetString (PyExc_ValueError, "The last element definition is incomplete");
      return NULL;
    }
  }

  if (PyList_Check (colors))
  {
    /* global color */
    mesh.gcolor = PyInt_AsLong (PyList_GetItem (colors, 0));

    TYPETEST (is_positive (mesh.gcolor, "[global] color"));

    /* read face colors */
    l = PyList_Size (colors);
    for (e = 0, i = 1; i < l; e ++)
    {
      k = PyInt_AsLong (PyList_GetItem (colors, i));

      if (!(k == 3 || k == 4))
      {
	PyErr_SetString (PyExc_ValueError, "A face must have 3 or 4 nodes");
	return NULL;
      }

      mesh.colors.push_back (k);

      for (j = i+1; j <= i+k; j ++) /* nodes */
      {
	m = PyInt_AsLong (PyList_GetItem (colors, j));

	if (m < 0 || m >= n)
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Node %d in face %d is outside of range [0, %d]", j-i, e, n-1);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}

	if (std::count(mesh.colors.end()-(j-i-1),mesh.colors.end(),m)>0)
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Node %d repeats in face %d", m, e);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}

	mesh.colors.push_back (m);
      }

      m = PyInt_AsLong (PyList_GetItem (colors, j)); /* color */

      TYPETEST (is_positive (m, "[face] color"));

      mesh.colors.push_back (m);

      i = ++j;

      if (i > l) /* incomplete */
      {
	PyErr_SetString (PyExc_ValueError, "The last color definition is incomplete");
	return NULL;
      }
    }

    mesh.nfaces = e;
  }
  else
  {
    mesh.gcolor = PyInt_AsLong (colors);

    TYPETEST (is_positive (mesh.gcolor, "[global] color"));
  }

  /* tranform nodes */
  if (transform)
  {
    if (PyTuple_Check (transform))
    {
      REAL trans[9];
      size_t type = PyTuple_Size (transform);
      for (i = 0; i < type; i ++)
	trans[i] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (transform, i)),
      dotransform (trans, type, mesh.nodes);
    }
    else
    {
      n = PyList_Size (transform);
      for (i = 0; i < n; i ++)
      {
	REAL trans[9];
	PyObject *tuple = PyList_GetItem (transform, i);
	size_t type = PyTuple_Size (tuple);
	for (j = 0; j < type; j ++)
	  trans[j] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (tuple, j)),
	dotransform (trans, type, mesh.nodes);
      }
    }
  }

  Py_RETURN_uint64_t (solfec::meshes_count-1);
}

/* print meshed body */
static PyObject* print_MESH (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("bodnum");
  size_t bodnum;

  if (solfec::notrun)
  {
    PARSEKEYS ("K", &bodnum);

    TYPETEST (is_in_map (bodnum, "MESH", solfec::meshes));

    struct mesh &mesh = solfec::meshes[bodnum];

    {
      std::cout << "MESH_" << bodnum << "_nodes = [";
      std::vector<REAL>::iterator it0 = mesh.nodes[0].begin();
      std::vector<REAL>::iterator it1 = mesh.nodes[1].begin();
      std::vector<REAL>::iterator it2 = mesh.nodes[2].begin();
      for (; it0 != mesh.nodes[0].end()-1; it0 ++, it1++, it2++)
      {
	 std::cout << (*it0) << "," << (*it1) << "," << (*it2) << "," << std::endl;
      }
      std::cout << (*it0) << "," << (*it1) << "," << (*it2) << "]" << std::endl;
    }
    {
      std::cout << "MESH_" << bodnum << "_elements = [";
      size_t ne = mesh.elements.size();
      std::vector<size_t>::iterator it = mesh.elements.begin();
      for (; it != mesh.elements.end(); )
      {
	switch (*it)
	{
	 case 8: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4] <<
	                   "," << it[5] << "," << it[6] << "," << it[7] << "," << it[8] << "," << it[9];
	 it += 10;
	 break;
	 case 6: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4] <<
	                   "," << it[5] << "," << it[6] << "," << it[7];
	 it += 8;
	 break;
	 case 5: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4] <<
	                   "," << it[5] << "," << it[6];
	 it += 7;
	 break;
	 case 4: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4] <<
	                   "," << it[5];
	 it += 6;
	 break;
	}

	if (it != mesh.elements.end()) std::cout << "," << std::endl;
	else std::cout << "]" << std::endl;
      }
    }
    {
      std::cout << "MESH_" << bodnum << "_colors = [" << mesh.gcolor << ",";
      std::vector<size_t>::iterator it = mesh.colors.begin();
      for (; it != mesh.colors.end(); )
      {
	switch (*it)
	{
	 case 4: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4] <<
	                   "," << it[5];
	 it += 6;
	 break;
	 case 3: std::cout << it[0] << "," << it[1] << "," << it[2] << "," << it[3] << "," << it[4];
	 it += 5;
	 break;
	}

	if (it != mesh.colors.end()) std::cout << "," << std::endl;
	else std::cout << "]" << std::endl;
      }
    }
  }

  Py_RETURN_NONE;
}
 
/* Restrain motion */
static PyObject* RESTRAIN (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("bodnum", "point", "color", "direction");
  size_t tuple_lengths[1] = {3}, i, j, n;
  PyObject *point, *direction;

  struct restrain &restrain = solfec::restrains[solfec::restrains_count++];

  restrain.color = 0;
  direction = NULL;
  point = NULL;

  PARSEKEYS ("K|OKO", &restrain.bodnum, &point, &restrain.color, &direction);

  TYPETEST (is_in_map (restrain.bodnum, "MESH", solfec::meshes) &&
            is_tuple_or_list_of_tuples (point, kwl[1], tuple_lengths, 1) &&
	    is_tuple (direction, kwl[3], 3));

  if (point)
  {
    std::array<REAL,3> temp;
    if (PyTuple_Check(point))
    {
      for (i = 0; i < 3; i ++)
	temp[i] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (point, i));
      restrain.points.push_back(temp);
    }
    else
    {
      n = PyList_Size (point);
      for (i = 0; i < n; i ++)
      {
	PyObject *tuple = PyList_GetItem (point, i);
	for (j = 0; j < 3; j ++)
	  temp[j] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (tuple, j));
        restrain.points.push_back(temp);
      }
    }
  }

  if (direction)
  {
    for (i = 0; i < 3; i ++)
    {
      restrain.direction[i] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (direction, i));
    }
  }

  solfec::meshes[restrain.bodnum].restrains.insert(solfec::restrains_count-1);

  Py_RETURN_uint64_t (solfec::restrains_count-1);
}

/* print restrain */
static PyObject* print_RESTRAIN (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("resnum");
  size_t resnum;

  if (solfec::notrun)
  {
    PARSEKEYS ("K", &resnum);

    TYPETEST (is_in_map (resnum, "RESTRAIN", solfec::restrains));

    struct restrain &restrain = solfec::restrains[resnum];

    std::cout << "RESTRAIN_" << resnum << "_bodnum = " << restrain.bodnum << std::endl;
    if (restrain.points.size() > 0)
    {
      std::cout << "RESTRAIN_" << resnum << "_points = [";
      std::vector<std::array<REAL,3>>::iterator it = restrain.points.begin();
      for (; it != restrain.points.end(); it ++)
      {
	std::cout << "(" << (*it)[0] << "," << (*it)[1] <<  "," << (*it)[2];
	if (it+1 != restrain.points.end()) std::cout << "),";
	else std::cout << ")]" << std::endl;
      }
    }
    if (restrain.color > 0)
    {
      std::cout << "RESTRAIN_" << resnum << "_color = " << restrain.color << std::endl;
    }
    if (DOT(restrain.direction,restrain.direction)>0.)
    {
      std::cout << "RESTRAIN_" << resnum << "_direction = (" << restrain.direction[0] <<
        "," << restrain.direction[1] << "," << restrain.direction[2] << ")" << std::endl;
    }
  }

  Py_RETURN_NONE;
}

/* Prescribe moion */
static PyObject* PRESCRIBE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("bodnum", "point", "color", "linear", "angular");
  size_t tuple_lengths[1] = {3}, i, j, n;
  PyObject *point, *linear, *angular;

  struct prescribe &prescribe = solfec::prescribes[solfec::prescribes_count++];

  prescribe.color = 0;
  point = NULL;
  linear = NULL;
  angular = NULL;

  PARSEKEYS ("K|OKOO", &prescribe.bodnum, &point, &prescribe.color, &linear, &angular);

  TYPETEST (is_in_map (prescribe.bodnum, "MESH", solfec::meshes) &&
            is_tuple_or_list_of_tuples (point, kwl[1], tuple_lengths, 1));

  if (point)
  {
    std::array<REAL,3> temp;
    if (PyTuple_Check(point))
    {
      for (i = 0; i < 3; i ++)
	temp[i] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (point, i));
      prescribe.points.push_back(temp);
    }
    else
    {
      n = PyList_Size (point);
      for (i = 0; i < n; i ++)
      {
	PyObject *tuple = PyList_GetItem (point, i);
	for (j = 0; j < 3; j ++)
	  temp[j] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (tuple, j));
        prescribe.points.push_back(temp);
      }
    }
  }

  if (linear)
  {
    prescribe.linear_applied = true;

    if (PyCallable_Check(linear))
    {
      prescribe.linear_callback = linear;
      prescribe.linear_values[0] = 0.;
      prescribe.linear_values[1] = 0.;
      prescribe.linear_values[2] = 0.;
      prescribe.linear_splines[0] = -1;
      prescribe.linear_splines[1] = -1;
      prescribe.linear_splines[2] = -1;
    }
    else if (PyTuple_Check(linear))
    {
      prescribe.linear_callback = NULL;

      if (PyTuple_Size(linear) != 3) 
      {
	PyErr_SetString (PyExc_ValueError, "'linear' tuple size != 3");
	return NULL;
      }

      PyObject *l[3] =  {PyTuple_GetItem(linear, 0), PyTuple_GetItem(linear, 1), PyTuple_GetItem(linear, 2)}; 

      for (int i = 0; i < 3; i ++)
      {
#if 0
	PyTypeObject* type = l[i]->ob_type;
        const char* p = type->tp_name;
	std::cout << "linear[i] Python Type = " << p << std::endl;
#endif
	if (PyFloat_Check(l[i]))
	{
          prescribe.linear_splines[i] = -1;
          prescribe.linear_values[i] = (REAL)PyFloat_AsDouble(l[i]);
	}
	else if (PyLong_Check(l[i]))
	{
	  size_t splnum = PyInt_AsLong(l[i]);

	  if (is_in_map(splnum, "SPLINE", solfec::splines))
	  {
            prescribe.linear_splines[i] = PyInt_AsLong(l[i]);
            prescribe.linear_values[i] = 0.0;
	  }
	  else
	  {
	    PyErr_SetString (PyExc_ValueError, "'linear' tuple item SPLINE number out of range");
	    return NULL;
	  }
	}
	else
	{
	  PyErr_SetString (PyExc_ValueError, "'linear' tuple item is neither a float nor an integer");
	  return NULL;
	}
      }
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "'linear' is neither a tuple nor a callback");
      return NULL;
    }
  }
  else prescribe.linear_applied = false;

  if (angular)
  {
    prescribe.angular_applied = true;

    if (PyCallable_Check(angular))
    {
      prescribe.angular_callback = angular;
      prescribe.angular_values[0] = 0.;
      prescribe.angular_values[1] = 0.;
      prescribe.angular_values[2] = 0.;
      prescribe.angular_splines[0] = -1;
      prescribe.angular_splines[1] = -1;
      prescribe.angular_splines[2] = -1;
    }
    else if (PyTuple_Check(angular))
    {
      prescribe.angular_callback = NULL;

      if (PyTuple_Size(angular) != 3) 
      {
	PyErr_SetString (PyExc_ValueError, "'angular' tuple size != 3");
	return NULL;
      }

      PyObject *l[3] =  {PyTuple_GetItem(angular, 0), PyTuple_GetItem(angular, 1), PyTuple_GetItem(angular, 2)}; 

      for (int i = 0; i < 3; i ++)
      {
	if (PyFloat_Check(l[i]))
	{
          prescribe.angular_splines[i] = -1;
          prescribe.angular_values[i] = (REAL)PyFloat_AsDouble(l[i]);
	}
	else if (PyLong_Check(l[i]))
	{
	  size_t splnum = PyInt_AsLong(l[i]);

	  if (is_in_map(splnum, "SPLINE", solfec::splines))
	  {
            prescribe.angular_splines[i] = PyInt_AsLong(l[i]);
            prescribe.angular_values[i] = 0.0;
	  }
	  else
	  {
	    PyErr_SetString (PyExc_ValueError, "'angular' tuple item SPLINE number out of range");
	    return NULL;
	  }
	}
	else
	{
	  PyErr_SetString (PyExc_ValueError, "'angular' tuple item is neither a float nor an integer");
	  return NULL;
	}
      }
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "'angular' is neither a tuple nor a callback");
      return NULL;
    }
  }
  else prescribe.angular_applied = false;

  solfec::meshes[prescribe.bodnum].prescribes.insert(solfec::prescribes_count-1);

  Py_RETURN_uint64_t (solfec::prescribes_count-1);
}

/* print prescribe */
static PyObject* print_PRESCRIBE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("prenum");
  size_t prenum;

  if (solfec::notrun)
  {
    PARSEKEYS ("K", &prenum);

    TYPETEST (is_in_map (prenum, "PRESCRIBE", solfec::prescribes));

    struct prescribe &prescribe = solfec::prescribes[prenum];

    std::cout << "PRESCRIBE_" << prenum << "_bodnum = " << prescribe.bodnum << std::endl;
    std::cout << "PRESCRIBE_" << prenum << "_points = ";
    std::vector<std::array<REAL,3>>::iterator it = prescribe.points.begin();
    if (it == prescribe.points.end()) std::cout << std::endl; else std::cout << "[";
    for (; it != prescribe.points.end(); it ++)
    {
      std::cout << "(" << (*it)[0] << "," << (*it)[1] <<  "," << (*it)[2];
      if (it+1 != prescribe.points.end()) std::cout << "),";
      else std::cout << ")]" << std::endl;
    }
    if (prescribe.color > 0)
    {
      std::cout << "PRESCRIBE_" << prenum << "_color = " << prescribe.color << std::endl;
    }
    if (prescribe.linear_applied)
    {
      std::cout << "PRESCRIBE_" << prenum << "_linear = ";
      if (prescribe.linear_callback) std::cout << "callback" << std::endl;
      else
      {
	std::cout << "(";
	for (int i = 0; i < 3; i ++)
	{
	  if (prescribe.linear_splines[i] >= 0) std::cout << prescribe.linear_splines[i];
	  else std::cout << FMT("%e") << prescribe.linear_values[i];
	  if (i == 0 || i == 1) std::cout << ",";
	}
	std::cout << ")" << std::endl;
      }
    }
    if (prescribe.angular_applied)
    {
      std::cout << "PRESCRIBE_" << prenum << "_angular = ";
      if (prescribe.angular_callback) std::cout << "callback" << std::endl;
      else
      {
	std::cout << "(";
	for (int i = 0; i < 3; i ++)
	{
	  if (prescribe.angular_splines[i] >= 0) std::cout << prescribe.angular_splines[i];
	  else std::cout << FMT("%e") << prescribe.angular_values[i];
	  if (i == 0 || i == 1) std::cout << ",";
	}
	std::cout << ")" << std::endl;
      }
    }
  }

  Py_RETURN_NONE;
}

/* Set rigid velocity */
static PyObject* VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("bodnum", "linear", "angular");
  PyObject *linear, *angular;
  struct velocity velocity; 

  linear = angular = NULL;

  PARSEKEYS ("K|OO", &velocity.bodnum, &linear, &angular);

  TYPETEST (is_in_map (velocity.bodnum, "MESH", solfec::meshes) &&
            is_tuple (linear, kwl[1], 3) && is_tuple (angular, kwl[2], 3));

  if (linear)
  {
    velocity.linear_values[0] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (linear, 0));
    velocity.linear_values[1] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (linear, 1));
    velocity.linear_values[2] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (linear, 2));
  }
  else
  {
    velocity.linear_values[0] = 
    velocity.linear_values[1] = 
    velocity.linear_values[2] = 0.0;
  }

  if (angular)
  {
    velocity.angular_values[0] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (angular, 0));
    velocity.angular_values[1] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (angular, 1));
    velocity.angular_values[2] = (REAL)PyFloat_AsDouble (PyTuple_GetItem (angular, 2));
  }
  else
  {
    velocity.angular_values[0] = 
    velocity.angular_values[1] = 
    velocity.angular_values[2] = 0.0;
  }

  solfec::velocities.push_back (velocity);

  Py_RETURN_NONE;
}

/* print velocities */
static PyObject* print_VELOCITIES (PyObject *self, PyObject *args, PyObject *kwds)
{
  if (solfec::notrun)
  {
    for (std::vector<velocity>::iterator it = solfec::velocities.begin(); it != solfec::velocities.end(); it++)
    {
      size_t velnum = it - solfec::velocities.begin();

      std::cout << "VELOCITY_" << velnum << "_bodnum = " << (*it).bodnum << std::endl;
      std::cout << "VELOCITY_" << velnum << "_linear = (" << (*it).linear_values[0] << ","
		<< (*it).linear_values[1] << "," << (*it).linear_values[2] << ")" << std::endl;
      std::cout << "VELOCITY_" << velnum << "_angular = (" << (*it).angular_values[0] << ","
		<< (*it).angular_values[1] << "," << (*it).angular_values[2] << ")" << std::endl;
    }
  }

  Py_RETURN_NONE;
}

/* Define friction parameters */
static PyObject* FRICTION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("color1", "color2", "static", "dynamic");
  double static_friction, dynamic_friction;
  struct friction friction;

  PARSEKEYS ("KKdd", &friction.color1, &friction.color2, &static_friction, &dynamic_friction);
  friction.static_friction = (REAL)static_friction;
  friction.dynamic_friction = (REAL)dynamic_friction;

  if (friction.color1 > friction.color2) std::swap(friction.color1, friction.color2);

  TYPETEST (is_non_negative (friction.static_friction, kwl[2]) && is_non_negative (friction.dynamic_friction, kwl[3]));

  std::pair<std::set<struct friction,friction_compare>::iterator,bool> ret = solfec::frictions.insert (friction);

  if (ret.second == false) /* redefine */
  {
    solfec::frictions.erase(ret.first);
    solfec::frictions.insert (friction);
  }

  Py_RETURN_NONE;
}

/* print frictions */
static PyObject* print_FRICTIONS (PyObject *self, PyObject *args, PyObject *kwds)
{
  size_t frinum = 0;

  if (solfec::notrun)
  {
    for (std::set<friction,friction_compare>::iterator it = solfec::frictions.begin(); it != solfec::frictions.end(); it++, frinum ++)
    {
      std::cout << "FRICTION_" << frinum << "_color1 = " << (*it).color1 << std::endl;
      std::cout << "FRICTION_" << frinum << "_color2 = " << (*it).color2 << std::endl;
      std::cout << "FRICTION_" << frinum << "_static = " << (*it).static_friction << std::endl;
      std::cout << "FRICTION_" << frinum << "_dynamic = " << (*it).dynamic_friction << std::endl;
    }
  }

  Py_RETURN_NONE;
}

/* Set gravity */
static PyObject* GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("gx", "gy", "gz");
  PyObject *g[3];

  PARSEKEYS ("OOO", &g[0], &g[1], &g[2]);

  for (int i = 0; i < 3; i ++)
  {
    if (PyCallable_Check(g[i]))
    {
      solfec::gravity.gcallback[i] = g[i];
    }
    else if (PyFloat_Check(g[i]))
    {
      solfec::gravity.gvalue[i] = (REAL)PyFloat_AsDouble(g[i]);
    }
    else if (PyLong_Check(g[i]))
    {
      TYPETEST (is_in_map(PyInt_AsLong(g[i]), "SPLINE", solfec::splines));

      solfec::gravity.gspline[i] = PyInt_AsLong(g[i]);
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "a gravity component is neither a float nor an integer nor a callback");
      return NULL;
    }
  }

  Py_RETURN_NONE;
}

/* print gravity */
static PyObject* print_GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("gx", "gy", "gz");

  if (solfec::notrun)
  {
    for (int i = 0; i < 3; i ++)
    {
      if (solfec::gravity.gcallback[i])
      {
	std::cout << "GRAVITY_" << kwl[i] << " = callback" << std::endl;
      }
      else if (solfec::gravity.gspline[i] >= 0)
      {
	std::cout << "GRAVITY_" << kwl[i] << " = " << solfec::gravity.gspline[i] << std::endl;
      }
      else
      {
	std::cout << "GRAVITY_" << kwl[i] << " = " << FMT("%e") << solfec::gravity.gvalue[i] << std::endl;
      }
    }
  }

  Py_RETURN_NONE;
}

/* Retrieve time history */
static PyObject* HISTORY (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Declare output entities */
static PyObject* OUTPUT (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Run simulation */
static PyObject* RUN (PyObject *self, PyObject *args, PyObject *kwds)
{
  solfec::notrun = false;

  /* TODO */

  solfec::velocities.clear();

  Py_RETURN_NONE;
}

/* Delete object */
static PyObject* DELETE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("objnum", "objkind");
  PyObject *objkind;
  size_t objnum;

  PARSEKEYS ("KO", &objnum, &objkind);

  IFIS (objkind,"MESH")
  {
    TYPETEST (is_in_map (objnum, "MESH", solfec::meshes));

    struct mesh &mesh = solfec::meshes[objnum];

    for (std::set<size_t>::iterator it = mesh.restrains.begin(); it != mesh.restrains.end(); it ++)
    {
      solfec::restrains.erase(*it);
    }

    for (std::set<size_t>::iterator it = mesh.prescribes.begin(); it != mesh.prescribes.end(); it ++)
    {
      solfec::prescribes.erase(*it);
    }

    solfec::meshes.erase(objnum);
  }
  ELIF (objkind,"RESTRAIN")
  {
    TYPETEST (is_in_map (objnum, "RESTRAIN", solfec::restrains));

    solfec::restrains.erase(objnum);
  }
  ELIF (objkind,"PRESCRIBE")
  {
    TYPETEST (is_in_map (objnum, "PRESCRIBE", solfec::prescribes));

    solfec::prescribes.erase(objnum);
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid object kind");
    return NULL;
  }

  Py_RETURN_NONE;
}

static PyMethodDef methods [] =
{
  {"ARGV", (PyCFunction)ARGV, METH_VARARGS|METH_KEYWORDS, "Command line arguments"},
  {"RESET", (PyCFunction)RESET, METH_NOARGS, "Reset simulation"},
  {"SPLINE", (PyCFunction)SPLINE, METH_VARARGS|METH_KEYWORDS, "Create linear spline"},
  {"print_SPLINE", (PyCFunction)print_SPLINE, METH_VARARGS|METH_KEYWORDS, "print linear spline"},
  {"MATERIAL", (PyCFunction)MATERIAL, METH_VARARGS|METH_KEYWORDS, "Create material"},
  {"print_MATERIAL", (PyCFunction)print_MATERIAL, METH_VARARGS|METH_KEYWORDS, "print material"},
  {"MESH", (PyCFunction)MESH, METH_VARARGS|METH_KEYWORDS, "Create meshed body"},
  {"print_MESH", (PyCFunction)print_MESH, METH_VARARGS|METH_KEYWORDS, "print meshed body"},
  {"RESTRAIN", (PyCFunction)RESTRAIN, METH_VARARGS|METH_KEYWORDS, "Restrain motion"},
  {"print_RESTRAIN", (PyCFunction)print_RESTRAIN, METH_VARARGS|METH_KEYWORDS, "print restrain"},
  {"PRESCRIBE", (PyCFunction)PRESCRIBE, METH_VARARGS|METH_KEYWORDS, "Prescribe motion"},
  {"print_PRESCRIBE", (PyCFunction)print_PRESCRIBE, METH_VARARGS|METH_KEYWORDS, "print prescribe"},
  {"VELOCITY", (PyCFunction)VELOCITY, METH_VARARGS|METH_KEYWORDS, "Set rigid velocity"},
  {"print_VELOCITIES", (PyCFunction)print_VELOCITIES, METH_NOARGS, "print velocities"},
  {"FRICTION", (PyCFunction)FRICTION, METH_VARARGS|METH_KEYWORDS, "Define friction parameters"},
  {"print_FRICTIONS", (PyCFunction)print_FRICTIONS, METH_NOARGS, "print frictions"},
  {"GRAVITY", (PyCFunction)GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity"},
  {"print_GRAVITY", (PyCFunction)print_GRAVITY, METH_NOARGS, "print gravity"},
  {"HISTORY", (PyCFunction)HISTORY, METH_VARARGS|METH_KEYWORDS, "Retrieve time history"},
  {"OUTPUT", (PyCFunction)OUTPUT, METH_VARARGS|METH_KEYWORDS, "Declare output entities"},
  {"RUN", (PyCFunction)RUN, METH_VARARGS|METH_KEYWORDS, "Run simulation"},
  {"DELETE", (PyCFunction)DELETE, METH_VARARGS|METH_KEYWORDS, "Delete object"},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static PyModuleDef solfec_module = {PyModuleDef_HEAD_INIT, "solfec", NULL, -1, methods, NULL, NULL, NULL, NULL};

static PyObject* PyInit_solfec_module(void)
{
    return PyModule_Create(&solfec_module);
}
#endif

/* interpret an input file (return 0 on success) */
int input (const char *path)
{
  int error, len;
  char *line;

  len = strlen (path);
 
  if (path[len-3] != '.' ||
      path[len-2] != 'p' ||
      path[len-1] != 'y')
  {
    fprintf (stderr, "ERROR: input file does not have '.py' extension!\n");
    fprintf (stderr, "       the input path reads: %s\n", path);
    return 1;
  }
  else solfec::outname.assign(path, len-3);

#if PY_MAJOR_VERSION >= 3
  if (PyImport_AppendInittab("solfec", &PyInit_solfec_module) < 0) return -1;

  Py_Initialize();

  wchar_t** wargv = (wchar_t**)PyMem_Malloc(sizeof(wchar_t*)*solfec::argc);
  for (int i = 0; i < solfec::argc; i++) {
    wchar_t* arg = Py_DecodeLocale(solfec::argv[i], NULL);
    wargv[i] = arg;
  }
  PySys_SetArgv(solfec::argc, wargv);

#else
  Py_Initialize();

  PySys_SetArgv(solfec::argc, solfec::argv);

  if (!Py_InitModule3 ("solfec", methods, "solfec module")) return -1;
#endif

  PyRun_SimpleString ("from solfec import ARGV\n"
                      "from solfec import RESET\n"
                      "from solfec import SPLINE\n"
                      "from solfec import print_SPLINE\n"
                      "from solfec import MATERIAL\n"
                      "from solfec import print_MATERIAL\n"
                      "from solfec import MESH\n"
                      "from solfec import print_MESH\n"
                      "from solfec import RESTRAIN\n"
                      "from solfec import print_RESTRAIN\n"
                      "from solfec import PRESCRIBE\n"
                      "from solfec import print_PRESCRIBE\n"
                      "from solfec import VELOCITY\n"
                      "from solfec import print_VELOCITIES\n"
                      "from solfec import FRICTION\n"
                      "from solfec import print_FRICTIONS\n"
                      "from solfec import GRAVITY\n"
                      "from solfec import print_GRAVITY\n"
                      "from solfec import HISTORY\n"
                      "from solfec import OUTPUT\n"
                      "from solfec import RUN\n"
                      "from solfec import DELETE\n");

  ERRMEM (line = new char [128 + strlen (path)]);
#if PY_MAJOR_VERSION >= 3
  sprintf (line, "exec(open('%s').read())", path);
#else
  sprintf (line, "execfile ('%s')", path);
#endif

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
                                    between WIN32 and UNIX while Python is often provided in binary form */
  delete line;

  return error;
}
