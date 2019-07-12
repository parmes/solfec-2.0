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
#include <string.h>
#include <stdio.h>
#include "real.h"
#include "err.h"
#include "solfec.h"

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
      #if PY_MAJOR_VERSION >= 3
      else PyList_Append (list, PyUnicode_FromString (argv[i]));
      #else
      else PyList_Append (list, PyString_FromString (argv[i]));
      #endif
    }
  }
  else
  {
    for (int i = 0; i < argc; i ++)
    {
      if (endswith (argv[i], "solfec4")) continue;
      else if (endswith (argv[i], "solfec8")) continue;
      #if PY_MAJOR_VERSION >= 3
      else PyList_Append (list, PyUnicode_FromString (argv[i]));
      #else
      else PyList_Append (list, PyString_FromString (argv[i]));
      #endif
    }
  }

  return list;
}

/* reset simulation */
static PyObject* RESET (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Create linear spline */
static PyObject* SPLINE (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Create material */
static PyObject* MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Create meshed body */
static PyObject* MESH (PyObject *self, PyObject *args, PyObject *kwds)
{
  unsigned long long bodnum = 0;

#if PY_MAJOR_VERSION >= 3
  return PyLong_FromUnsignedLongLong(bodnum);
#else
  return Py_BuildValue("K", bodnum);
#endif
  Py_RETURN_NONE;
}

/* Define contact parameters */
static PyObject* CONTACT (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Restrain motion */
static PyObject* RESTRAIN (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Prescribe moion */
static PyObject* PRESCRIBE (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Set rigid velocity */
static PyObject* VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  Py_RETURN_NONE;
}

/* Set gravity */
static PyObject* GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
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
  Py_RETURN_NONE;
}

static PyMethodDef methods [] =
{
  {"ARGV", (PyCFunction)ARGV, METH_VARARGS|METH_KEYWORDS, "Command line arguments"},
  {"RESET", (PyCFunction)RESET, METH_NOARGS, "Reset simulation"},
  {"SPLINE", (PyCFunction)SPLINE, METH_VARARGS|METH_KEYWORDS, "Create linear spline"},
  {"MATERIAL", (PyCFunction)MATERIAL, METH_VARARGS|METH_KEYWORDS, "Create material"},
  {"MESH", (PyCFunction)MESH, METH_VARARGS|METH_KEYWORDS, "Create meshed body"},
  {"CONTACT", (PyCFunction)CONTACT, METH_VARARGS|METH_KEYWORDS, "Define contact parameters"},
  {"RESTRAIN", (PyCFunction)RESTRAIN, METH_VARARGS|METH_KEYWORDS, "Restrain motion"},
  {"PRESCRIBE", (PyCFunction)PRESCRIBE, METH_VARARGS|METH_KEYWORDS, "Prescribe motion"},
  {"VELOCITY", (PyCFunction)VELOCITY, METH_VARARGS|METH_KEYWORDS, "Set rigid velocity"},
  {"GRAVITY", (PyCFunction)GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity"},
  {"HISTORY", (PyCFunction)HISTORY, METH_VARARGS|METH_KEYWORDS, "Retrieve time history"},
  {"OUTPUT", (PyCFunction)OUTPUT, METH_VARARGS|METH_KEYWORDS, "Declare output entities"},
  {"RUN", (PyCFunction)RUN, METH_VARARGS|METH_KEYWORDS, "Run simulation"},
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
  using namespace solfec; /* output_path */
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
  else output_path.assign(path, len-3);

#if PY_MAJOR_VERSION >= 3
  if (PyImport_AppendInittab("solfec", &PyInit_solfec_module) < 0) return -1;

  Py_Initialize();
#else
  Py_Initialize();

  if (!Py_InitModule3 ("solfec", methods, "solfec module")) return -1;
#endif

  PyRun_SimpleString ("from solfec import ARGV\n"
                      "from solfec import RESET\n"
                      "from solfec import SPLINE\n"
                      "from solfec import MATERIAL\n"
                      "from solfec import MESH\n"
                      "from solfec import CONTACT\n"
                      "from solfec import RESTRAIN\n"
                      "from solfec import PRESCRIBE\n"
                      "from solfec import VELOCITY\n"
                      "from solfec import GRAVITY\n"
                      "from solfec import HISTORY\n"
                      "from solfec import OUTPUT\n"
                      "from solfec import RUN\n");

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
