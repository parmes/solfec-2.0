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
#include <string.h>
#include <float.h>
#include "real.h"
#include "err.h"
#include "solfec.h"

static size_t findmarker (std::vector<std::array<REAL,2>>::iterator begin, std::vector<std::array<REAL,2>>::iterator end, REAL xval)
{
  std::vector<std::array<REAL,2>>::iterator low = begin, high = end-1, mid;

  while (low <= high)
  {
    mid = low + (high-low) / 2;

    if (xval >= (*mid)[0] &&
        xval < (*(mid+1))[0])
    {
      return (mid - begin);
    }
    else if (xval < (*mid)[0])
    {
      high = mid - 1;
    }
    else
    {
      low = mid + 1;
    }
  }

  return 0;
}

static size_t findoffset (std::vector<REAL>::iterator begin, std::vector<REAL>::iterator end, REAL xval)
{
  std::vector<REAL>::iterator low = begin, high = end-1, mid = low + (high-low) / 2;

  while (low < mid && mid < high)
  {
    if (xval < (*mid))
    {
      high = mid;
    }
    else if (xval >= (*mid))
    {
      low = mid;
    }

    mid = low + (high-low) / 2;
  }

  return mid - begin;
}

static REAL linterp (std::vector<std::array<REAL,2>>::iterator point, REAL xval)
{
  REAL dt = xval - point[0][0];
  return point[0][1] + (point[1][1]-point[0][1]) * (dt / (point[1][0] - point[0][0]));
}

/* read spline from file */
void spline_from_file (const char *path, int cache, struct spline *spline)
{
  char buf [4096];
  REAL x0, x, y;
  size_t size;
  off_t off0;
  FILE *fp;

  ASSERT ((fp = fopen (path, "r")), "opening SPLINE file [%s] has failed", path);

  spline->cache = cache;

  if (cache) spline->path.assign(path);

  size = 0;
  x = -REAL_MAX;
  while ((off0 = ftello(fp)) >= 0 && fgets(buf, sizeof(buf), fp) != NULL) /* read line */
  {
    if (buf[0] == '#') continue; /* skip comments */

    if (cache == 0 || spline->points.size() <= cache)
    {
      x0 = x;

#if   REALSIZE==4
      if (sscanf (buf, "%f%f", &x, &y) == EOF) break;
#else
      if (sscanf (buf, "%lf%lf", &x, &y) == EOF) break;
#endif

      ASSERT (x > x0, "spline argument must be increasing");

      std::array<REAL,2> temp = {x, y};
      spline->points.push_back (temp);

      size ++;
    }
    else
    {
      x0 = x;

#if   REALSIZE==4
      if (sscanf (buf, "%f%f", &x, &y) == EOF) break;
#else
      if (sscanf (buf, "%lf%lf", &x, &y) == EOF) break;
#endif

      ASSERT (x > x0, "spline argument must be increasing");

      size ++;
    }

    if (spline->points.size() == 1 || size == cache) /* mark cache sized offests in file */
    {
      spline->offset.push_back(off0);
      spline->xval.push_back(x);

      if (size == cache) size = 0;
    }
  }

  fclose (fp);

  ASSERT (spline->points.size() > 0, "spline definition is empty");

  if (cache)
  {
    if (spline->xval.back() < x) /* finished before cache limit */
    {
      spline->offset.push_back(off0);
      spline->xval.push_back(x);
    }
  }
}

/* calculate splane value */
REAL spline_value (struct spline *spline, REAL xval)
{
  REAL lo, hi;

  /* if the spline is partially cached - when out of bounds - reload cache */
  if ((xval < spline->points[0][0] || xval > spline->points.back()[0]) && spline->path.length() > 0)
  {
    char buf [4096];
    int off = findoffset (spline->xval.begin(), spline->xval.end(), xval);
    FILE *fp = fopen (spline->path.c_str(), "r");
    ASSERT (fp, "opening file [%s] for a partially cached SPLINE has failed", spline->path);
    fseeko (fp, spline->offset[off], SEEK_SET);
    spline->points.resize(0);
    while (spline->points.size() <= spline->cache && fgets(buf, sizeof(buf), fp) != NULL) /* read line */
    {
      if (buf[0] == '#') continue; /* skip comments */

      REAL x, y;

#if   REALSIZE==4
      if (sscanf (buf, "%f%f", &x, &y) == EOF) break;
#else
      if (sscanf (buf, "%lf%lf", &x, &y) == EOF) break;
#endif
      std::array<REAL,2> temp = {x, y};
      spline->points.push_back (temp);
    }
    fclose (fp);
  }

  if (xval < spline->points[0][0]) return spline->points[0][1];
  else if (xval > spline->points.back()[0]) return spline->points.back()[1];

  lo = spline->points[spline->marker > 0 ? spline->marker - 1 : spline->marker][0];
  hi = spline->points[spline->marker < spline->points.size() - 1 ? spline->marker + 1 : spline->marker][0];

  if (xval < lo || xval > hi)
  {
    spline->marker = findmarker (spline->points.begin(), spline->points.end(), xval);
  }
  else if (xval >= lo && spline->marker &&
	   xval < spline->points[spline->marker][0]) spline->marker --;
  else if (xval >= spline->points[spline->marker+1][0] &&
	   xval < hi && spline->marker < spline->points.size() - 1) spline->marker ++;

  return linterp (spline->points.begin()+spline->marker, xval);
}
