#ifndef _SEISPP_H_
#define _SEISPP_H_
//@{
// Seismic C++ Library (SEISPP).
// This is a collection of C++ objects used for processing seismic data.  
// Objects define some common seismological concepts including time series
// data, three component seismograms, ensembles of seismograms, velocity
// models, hypocenters, slowness vectors, and others.  
// The library has a strong dependence on Antelope in implementation, but
// the API design was intended to be more general.  
//@}
#ifdef sun
#include <sunmath.h>
#endif
#include <limits.h>
//
// These STL includes are needed for templates in this file
//
#include <vector>
#include <algorithm>
//
// Antelope includes required
//
#include "stock.h"
#include "db.h"
#include "tr.h"
#include "pf.h"
//
// These are glp stuff external to this library
// 
#include "pfstream.h"
//
// Library internal includes files required for the seispp.h  include
// 
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "gclgrid.h"
//@{
// The SEISPP namespace encapsulates the library functions and 
// classes that defined the SEISPP seismic processing library in C++.
// Almost all applications using this library will need a 
// "using namespace SEISPP" line to make the package visible to 
// the compiler and linker.
//@}
namespace SEISPP 
{
using namespace SEISPP;
using namespace std;
//@{
// Turns verbose mode on and off.  
//@}
extern bool SEISPP_verbose;
//
//now the long list of includes for this directory containing
//include files for main object definitions.  This is in keeping
//with recommended style in OOP.
//

//@{
// Applies a time shift called a geometric static to TimeSeries object ts 
// using velocity vel and elevation elev.  This is a simple elev/vel correction.
//@}
void ApplyGeometricStatic(TimeSeries *ts, double vel, double elev);
//@{
// Applies a time shift called a geometric static to TimeSeries object ts 
// using velocity and elevation extracted from the metadata area of ts.  
// The function requires that attributes "elevation" and "surface_velocity"
// be defined in the object.  If these attributes are not defined the 
// data are not altered but a diagnostic is issued to stderr.
//@}
void ApplyGeometricStatic(TimeSeries *ts);
//@{
// Pfstream method for getting a time series object from an input stream.
// Useful only in a pfstream environment which is currently not well developed
// for this application.
//@}
TimeSeries* GetNextTimeSeries(Pfstream_handle *pfh);
//@{
// Companion to GetNextTimeSeries.  
// Useful only in a pfstream environment which is currently not well developed
// for this application.
//@}
TimeSeries *LoadTimeSeriesUsingPf(Pf *pf);
//@{
// Load an ensemble through a pfstream.
//@}
TimeSeriesEnsemble *GetNextEnsemble(Pfstream_handle *pfh,
	 char *tag,MetadataList& mdlist) throw(SeisppError);
//@{
// Load a 3c ensemble through a pfstream.
//@}
ThreeComponentEnsemble *GetNext3cEnsemble(Pfstream_handle *pfh,
	 char *tag,MetadataList& mdlist) throw(SeisppError);
//@{
// Save a 3c seismgram using a pfstream output.
//@}
void PfstreamSave3cseis(ThreeComponentSeismogram *seis,string tag,
	string dir, string dfile, Pfstream_handle *pfh) throw(SeisppError);
//@{
//  Used by TimeSeries constructors to set data gaps using Antelope methods.
//@}
void SetGaps(TimeSeries&,Trsample *,int, string)
		throw(SeisppError);
//@{
// Return direction of particle motion for a P wave with 
// slowness (ux,uy) at a surface with P velocity vp0 and 
// S velocity vs0.
//@}
SphericalCoordinate PMHalfspaceModel(double vp0,double vs0,
	double ux,double uy);
//@{
// Extract one component from a ThreeComponentSeismogram and 
// create a TimeSeries object from it.  
//
//@param tcs is the ThreeComponentSeismogram to convert.
//@param component is the component to extract (0, 1, or 2)
//@param mdl list of metadata to copy to output from input object.
//@}
TimeSeries *ExtractComponent(ThreeComponentSeismogram& tcs,int component,
        MetadataList& mdl);
//@{
// Extract one component from a ThreeComponentSeismogram and 
// create a TimeSeries object from it.  
// Similar to overloaded function of same name, but all metadata from
// parent is copied to the output.
//
//@param tcs is the ThreeComponentSeismogram to convert.
//@param component is the component to extract (0, 1, or 2)
//@}
TimeSeries *ExtractComponent(ThreeComponentSeismogram& tcs,int component);
//@{
// Returns a new seismogram in an arrival time reference.
// An arrival time reference means that the time is set to relative and 
// zero is defined as an arrival time extracted from the metadata area of
// the object.  The key used to extract the arrival time used for the
// conversion is passed as a variable as this requires some flexibility.
// To preserve the absolute time standard in this conversion the 0 time
// computed from the arrival time field is used to compute the absolute
// time of the start of the output seismogram as atime+t0.  This result
// is stored in the metadata field keyed by the word "time".  This allows
// one to convert the data back to an absolute time standard if they so
// desire, but it is less flexible than the input key method.  
//
//@throws SeisppError for errors in extracting required information from metadata area.
//
//@param din  is input seismogram
//@param key is the metadata key used to find the arrival time to use as a reference.
//@param tw is a TimeWindow object that defines the window of data to extract around
//    the desired arrival time.
//@}
auto_ptr<ThreeComponentSeismogram> ArrivalTimeReference(ThreeComponentSeismogram& din,
	string key, TimeWindow tw);
//@{
// Returns a gather of ThreeComponentSeismograms in an arrival time reference fram.
// An arrival time refernce means that the time is set to relative and 
// zero is defined as an arrival time extracted from the metadata area of
// each member object.
//
//@throws SeisppError for errors in extracting required information from metadata area.
//
//@param din  is input gather
//@param key is the metadata key used to find the arrival time to use as a reference.
//@param tw is a TimeWindow object that defines the window of data to extract around
//    the desired arrival time.
//@}
auto_ptr<ThreeComponentEnsemble> ArrivalTimeReference(ThreeComponentEnsemble& din,
	string key, TimeWindow tw);
//@{
// Returns a new TimeSeries seismogram in an arrival time reference.
// An arrival time reference means that the time is set to relative and 
// zero is defined as an arrival time extracted from the metadata area of
// the object.  The key used to extract the arrival time used for the
// conversion is passed as a variable as this requires some flexibility.
// To preserve the absolute time standard in this conversion the 0 time
// computed from the arrival time field is used to compute the absolute
// time of the start of the output seismogram as atime+t0.  This result
// is stored in the metadata field keyed by the word "time".  This allows
// one to convert the data back to an absolute time standard if they so
// desire, but it is less flexible than the input key method.  
//
//@throws SeisppError for errors in extracting required information from metadata area.
//
//@param din  is input seismogram
//@param key is the metadata key used to find the arrival time to use as a reference.
//@param tw is a TimeWindow object that defines the window of data to extract around
//    the desired arrival time.
//@}
auto_ptr<TimeSeries> ArrivalTimeReference(TimeSeries& din,
	string key, TimeWindow tw);
//@{
// Returns a gather of TimeSeries objects in an arrival time reference frame.
// An arrival time refernce means that the time is set to relative and 
// zero is defined as an arrival time extracted from the metadata area of
// each member object.
//
//@throws SeisppError for errors in extracting required information from metadata area.
//
//@param din  is input gather
//@param key is the metadata key used to find the arrival time to use as a reference.
//@param tw is a TimeWindow object that defines the window of data to extract around
//    the desired arrival time.
//@}
auto_ptr<TimeSeriesEnsemble> ArrivalTimeReference(TimeSeriesEnsemble& din,
	string key, TimeWindow tw);


//@{
// Bombproof low level write routine for a vector of doubles.  
// Uses fwrite to write vector x to the file dir+"/"+dfile.
//
//@throws SeisppError object if there are problems saving data to requested file.
//@param x vector of data to be saved.
//@param n length of vector x
//@param dir directory to place file.  If empty assumes current directory.
//@param dfile file name 
//@}
long int vector_fwrite(double *x,int n, string dir, string dfile) throw(SeisppError);
//@{
// Bombproof low level write routine for a vector of doubles.  
// Uses fwrite to write vector x to the file fname
//
//@throws SeisppError object if there are problems saving data to requested file.
//@param x vector of data to be saved.
//@param n length of vector x
//@param fname file name 
//@}
long int vector_fwrite(double *x,int n, string fname) throw(SeisppError);
//@{
// Bombproof low level write routine for a vector of floats.  
// Uses fwrite to write vector x to the file dir+"/"+dfile.
//
//@throws SeisppError object if there are problems saving data to requested file.
//@param x vector of data to be saved.
//@param n length of vector x
//@param dir directory to place file.  If empty assumes current directory.
//@param dfile file name 
//@}
long int vector_fwrite(float *x,int n, string dir, string dfile) throw(SeisppError);
//@{
// Bombproof low level write routine for a vector of floats.  
// Uses fwrite to write vector x to the file fname
//
//@throws SeisppError object if there are problems saving data to requested file.
//@param x vector of data to be saved.
//@param n length of vector x
//@param fname file name 
//@}
long int vector_fwrite(float *x,int n, string fname) throw(SeisppError);
//@{
// Save the data in a TimeSeries object to a database.
// This function works only with an Antelope (Datascope) database but the
// design is aimed to be schema independent.  That is, most raw 
// earthquake seismology data is indexed with a table defined in css3.0
// called wfdisc.  This function will work with a css3.0 wfdisc, but it
// will work with any other table as well provided you set up the
// interface correctly.  This is done through the MetadataList object
// which tells the function what attributes are to be saved to the
// output database along with the time series data.  
//
// A TimeSeries object contains a Metadata object it acquires by 
// inheritance.  The Metadata area is assumed to contain attributes
// listed in the MetadataList object passed to this function.  The
// basic algorithm is that the list of metadata in mdl are processed in order.
// They are translated to the database namespace using the AttributeMap
// object am and pushed to an output record using the Datascope dbputv
// function one attribute at a time.  The data are saved to files
// whose name and location are driven by two (frozen) standard names
// extracted from the metadata area:  dir and dfile.  The filename
// for output is created as dir+"/"+dfile or simply dfile if dir
// is not defined (assumes current directory).  
//
// This function is dogmatic about four database output names.  
// It always translates it's internal t0 to be a "time" database
// attribute,  the ns variable is saved as "nsamp", the sample
// interval (dt) is always converted to 1/dt and called "samprate",
// and an endtime (computed as this->endtime()) is computed and
// saved as the database attribute "endtime".   These are the css3.0
// name conventions and I chose to leave them as frozen names.  
// Note also that if the "live" boolean in the object is set false
// this function silently returns immediately doing nothing.
//
//@throws SeisppError object if there are any problems saving the data or 
//    writing attributes into the database.
//
//@returns -1 if live is false, record number of added row otherwise
//
//@param ts is the TimeSeries object to be saved.
//@param db is a Datascope database pointer.  It need only point at a valid
//   open database.
//@param table is the name of the table to index this time series data
//   (e.g. "wfdisc").
//@param md  is the list of metadata to be dumped to the database as described above.
//@param am is a mapping operator that defines how internal names are to be mapped
//    to database attribute names and tables.  
//@}
int dbsave(TimeSeries& ts,Dbptr db,string table, MetadataList& md, AttributeMap& am)
		throw(SeisppError);
//@{
// Save the data in a ThreeComponentSeismogram object to a database.
// This function works only with an Antelope (Datascope) database but the
// design is aimed to be schema independent.  That is, most raw 
// earthquake seismology data is indexed with a table defined in css3.0
// called wfdisc.  This function will work with a css3.0 wfdisc, but it
// will work with any other table as well provided you set up the
// interface correctly.  This is done through the MetadataList object
// which tells the function what attributes are to be saved to the
// output database along with the time series data.  
//
// A ThreeComponentSeismogram object contains a Metadata object it acquires by 
// inheritance.  The Metadata area is assumed to contain attributes
// listed in the MetadataList object passed to this function.  The
// basic algorithm is that the list of metadata in mdl are processed in order.
// They are translated to the database namespace using the AttributeMap
// object am and pushed to an output record using the Datascope dbputv
// function one attribute at a time.  The data are saved to files
// whose name and location are driven by two (frozen) standard names
// extracted from the metadata area:  dir and dfile.  The filename
// for output is created as dir+"/"+dfile or simply dfile if dir
// is not defined (assumes current directory).  
//
// This function is dogmatic about four database output names.  
// It always translates it's internal t0 to be a "time" database
// attribute,  the ns variable is saved as "nsamp", the sample
// interval (dt) is always converted to 1/dt and called "samprate",
// and an endtime (computed as this->endtime()) is computed and
// saved as the database attribute "endtime".   These are the css3.0
// name conventions and I chose to leave them as frozen names.  
// Note also that if the "live" boolean in the object is set false
// this function silently returns immediately doing nothing.
//
// This function differs in a significant way from an overloaded function
// with the same name.  The other has a "chanmap" argument to tell the 
// function how to split up the 3 components into channel codes.  This
// function takes a very different approach and saves data by dumping
// the internal 3xns matrix as the basic output data series. As a result
// this function will write ONE AND ONLY ONE DATABASE ROW PER OBJECT.
// This means somewhat by definition that the output table CANNOT be
// wfdisc if this function is called.  Consequently, this routine will
// throw an exception and do nothing if table=="wfdisc".
//
//@throws SeisppError object if there are any problems saving the data or 
//    writing attributes into the database.
//
//@returns -1 if live is false, record number of added row otherwise
//
//@param ts is the TimeSeries object to be saved.
//@param db is a Datascope database pointer.  It need only point at a valid
//    open database.
//@param table is the name of the table to index this time series data
//   (e.g. "wfdisc").
//@param md  is the list of metadata to be dumped to the database as described above.
//@param am is a mapping operator that defines how internal names are to be mapped
//    to database attribute names and tables.  
//@}
int dbsave(ThreeComponentSeismogram& ts,Dbptr db,string table, 
	MetadataList& md, AttributeMap& am);
//@{
// Save the data in a ThreeComponentSeismogram object to a database.
// This function works only with an Antelope (Datascope) database but the
// design is aimed to be schema independent.  That is, most raw 
// earthquake seismology data is indexed with a table defined in css3.0
// called wfdisc.  This function will work with a css3.0 wfdisc, but it
// will work with any other table as well provided you set up the
// interface correctly.  This is done through the MetadataList object
// which tells the function what attributes are to be saved to the
// output database along with the time series data.  
//
// A ThreeComponentSeismogram object contains a Metadata object it acquires by 
// inheritance.  The Metadata area is assumed to contain attributes
// listed in the MetadataList object passed to this function.  The
// basic algorithm is that the list of metadata in mdl are processed in order.
// They are translated to the database namespace using the AttributeMap
// object am and pushed to an output record using the Datascope dbputv
// function one attribute at a time.  The data are saved to files
// whose name and location are driven by two (frozen) standard names
// extracted from the metadata area:  dir and dfile.  The filename
// for output is created as dir+"/"+dfile or simply dfile if dir
// is not defined (assumes current directory).  
//
// This function is dogmatic about four database output names.  
// It always translates it's internal t0 to be a "time" database
// attribute,  the ns variable is saved as "nsamp", the sample
// interval (dt) is always converted to 1/dt and called "samprate",
// and an endtime (computed as this->endtime()) is computed and
// saved as the database attribute "endtime".   These are the css3.0
// name conventions and I chose to leave them as frozen names.  
// Note also that if the "live" boolean in the object is set false
// this function silently returns immediately doing nothing.
//
// The chanmap and output_to_standard variables control how the 
// data are saved externally.  If output_as_standard is set true
// (highly recommended in general)  the data are restored (if necessary)
// to standard 3c data geometry (ew,ns,z) before being written to 
// output.  In that case vang and hang are set accordingly in 
// case the output algorithm requires these to be stored.  
// The components are then extraced from the 3c object one by 
// one and written in three successive database rows with the
// channel code ("chan" attribute in css3.0) derived from the
// chanmap array (chanmap[0]=channel name for component 0,
// chanmap[1]=component 1, and chanmap[2]=component 2).
//
//@throws SeisppError object if there are any problems saving the data or 
//    writing attributes into the database.
//
//@returns -1 if live is false, record number of added row otherwise
//
//@param ts is the TimeSeries object to be saved.
//@param db is a Datascope database pointer.  It need only point at a valid
//    open database.
//@param table is the name of the table to index this time series data
//   (e.g. "wfdisc").
// @param md  is the list of metadata to be dumped to the database as described above.
// @param am is a mapping operator that defines how internal names are to be mapped
//    to database attribute names and tables.  
//@param chanmap is a set of channel names to map each component to channel code (see above)
//@param output_as_standard when true forces data to be converted to ew,ns, z system
//@}
int dbsave(ThreeComponentSeismogram& ts,Dbptr db,
	string table, MetadataList& md, 
	AttributeMap& am, vector<string>chanmap,bool output_as_standard);
//@{
// Builds a new ensemble of members that satisfy unix regular expression
// for sta and chan attributes passed as sta_expr and chan_expr.
//
// @param parent original ensemble to be subsetted
// @param sta_expr unix regular expression to apply to sta Metadata
//    attribute
// @param chan_expr unix regular expression to apply to chan Metadata 
//    attribute
//
//@author Gary L. Pavlis
//@}
TimeSeriesEnsemble *StaChanRegExSubset(TimeSeriesEnsemble& parent,
        string sta_expr, string chan_expr);

//@{
// Extracts a requested time window of data from a parent TimeSeries object.
//
// It is common to need to extract a smaller segment of data from a larger
// time window of data.  This function accomplishes this in a nifty method that
// takes advantage of the methods contained in the BasicTimeSeries object for
// handling time and data gaps.
//
//@returns new TimeSeries object derived from  parent but windowed by input
//      time window range.
//
//@throws SeisppError object if the requested time window does not overlap data
//
//@param parent is the larger TimeSeries object to be windowed
//@param tw defines the data range to be extracted from parent.
//@author Gary L. Pavlis
//@}
TimeSeries WindowData(TimeSeries& parent, TimeWindow& tw);

//@{
// Extracts a requested time window of data from a parent ThreeComponentSeismogram object.
//
// It is common to need to extract a smaller segment of data from a larger
// time window of data.  This function accomplishes this in a nifty method that
// takes advantage of the methods contained in the BasicTimeSeries object for
// handling time and data gaps.
//
//@returns new ThreeComponentSeismogram object derived from  parent but windowed by input
//      time window range.
//
//@throws SeisppError object if the requested time window does not overlap data
//
//@param parent is the larger ThreeComponentSeismogram object to be windowed
//@param tw defines the data range to be extracted from parent.
//@author Gary L. Pavlis
//@}
ThreeComponentSeismogram WindowData(ThreeComponentSeismogram& parent, TimeWindow& tw);

//@{ Extract a specified time window from an ensemble.
// The seispp library defines a fairly generic ensemble object that
// uses an STL vector container to hold an array of objects 
// (currently TimeSeries or ThreeComponentSeismogram objects) with the
// generic symbol "member".  This template applies a comparable 
// WindowData function to each member of the ensemble returning a
// new ensemble cut to the specified window.
//
//@throws SeisppError exception if TimeWindow is not consistent
// with input data.
//
//@param  parent input ensemble 
//@param tw TimeWindow to cut parent to produce output.
//
//@returns new ensemble T as an auto_ptr cut to desired window.
//@}
//template <class T> auto_ptr<T>WindowData(T& parent, TimeWindow& tw);

//@{
// Sorts an ensemble by station:channel.  
// In earthquake seismic data processing sorting data by station name
// and channel code is a very common operation.  This procedure implements
// this using the STL sort algorithm.  It depends upon the metadata
// fields keyed by "sta" and "chan" being defined in each member of the 
// input ensemble.
//
//@param ensemble is the input data to be sorted.  The STL algorithm
// sorts this in place so the result is altered.  Any indices using
// subscripts will no longer be valid on exit.
//@}
void StaChanSort(TimeSeriesEnsemble& ensemble);
//@{ Generic routine to compute a median.
// This template can be used to compute the median of a vector
// of objects for any class which has the default comparison
// operator needed by STL sort.  It is most likely to be used
// for standard numerical types (int, float, double) where it
// is guaranteed to work.  If T is more exotic, you need to understand
// the rules of what sort expects.  
//
//@param x - input STL vector of data to be compute median.
//   x is not altered.  We make a copy of this input and sort it.
//   Not the best algorithm if the sorted output is desired
//@}
template <class T> T median(vector<T>& x)
{
	int count=x.size();
	if(count<=1)return(x[0]);
	int medposition=count/2;
	T result;
	vector<T> copyx(x);
	sort(copyx.begin(),copyx.end());
	if(count%2)
		result=copyx[medposition];
	else
		result=(copyx[medposition]+copyx[medposition+1])/2.0;
	return (result);
}
//@{
// Aligns an ensemble of data by moveout.  
//
// A common theme in multichannel processing is alignment of data
// by a static shift.  That is, the entire trace is shifted by a fixed
// time amount.  This generic method does this for ensembles with a
// vector of objects held in a variable called "member".  It uses
// a special keyword to extract a "moveout" value from the data objects
// metadata and applies this to forever shift the 0 time reference for
// that data object.  For this reason it always returns a copy of the
// original data.
//
//@param ensemble - input data ensemble.
//
//@returns copy of ensemble but with t0 values modified by moveout.
//@}
template <class Tensemble> Tensemble MoveoutTimeShift(Tensemble& d)
{
	Tensemble dshift(d);
	try {
		for(int i=0;i<d.member.size();++i)
		{
			double tshift;
			if(d.member[i].live)
			{
				tshift=d.member[i].get_double(moveout_keyword);
				dshift.member[i].t0-=tshift;
			}
		}
	} catch (MetadataGetError mde)
	{
		throw SeisppError(string("MoveoutTimeShift:  moveout keyword")
			+ moveout_keyword
			+string(" not defined for one or members of input ensemble") );
	}
	return (dshift);
}

//@{
// Convert a velocity model in flattened earth geometry to true geometry.
//
// Multiple programs exist in seismology that use a flattening transformation
// to convert spherical geometry into a regular grid.  The main reason
// for this in all examples I know of is to allow standard finite difference
// algorithms to used.  This is essential since standard finite differences
// demand regular grid geometries and the flattening transformation provides
// a convenient way to turn a spherical shell problem into a standard 
// finite different grid.  A GCLgrid3d object is aimed at getting the
// true geometry correct, however, so using velocity models defined on 
// a flattened coordinate system requires a conversion.  This function
// does this for a complete grid of points.  This algorithm assumes 
// a 3D model has been loaded into a GCLgrid3d object, BUT the coordinates
// and velocities stored in the grid have not been corrected by the
// inverse flattening transformation.  That is, the assumption is that
// the original model was defined in on flattened coordinate system but
// the points were converted to lat,lon, and (flattened) depth in building
// the input GCLscalarfield3d object.  The algorithm used here then just
// passes through the grid doing two things:  (1) the depth() method is 
// called for the parent grid and the uflatz (see gclgrid.h) function is
// called to convert this to true depth; and (2) the velocity value at each
// point is altered by in inverse flattening transformation.  Thus the 
// net effect is a depth dependent change in velocity values and a nonlinear
// distortion of the depth of each point in the grid.  
//
//@param vmodel is the parent grid.  It is altered in place to convert
//  flattened to true geometry as described above.
//
//@author Gary L. Pavlis
//@}
void ConvertFlatModelToSpherical(GCLscalarfield3d& vmodel);
		
}
#endif
