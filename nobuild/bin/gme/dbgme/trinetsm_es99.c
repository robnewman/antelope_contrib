/* Copyright (c) 2004 Boulder Real Time Technologies, Inc. */
/* All rights reserved. */
/*                                                                     
/* Written by Dr. Kent Lindquist, Lindquist Consulting, Inc. */
/*
/* This software may be used freely in any way as long as */
/* the copyright statement above is not removed. */

#include "dbgme.h"

#define DEFAULT_PGA_CALC_UNITS "gravity"
#define DEFAULT_PGV_CALC_UNITS "cm/sec"
#define DEFAULT_PGA_OUTPUT_UNITS "centigravity"
#define DEFAULT_PGV_OUTPUT_UNITS "cm/sec"
#define DEFAULT_WFMGME_PVA_UNITS "milligravity"
#define DEFAULT_WFMGME_PVV_UNITS "nm/sec"
#define DEFAULT_MMI_PGA_UNITS "cm/sec/sec"
#define DEFAULT_MMI_PGV_UNITS "cm/sec"

enum Wfmode { WFMGME, WFMEAS };

struct mypdata_ {
	char	chan[STRSZ];
	char	units[STRSZ];
	char	sta[STRSZ];
};

struct mmi_params_ {
	CGGrid	*cgg_pga;
	CGGrid	*cgg_pgv;
	double	low_coeff;
	double	low_offset;
	double	low_cutoff_mmi;
	double	pga_coeff;
	double	pga_offset;
	double	pgv_coeff;
	double	pgv_offset;
	double	pga_cutoff_mmi;
};

static double  *Ampfactor_cutoffs_g = 0;
static double **Ampfactors = 0;
static int Ncutoffs = 0;
static int Nvelocities = 0;

void
free_ampfactors( void )
{
	int	i;

	if( Ampfactor_cutoffs_g ) {

		free( Ampfactor_cutoffs_g );

		Ampfactor_cutoffs_g = 0;
	}

	if( Ampfactors ) {

		for( i = 0; i < Nvelocities; i++ ) {
		
			if( Ampfactors[i] ) {

				free( Ampfactors[i] );
			}
		}

		free( Ampfactors );

		Ampfactors = 0;
	}

	Nvelocities = 0;

	return;
}

int 
setup_ampfactors( Pf *pf )
{
	char	*pga_sitecorr_cutoffs_g;
	char	*arow;
	Tbl	*pga_sitecorr_table;
	Tbl	*cutoffs;
	Tbl	*amps;
	int	icut, ivel;

	pga_sitecorr_cutoffs_g = pfget_string( pf, "pga_sitecorr_cutoffs_g" );
	pga_sitecorr_table = pfget_tbl( pf, "pga_sitecorr_table" );

	if( pga_sitecorr_cutoffs_g == NULL ) {
		elog_complain( 0, "trinetsm_es99: need value for "
			     "pga_sitecorr_cutoffs_g\n" );
		return -1;
	}

	if( pga_sitecorr_table == NULL ) {
		elog_complain( 0, "trinetsm_es99: need value for "
			     "pga_sitecorr_table\n" );
		return -1;
	}

	pga_sitecorr_cutoffs_g = strdup( pga_sitecorr_cutoffs_g );
	cutoffs = split( pga_sitecorr_cutoffs_g, ' ' );

	Ncutoffs = maxtbl( cutoffs );

	if( Ncutoffs < 1 ) {
		elog_complain( 0, "trinetsm_es99: need at least one value in "
			     "pga_sitecorr_cutoffs_g\n" );
		free( pga_sitecorr_cutoffs_g );
		freetbl( cutoffs, 0 );
		return -1;
	}

	allot( double *, Ampfactor_cutoffs_g, Ncutoffs );

	for( icut = 0; icut < Ncutoffs; icut++ ) {

		Ampfactor_cutoffs_g[icut] = atof( gettbl( cutoffs, icut ) );
	}

	free( pga_sitecorr_cutoffs_g );
	freetbl( cutoffs, 0 );

	Nvelocities = maxtbl( pga_sitecorr_table );

	if( Nvelocities < 1 ) {
		elog_complain( 0, "trinetsm_es99: need at least one row in "
			     "pga_sitecorr_table\n" );
		free_ampfactors();
		return -1;
	}

	allot( double **, Ampfactors, Nvelocities );

	for( ivel = 0; ivel < Nvelocities; ivel++ ) {
		
		allot( double *, Ampfactors[ivel], Ncutoffs + 1 );

		arow = gettbl( pga_sitecorr_table, ivel );

		arow = strdup( arow );
		amps = split( arow, ' ' );

		if( maxtbl( amps ) != Ncutoffs + 1 ) {
			elog_complain( 0,
				"trinetsm_es99: wrong number of entries "
				"in pga_sitecorr_table row <%s>: should be "
				"\"minvel ampfactor ampfactor ....\" "
				"(need %d + minvel = %d entries, have %d)\n",
				arow, Ncutoffs, Ncutoffs+1, maxtbl( amps ));
			free( arow );
			freetbl( amps, 0 );
			free_ampfactors();
			return -1;
		} 
			
	   	for( icut = 0; icut < Ncutoffs + 1; icut++ ) {
			Ampfactors[ivel][icut] = atof( gettbl( amps, icut ) );	
	   	}

		free( arow );
		freetbl( amps, 0 );
	}

	freetbl( pga_sitecorr_table, 0 );

	return 0;
}

double
pga_ampfactor( double sitevel, double accel )
{
	int	icutoff;
	int	ivel;

	for( ivel = Nvelocities-1; ivel >= 0; ivel-- ) {

		if( sitevel >= Ampfactors[ivel][0] ) break;
	}

	for( icutoff = Ncutoffs; icutoff >= 1; icutoff-- ) {

		if( ( accel > Ampfactor_cutoffs_g[icutoff-1] ) ||
		    icutoff == 1 ) {

			break;
		}
	}

	return Ampfactors[ivel][icutoff];
}

double 
apply_ampfactor( CGGrid *cgg, int ix, int iy, void *sitegridptr )
{
	CGGrid	*sitegrid = (CGGrid *) sitegridptr;
	double	ampfactor;

	if( ! strcmp( cgg->units, DEFAULT_PGA_CALC_UNITS ) ) {

		ampfactor = pga_ampfactor( sitegrid->data[ix][iy], 
					   cgg->data[ix][iy] );

	} else {

		ampfactor = pga_ampfactor( sitegrid->data[ix][iy], 0 );
	}

	return cgg->data[ix][iy] * ampfactor;
}

void *
fill_pdata( Dbptr db )
{
	struct mypdata_ *mp;

	allot( struct mypdata_ *, mp, 1 );

	dbgetv( db, 0, "wfmeas.sta", mp->sta, 0 );
	dbgetv( db, 0, "wfmeas.chan", mp->chan, 0 );
	dbgetv( db, 0, "units1", &mp->units, 0 );

	return mp;
}

void *
dup_pdata( void *pdata_in ) 
{
	struct mypdata_ *mp = (struct mypdata_ *) pdata_in;
	struct mypdata_ *dup;

	allot( struct mypdata_ *, dup, 1 );

	strncpy( dup->sta, mp->sta, STRSZ );
	strncpy( dup->units, mp->units, STRSZ );

	return dup;
}

int
cmp_string( char **a, char **b, void *private )
{
	return strcmp( *a, *b );
}

double 
compute_mmi( CGGrid *cgg, int ix, int iy, void *mp_ptr )
{
	struct mmi_params_ *mp = (struct mmi_params_ *) mp_ptr;
	double	mmi;
	double	mmi_pga;
	double	mmi_pgv;
	double	weight;

	if( is_nan( mp->cgg_pga->data[ix][iy] ) ) {

		elog_complain( 0, "compute_mmi: NaN PGA value at %f %f, "
				  "setting MMI to 0\n",
				  CGGRID_X( mp->cgg_pga, ix ),
				  CGGRID_Y( mp->cgg_pga, iy ) );
		return 0;
	}

	mmi = mp->low_coeff * log10( mp->cgg_pga->data[ix][iy] ) + 
	                                                mp->low_offset;

	if( mmi < 0 ) {

		return 0;

	} else if( mmi <= mp->low_cutoff_mmi ) {
		
		return mmi;
	}

	if( is_nan( mp->cgg_pgv->data[ix][iy] ) ) {

		elog_complain( 0, "compute_mmi: NaN PGV value at %f %f, "
				  "setting MMI to 0\n",
				  CGGRID_X( mp->cgg_pgv, ix ),
				  CGGRID_Y( mp->cgg_pgv, iy ) );
		return 0;
	}

	mmi_pga = mp->pga_coeff * log10( mp->cgg_pga->data[ix][iy] ) + 
							mp->pga_offset;

	mmi_pgv = mp->pgv_coeff * log10( mp->cgg_pgv->data[ix][iy] ) + 
							mp->pgv_offset;

	/* Just in case: */
	if( mmi_pga <= mp->low_cutoff_mmi ) {
		
		return mmi_pga;
	}

	if( mmi_pga > mp->pga_cutoff_mmi ) {
		
		return mmi_pgv;
	}

	weight = ( mmi_pga - mp->low_cutoff_mmi ) / 
		 ( mp->pga_cutoff_mmi - mp->low_cutoff_mmi );

	mmi = (1-weight) * mmi_pga + weight * mmi_pgv;

	return mmi; 
}

/* Note: cggridval_mult has been moved to libcgeom and is being
 * kept here temporarily as scaffolding to sidestep version skew */
static double
cggridval_mult_static( CGGrid *cgg, int ix, int iy, void *factor_ptr )
{
	double	factor = *( (double *) factor_ptr );

	cgg->data[ix][iy] *= factor;

	return cgg->data[ix][iy];
}

double
nearest_deg( CGPointset *cgps, double lon, double lat )
{
	CGPoint	*cgpt;
	int	ipt;
	double	nearest = 99999999999;
	double	mydist, myaz;

	for( ipt = 0; ipt < cgpointset_cnt( cgps ); ipt++ ) {

		cgpt = cgpointset_getpoint( cgps, ipt );

		dist( rad( cgpt->y ), rad( cgpt->x ), 
		      rad( lat ), rad( lon ),
		      &mydist, &myaz );

		mydist = deg( mydist );

		if( mydist < nearest ) {
			nearest = mydist;
		}
	}

	return nearest;
}

CGPoint *
new_phantom( double lon, double lat, double dist_km, double mag, 
	     double const_coeff, double r_offset, double r_coeff, 
	     double mag_coeff, double P_coeff, double P,
	     double S_coeff, double S )
{
	CGPoint *cgpt;
	double	r;
	double	dataval;

	r = sqrt( dist_km * dist_km + r_offset * r_offset );

	dataval = const_coeff;
	dataval += mag_coeff * mag;
	dataval -= log10( r );
	dataval -= r_coeff * r;
	dataval += P_coeff * P;
	dataval += S_coeff * S;
	dataval = pow( 10, dataval );

	cgpt = cgpoint_new( lon, lat, dataval ); 

	return cgpt;
}

int
trinetsm_es99_mmi( Dbptr db, Pf *pf )
{
	CGGrid	*cgg;
	CGGrid	*cgg_pgv_orig;
	Dbptr	dbqgrid;
	Dbptr	dbpga;
	Dbptr	dbpgv;
	double	latc;
	double	lonc;
	double	ml;
	double	conversion_factor;
	struct 	mmi_params_ mp;
	int	orid;
	int	nrecs;
	char 	*recipe_name;
	char 	*pga_recipe;
	char 	*pgv_recipe;
	char	pga_filename[FILENAME_MAX];
	char	pgv_filename[FILENAME_MAX];
	char	*output_units = "mmi";
	char	*output_file;
	char	*qgridfmt;
	char	*qgridtype;
	char	*auth;
	char	grid_name[STRSZ];
	char	cmd[STRSZ];
	char	expr[STRSZ];
	char	returned_units[STRSZ];
	Tbl	*commands;
	FILE	*fp;
	int	flags = 0;
	int	rc;

	qgridtype = pfget_string( pf, "qgridtype" );

	if( strcmp( qgridtype, "mmi" ) ) {

		elog_complain( 0, 
			  "trinetsm_es99_mmi: unknown qgridtype '%s' " 
			  "requested in parameter file, Can't continue\n",
			  qgridtype );
		return -1;
	}

	if( Verbose ) {
		elog_notify( 0, "trinetsm_es99: Calculating an " 
				"%s grid\n", qgridtype );
	}

	mp.low_coeff = pfget_double( pf, "low_coeff" );
	mp.low_offset = pfget_double( pf, "low_offset" );
	mp.low_cutoff_mmi = pfget_double( pf, "low_cutoff_mmi" );
	mp.pga_coeff = pfget_double( pf, "pga_coeff" );
	mp.pga_offset = pfget_double( pf, "pga_offset" );
	mp.pgv_coeff = pfget_double( pf, "pgv_coeff" );
	mp.pgv_offset = pfget_double( pf, "pgv_offset" );
	mp.pga_cutoff_mmi = pfget_double( pf, "pga_cutoff_mmi" );

	recipe_name = pfget_string( pf, "recipe_name" );
	pga_recipe = pfget_string( pf, "pga_recipe" );
	if( pga_recipe == NULL ) {
		elog_complain( 0, 
			"trinetsm_es99_mmi: Need pga_recipe in "
			"parameter file\n" );
		return -1;
	}
	pgv_recipe = pfget_string( pf, "pgv_recipe" );
	if( pgv_recipe == NULL ) {
		elog_complain( 0, 
			"trinetsm_es99_mmi: Need pgv_recipe in "
			"parameter file\n" );
		return -1;
	}
	output_file = pfget_string( pf, "output_file" );
	qgridfmt = pfget_string( pf, "qgridfmt" );
	auth = pfget_string( pf, "auth" );

	db.record = 0;
	dbgetv( db, 0, "origin.lat", &latc, 
		       "origin.lon", &lonc,
		       "ml", &ml, 
		       "orid", &orid,
		       0 );

	sprintf( grid_name, "orid_%d", orid );

	sprintf( cmd, "dbsubset orid == %d", orid );

	commands = strtbl( "dbopen qgrid", cmd, 0 );
	dbqgrid = dbprocess( db, commands, 0 );
	freetbl( commands, 0 );

	dbquery( dbqgrid, dbRECORD_COUNT, &nrecs );

	if( nrecs < 2 ) {
		elog_complain( 0, 
			"trinetsm_es99_mmi: Couldn't find enough input qgrids\n" );
		return -1;
	}

	dbpga = dbqgrid;
	dbpga.record = -1;
	sprintf( expr, "recipe == \"%s\"", pga_recipe );
	rc = dbfind( dbpga, expr, 0, 0 );

	if( rc < 0 ) {

		elog_complain( 0, 
			"trinetsm_es99_mmi: Couldn't find a '%s' "
			"grid for orid %d\n",
			pga_recipe, orid );

		return -1;

	} else {

		dbpga.record = rc;

		rc = dbfilename( dbpga, pga_filename );

		if( rc < 0 ) {

			elog_complain( 0, 
				"trinetsm_es99_mmi: pga grid file '%s' "
				"for orid %d does not exist\n",
				pga_recipe, orid );

			return -1;
		}
	}

	dbpgv = dbqgrid;
	dbpgv.record = -1;
	sprintf( expr, "recipe == \"%s\"", pgv_recipe );
	rc = dbfind( dbpga, expr, 0, 0 );

	if( rc < 0 ) {

		elog_complain( 0, 
			"trinetsm_es99_mmi: Couldn't find a '%s' "
			"grid for orid %d\n",
			pgv_recipe, orid );

		return -1;

	} else {

		dbpgv.record = rc;

		rc = dbfilename( dbpgv, pgv_filename );

		if( rc < 0 ) {

			elog_complain( 0, 
				"trinetsm_es99_mmi: pgv grid file '%s' "
				"for orid %d does not exist\n",
				pgv_recipe, orid );

			return -1;
		}
	}

	if( Verbose ) {

		elog_notify( 0, 
			"trinetsm_es99_mmi: Combining pga grid %s and "
			"pgv grid %s\n", pga_filename, pgv_filename );
	}

	fp = fopen( pga_filename, "r" );
	mp.cgg_pga = cggrid_read( fp );	
	fclose( fp );

	if( mp.cgg_pga == (CGGrid *) NULL ) {

		elog_complain( 0,
			  "trinetsm_es99: Error reading pga file "
			  "%s.\n", pga_filename );
		return -1;
	}

	rc = units_convert( 1.0, mp.cgg_pga->units, DEFAULT_MMI_PGA_UNITS, 
			    &conversion_factor, returned_units );

	if( rc == 1 ) {

		elog_complain( 0, "units_convert result: '%s' did not match "
			     "request for '%s'\n",
			     mp.cgg_pga->units, DEFAULT_MMI_PGA_UNITS );

	} else if( rc == -1 ) {

		elog_complain( 0, "units_convert result: '%s' not recognized\n",
			     mp.cgg_pga->units );
	}

	cggrid_apply( &mp.cgg_pga, cggridval_mult_static, &conversion_factor  );

	fp = fopen( pgv_filename, "r" );
	cgg_pgv_orig = cggrid_read( fp );	
	fclose( fp );

	if( cgg_pgv_orig == (CGGrid *) NULL ) {

		elog_complain( 0,
			  "trinetsm_es99: Error reading pgv file "
			  "%s.\n", pgv_filename );
		cggrid_free( &mp.cgg_pga );

		/* Strictly speaking it's would only be necessary to 
		   fail here if the MMI is going to be greater than
		   the mp.low_cutoff_mmi value */

		return -1;
	}

	/* Create a PGV grid co-registered with the PGA grid: */

	mp.cgg_pgv = cggrid_new( mp.cgg_pga->minx, 
			  	 mp.cgg_pga->maxx, 
			  	 mp.cgg_pga->miny, 
			  	 mp.cgg_pga->maxy, 
				 mp.cgg_pga->dx, 
				 mp.cgg_pga->dy );

	cggrid_reregister( cgg_pgv_orig, mp.cgg_pgv );
	cggrid_free( &cgg_pgv_orig );

	rc = units_convert( 1.0, mp.cgg_pgv->units, DEFAULT_MMI_PGV_UNITS, 
			    &conversion_factor, returned_units );

	if( rc == 1 ) {

		elog_complain( 0, "units_convert result: '%s' did not match "
			     "request for '%s'\n",
			     mp.cgg_pgv->units, DEFAULT_MMI_PGV_UNITS );

	} else if( rc == -1 ) {

		elog_complain( 0, "units_convert result: '%s' not recognized\n",
			     mp.cgg_pgv->units );
	}

	cggrid_apply( &mp.cgg_pgv, cggridval_mult_static, &conversion_factor  );

	cgg = cggrid_new( mp.cgg_pga->minx, 
			  mp.cgg_pga->maxx, 
			  mp.cgg_pga->miny, 
			  mp.cgg_pga->maxy, 
			  mp.cgg_pga->dx, 
			  mp.cgg_pga->dy );

	cggrid_apply( &cgg, compute_mmi, &mp );

	strcpy( cgg->units, output_units );

	if( Force ) {

		flags |= CG_OVERWRITE;
	}
	
	rc = cggrid2db( db, cgg, recipe_name, grid_name, 
		        output_file, qgridfmt, output_units, 
		        qgridtype, auth, flags );

	if( rc < 0 ) {
		elog_clear_register( 1 );
	}

	cggrid_free( &mp.cgg_pga );
	cggrid_free( &mp.cgg_pgv );
	cggrid_free( &cgg );

	return 0;
}

int
trinetsm_es99( Dbptr db, Pf *pf )

	/* Wald et al. Earthquake Spectra 15, pp 537ff. (1999) */
{
	CGGrid	*cgg;
	CGGrid	*sitegrid;
	CGGrid	*sitecorr_orig;
	CGPoint *cgpt;
	CGPointset *cgps;
	CGPointset *phantom;
	Dbptr	dbg;
	Dbptr	dbbundle;
	double	latc;
	double	lonc;
	double	alat;
	double	alon;
	double	wdellon;
	double	edellon;
	double	sdellat;
	double	ndellat;
	double	qdlat;
	double	qdlon;
	double	gridval;
	double	ml;
	double	mag;
	double	minx;
	double	maxx;
	double	miny;
	double	maxy;
	double 	vs30_default = 0;
	double	phantom_spacing_deg;
	double	phantom_mindist_deg;
	double	centroid_mindist_deg;
	double	const_coeff;
	double 	mag_coeff;
	double	r_coeff;
	double	P_coeff;
	double	S_coeff;
	double	S_hardrock_cutoff_mps;
	double	P;
	double	S = 0;
	double	r_offset;
	double	d;
	double	az;
	double	tension;
	double	overrelaxation;
	double	convergence; 
	double	aval;
	double	aval_converted;
	double	pva_val;
	double	pva_val_converted;
	double	pvv_val;
	double	pvv_val_converted;
	char	*output_units;
	char	aunit[STRSZ];
	char	expr[STRSZ];
	char	grid_name[STRSZ];
	char 	*recipe_name;
	char	*output_file;
	char	*qgridfmt;
	char	*qgridtype;
	char	*auth;
	char	*sitecorr_file = 0;
	int	orid;
	int	max_iterations;
	int	ilat;
	int	ilon;
	int	ipt;
	int	start;
	int	end;
	int	flags = 0;
	int	nrecs;
	int	nchannels;
	Tbl	*s;
	Tbl	*views;
	Tbl	*view_tables;
	int	ns, ne;
	FILE	*fp;
	struct mypdata_ *mp;
	enum Wfmode wfmode;
	void	*private = NULL;
	char	*wfmgme = "wfmgme";
	double	sitevel;
	double	ampfactor;
	double	conversion_factor;
	enum Calcmode { PGA, PGV } calcmode;
	char	returned_units[STRSZ];
	char	input_info[STRSZ];
	int	rc;

	qgridtype = pfget_string( pf, "qgridtype" );

	if( ! strcmp( qgridtype, "pga" ) ) {

		calcmode = PGA;

	} else if( ! strcmp( qgridtype, "pgv" ) ) {

		calcmode = PGV;

	} else {

		elog_complain( 0, 
			  "trinetsm_es99: unknown qgridtype '%s' " 
			  "requested in parameter file, Can't continue\n",
			  qgridtype );
		return -1;
	}

	if( Verbose ) {
		elog_notify( 0, "trinetsm_es99: Calculating a " 
				"%s grid\n", qgridtype );
	}

	recipe_name = pfget_string( pf, "recipe_name" );
	output_file = pfget_string( pf, "output_file" );
	qgridfmt = pfget_string( pf, "qgridfmt" );
	auth = pfget_string( pf, "auth" );

	views = newtbl( 0 );

	dbquery( db, dbVIEW_TABLES, &view_tables );
	view_tables = duptbl( view_tables, (void *(*)()) strdup );
	sorttbl( view_tables, (int (*)(void *, void *, void *)) cmp_string, 0 );
	searchtbl( (char *) &wfmgme, view_tables, 
		   (int (*)(void *, void *, void *)) cmp_string,
		   private, &ns, &ne );
	if( ns > ne ) { 
		wfmode = WFMEAS;
	} else {
		wfmode = WFMGME;
	}

	if( Verbose ) {

		elog_notify( 0, "trinetsm_es99: wfmode is %s\n", 
			 wfmode == WFMEAS ? "wfmeas" : "wfmgme" );
	}

	if( wfmode == WFMEAS && calcmode == PGV ) {

		elog_complain( 0, "trinetsm_es99: pgv calculation not "
				  "currently supported for wfmeas rows "
				  "(use wfmgme rows)\n" );
		free_views( db, views );
		return -1;
	}

	freetbl( view_tables, free );

	wdellon = pfget_double( pf, "wdellon" );
	edellon = pfget_double( pf, "edellon" );
	sdellat = pfget_double( pf, "sdellat" );
	ndellat = pfget_double( pf, "ndellat" );
	qdlat = pfget_double( pf, "qdlat" );
	qdlon = pfget_double( pf, "qdlon" );

	sitecorr_file = pfget_string( pf, "sitecorr_file" );
	vs30_default = pfget_double( pf, "vs30_default_mps" );

	phantom_spacing_deg = pfget_double( pf, "phantom_spacing_deg" );
	phantom_mindist_deg = pfget_double( pf, "phantom_mindist_deg" );
	centroid_mindist_deg = pfget_double( pf, "centroid_mindist_deg" );
	const_coeff = pfget_double( pf, "const_coeff" );
	mag_coeff = pfget_double( pf, "mag_coeff" );
	r_coeff = pfget_double( pf, "r_coeff" );
	P_coeff = pfget_double( pf, "P_coeff" );
	S_coeff = pfget_double( pf, "S_coeff" );
	S_hardrock_cutoff_mps = pfget_double( pf, "S_hardrock_cutoff_mps" );
	P = pfget_double( pf, "P" );
	r_offset = pfget_double( pf, "r_offset" );
	tension = pfget_double( pf , "tension" );
	overrelaxation = pfget_double( pf, "overrelaxation" );
	convergence = pfget_double( pf, "convergence" );
	max_iterations = pfget_int( pf, "max_iterations" );

	if( setup_ampfactors( pf ) < 0 ) {
		elog_complain( 0, 
			  "trinetsm_es99: failed to setup " 
			  "amplification factors, Can't continue\n" );
		free_views( db, views );
		return -1;
	}

	output_units = pfget_string( pf, "output_units" );

	if( output_units == NULL ) {

		allot( char *, output_units, 20 );

		if( calcmode == PGA ) {
			strcpy( output_units, DEFAULT_PGA_OUTPUT_UNITS );
		} else {
			strcpy( output_units, DEFAULT_PGV_OUTPUT_UNITS );
		}

		elog_complain( 0, "trinetsm_es99: No output_units specified"
			     " in parameter file! Defaulting to %s\n",
			     output_units );
	}

	db.record = 0;
	dbgetv( db, 0, "origin.lat", &latc, 
		       "origin.lon", &lonc,
		       "ml", &ml, 
		       "orid", &orid,
		       0 );

	sprintf( grid_name, "orid_%d", orid );

	if( ml == -999.00 ) {

		elog_complain( 0, 
			  "trinetsm_es99: Null local magnitude for " 
			  "orid %d, Can't continue\n", 
			  orid );
		free_views( db, views );
		free_ampfactors();
		return -1;

	} else {

		mag = ml;
	}

	minx = lonc + wdellon;
	maxx = lonc + edellon;
	miny = latc + sdellat;
	maxy = latc + ndellat;

	cgg = cggrid_new( minx, maxx, miny, maxy, qdlon, qdlat );
	sitegrid = cggrid_new( minx, maxx, miny, maxy, qdlon, qdlat );

	if( calcmode == PGA ) {
		strcpy( cgg->units, DEFAULT_PGA_CALC_UNITS );
	} else {
		strcpy( cgg->units, DEFAULT_PGV_CALC_UNITS );
	}
	
	if( Verbose ) {

		elog_notify( 0, "trinetsm_es99: Reading "
				 "site-corrections grid\n" );
	}

	if( sitecorr_file == 0 || ( ! strcmp( sitecorr_file, "" ) ) ) {

		if( Verbose ) {

			elog_notify( 0, 
				  "trinetsm_es99: No sitecorr file: "
				  "setting Vs30 to specified "
				  "default %.2f m/s\n", vs30_default );
		}

	} else {

		if( ( fp = fopen( sitecorr_file, "r" ) ) != NULL ) {

			sitecorr_orig = cggrid_read( fp );	
			
			if( sitecorr_orig == (CGGrid *) NULL ) {

				elog_complain( 0,
				  	"trinetsm_es99: Error reading sitecorr file "
				  	"%s. Setting Vs30 to specified "
				  	"default %.2f m/s\n",
				  	sitecorr_file, vs30_default );

				cggrid_apply( &sitegrid, cggridval_set, 
				      	&vs30_default );

			} else {

				cggrid_reregister( sitecorr_orig, sitegrid );
				cggrid_free( &sitecorr_orig );
			}

			fclose( fp );

		} else {

			elog_complain( 1,
			  	"trinetsm_es99: Cannot find sitecorr file "
			  	"%s. Setting Vs30 to specified "
			  	"default %.2f m/s\n",
			  	sitecorr_file, vs30_default );

			cggrid_apply( &sitegrid, cggridval_set, 
			      	&vs30_default );
		}
	}

	cggrid_apply( &sitegrid, cggridval_resetnan, &vs30_default );

	cgps = cgpointset_new();

	s = strtbl( "sta", "arrival.time", 0 );
	dbg = dbsort( db, s, 0, 0 );
	pushtbl( views, (char *) dbg.table );
	dbg = dbgroup( dbg, s, 0, 1 );
	pushtbl( views, (char *) dbg.table );
	freetbl( s, 0 );
	dbquery( dbg, dbRECORD_COUNT, &nrecs );

	for( dbg.record = 0; dbg.record < nrecs; dbg.record++ ) {

		dbgetv( dbg, 0, "bundle", &dbbundle, 0 );

		if( wfmode == WFMEAS ) {
			dbex_evalstr( dbg, "count()", dbINTEGER, &nchannels );
		
			if( nchannels > 3 ) {
				elog_complain( 0,
				   "not expecting more than 3 channels "
				   "per station, skipping sta %s\n", mp->sta );
				continue;
			}
		}

		dbget_range( dbbundle, &start, &end );

		allot( struct mypdata_ *, mp, 1 );

		db.record = start;
		dbgetv( db, 0, "site.lon", &alon, "site.lat", &alat, 
				"sta", mp->sta, 0 );

		if( wfmode == WFMEAS ) {

			pva_val = 0;
			for( db.record = start; db.record < end; db.record++ ) {
	
				dbgetv( db, 0, "val1", &aval, 
				       	"units1", aunit, 0 );

				sprintf( input_info, "%.8f %s", aval, aunit );
						
				rc = units_convert( 1.0, input_info,
						    DEFAULT_PGA_CALC_UNITS, 
						    &aval_converted, 
						    returned_units ); 
				aval = aval_converted;
				if( rc == 1 ) {

					elog_complain( 0, 
						"trinetsm_es99: '%s' "
						"did not match "
			     			"request for '%s'\n",
			     			DEFAULT_PGA_CALC_UNITS, 
						output_units );

				} else if( rc == -1 ) {

					elog_complain( 0,
						"trinetsm_es99: '%s' not "
						"recognized\n",
			     			DEFAULT_PGA_CALC_UNITS );
				}

				pva_val += aval * aval;
			}

			pva_val = sqrt( pva_val );

		} else {

			/* wfmode WFMGME: expecting only one row */	

			db.record = start;
			dbgetv( db, 0, "pva", &pva_val, 0 );
			dbgetv( db, 0, "pvv", &pvv_val, 0 );

			sprintf( input_info, "%.8f %s", 
					pva_val, DEFAULT_WFMGME_PVA_UNITS );
						
			rc = units_convert( 1.0, input_info,
				       DEFAULT_PGA_CALC_UNITS, 
				       &pva_val_converted,
				       returned_units );
			pva_val = pva_val_converted;
			if( rc == 1 ) {

				elog_complain( 0, 
					"trinetsm_es99: '%s' did not match "
			     		"request for '%s'\n",
			     		DEFAULT_PGA_CALC_UNITS, output_units );

			} else if( rc == -1 ) {

				elog_complain( 0,
					"trinetsm_es99: '%s' not recognized\n",
			     		DEFAULT_PGA_CALC_UNITS );
			}

			sprintf( input_info, "%.8f %s", 
					pvv_val, DEFAULT_WFMGME_PVV_UNITS );
						
			rc = units_convert( 1.0, input_info,
				       DEFAULT_PGV_CALC_UNITS, 
				       &pvv_val_converted,
				       returned_units );
			pvv_val = pvv_val_converted;
			if( rc == 1 ) {

				elog_complain( 0, 
					"trinetsm_es99: '%s' did not match "
			     		"request for '%s'\n",
			     		DEFAULT_PGV_CALC_UNITS, output_units );

			} else if( rc == -1 ) {

				elog_complain( 0,
					"trinetsm_es99: '%s' not recognized\n",
			     		DEFAULT_PGV_CALC_UNITS );
			}
		}

		if( calcmode == PGA ) {

			elog_notify( 0,
				"trinetsm_es99: Station %s at %.3f,%.3f "
				"%.3g %s\n",
				mp->sta, alon, alat, pva_val, 
				DEFAULT_PGA_CALC_UNITS );
		} else {

			elog_notify( 0,
				"trinetsm_es99: Station %s at %.3f,%.3f "
				"%.3g %s\n",
				mp->sta, alon, alat, pvv_val, 
				DEFAULT_PGV_CALC_UNITS );
		}

		sitevel = cggrid_probe( sitegrid, alon, alat );

		if( is_nan( sitevel ) ) {

			elog_notify( 0, 
			   "trinetsm_es99: Station %s outside site-correction "
			   "grid; not correcting to rock value\n", mp->sta );

		} else {

			/* Acceleration amp factors are used to correct
			   both PGA and PGV values (albeit in different 
			   frequency bands, hopefully set up correctly
			   in parameter file): */

			if( calcmode == PGA ) {

				ampfactor = pga_ampfactor( sitevel, pva_val );

			} else {
				/* Ignore nonlinear response for 
				   PGV calculation: */

				ampfactor = pga_ampfactor( sitevel, 0 );
			}

			if( ! is_nan( ampfactor ) && ampfactor != 0 ) {

				elog_notify( 0, 
			   	   "trinetsm_es99: Correcting %s to rock value "
				   "with site amplification factor %.3f\n",
				   mp->sta, ampfactor );

				if( calcmode == PGA ) {

					pva_val /= ampfactor;

				} else {

					pvv_val /= ampfactor;
				}

			} else {
				elog_notify( 0, 
			   	   "trinetsm_es99: No site-amplification factor for %s; "
			   	   "not correcting to rock value\n", mp->sta );
			}
		}

		if( calcmode == PGA ) {
			cgpt = cgpoint_init( alon, alat, 0, 
					     pva_val, 0, dup_pdata, free );
		} else {
			cgpt = cgpoint_init( alon, alat, 0, 
					     pvv_val, 0, dup_pdata, free );
		}

		cgpt->pdata = (void *) mp;

		cgpointset_addpoint( cgps, cgpt );
		cgpoint_free( &cgpt );
	}

	/* Add phantom stations */

	phantom = cgpointset_new();

	if( nearest_deg( cgps, lonc, latc ) > centroid_mindist_deg ) { 

		sitevel = cggrid_probe( sitegrid, lonc, latc );

		if( calcmode == PGV && sitevel < S_hardrock_cutoff_mps ) {
			S = 1;
		} else {
			S = 0;
		}

		cgpt = new_phantom( lonc, latc, 0, mag, const_coeff, 
				    r_offset, r_coeff, mag_coeff, 
				    P_coeff, P, S_coeff, S );

		cgpointset_addpoint( phantom, cgpt );
		cgpoint_free( &cgpt );
	}

	for( alon = cgg->minx, ilon = 0;
	     alon <= cgg->maxx; 
	     alon = CGGRID_SAMP2COORD(cgg->minx,phantom_spacing_deg,++ilon )) {

	  for( alat = cgg->miny, ilat = 0;
	       alat <= cgg->maxy; 
	       alat = CGGRID_SAMP2COORD(cgg->miny,phantom_spacing_deg,++ilat )) {

		if( nearest_deg( cgps, alon, alat ) < phantom_mindist_deg ) {

			continue;
		}

		db.record = 0;
		sprintf( expr, 
			 "distance( lat, lon, %f, %f ) * 111.195", 
			 alat, alon );
		dbex_evalstr( db, expr, dbREAL, (void *) &d );

		sitevel = cggrid_probe( sitegrid, alon, alat );

		if( calcmode == PGV && sitevel < S_hardrock_cutoff_mps ) {
			S = 1;
		} else {
			S = 0;
		}

		cgpt = new_phantom( alon, alat, d, mag, const_coeff, 
				    r_offset, r_coeff, mag_coeff, 
				    P_coeff, P, S_coeff, S );

		cgpointset_addpoint( phantom, cgpt );
		cgpoint_free( &cgpt );
	  }
	}

	for( ipt = 0; ipt < cgpointset_cnt( phantom ); ipt++ ) {

		cgpt = cgpointset_getpoint( phantom, ipt );
		cgpointset_addpoint( cgps, cgpt );
	}

	cgpointset_free( &phantom );

	if( Verbose ) {
		
		elog_notify( 0, "trinetsm_es99: fitting measurements\n" );
	}

	cggrid_splinefit_pointset( cgg, cgps, tension, tension, 1,
				   overrelaxation, convergence, 
				   max_iterations, 0, a_quiet_nan() );

	cggrid_apply( &cgg, apply_ampfactor, sitegrid );

	rc = units_convert( 1.0, cgg->units, output_units, 
			    &conversion_factor, returned_units );

	if( rc == 1 ) {

		elog_complain( 0, "units_convert result: '%s' did not match "
			     "request for '%s'\n",
			     cgg->units, output_units );

	} else if( rc == -1 ) {

		elog_complain( 0, "units_convert result: '%s' not recognized\n",
			     cgg->units );
	}

	cggrid_apply( &cgg, cggridval_mult_static, &conversion_factor  );

	strcpy( cgg->units, output_units );

	if( Force ) {

		flags |= CG_OVERWRITE;
	}
	
	rc = cggrid2db( db, cgg, recipe_name, grid_name, 
		        output_file, qgridfmt, output_units, 
		        qgridtype, auth, flags );
	
	if( rc < 0 ) {

		elog_clear_register( 1 );
	}

	cggrid_free( &cgg );
	cggrid_free( &sitegrid );
	cgpointset_free( &cgps );
	free_ampfactors();
	free_views( db, views );

	return 0;
}
