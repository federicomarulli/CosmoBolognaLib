/*------------------------------------------------------------------------------
© A J S Hamilton 2003
------------------------------------------------------------------------------*/
#ifndef MANGLEFN_H
#define MANGLEFN_H

#include "defines.h"
#include "format.h"
#include "harmonics.h"
#include "logical.h"
#include "polygon.h"
#include "vertices.h"
#include "polysort.h"

void	advise_fmt(format *);

void	azel_(long double *, long double *, long double *, long double *, long double *, long double *, long double *);
void	azell_(long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *);

void	braktop(long double, int *, long double [], int, int);
void	brakbot(long double, int *, long double [], int, int);
void	braktpa(long double, int *, long double [], int, int);
void	brakbta(long double, int *, long double [], int, int);
void	braktop_(long double *, int *, long double [], int *, int *);
void	brakbot_(long double *, int *, long double [], int *, int *);
void	braktpa_(long double *, int *, long double [], int *, int *);
void	brakbta_(long double *, int *, long double [], int *, int *);

void	cmminf(polygon *, int *, long double *);

void	vert_to_poly(vertices *, polygon *);
void	edge_to_poly(vertices *, int, int *, polygon *);
void    rect_to_poly(long double [4], polygon *);

#ifdef	GCC
void	rps_to_vert(int nv, vec [nv], vertices *);
#else
void	rps_to_vert(int nv, vec [/*nv*/], vertices *);
#endif
void	rp_to_azel(vec, azel *);
void	azel_to_rp(azel *, vec);
void	azel_to_gc(azel *, azel *, vec, long double *);
void	rp_to_gc(vec, vec, vec, long double *);
void	edge_to_rpcm(azel *, azel *, azel *, vec, long double *);
void	rp_to_rpcm(vec, vec, vec, vec, long double *);
void	circ_to_rpcm(long double [3], vec, long double *);
void	rpcm_to_circ(vec, long double *, long double [3]);
void	az_to_rpcm(long double, int, vec, long double *);
void	el_to_rpcm(long double, int, vec, long double *);
long double	thij(vec, vec);
long double	cmij(vec, vec);
int	poly_to_rect(polygon *, long double *, long double *, long double *, long double *);
int	antivert(vertices *, polygon *);

void	copy_format(format *, format *);

void	copy_poly(polygon *, polygon *);
void	copy_polyn(polygon *, int, polygon *);
void	poly_poly(polygon *, polygon *, polygon *);
void	poly_polyn(polygon *, polygon *, int, int, polygon *);
#ifdef	GCC
void	group_poly(polygon *poly, int [poly->np], int, polygon *);
#else
void	group_poly(polygon *poly, int [/*poly->np*/], int, polygon *);
#endif

void    assign_parameters();
void    pix2ang(int, unsigned long, long double *, long double *);
void    ang2pix(int, long double, long double, unsigned long *);
void    pix2ang_radec(int, unsigned long, long double *, long double *);
void    ang2pix_radec(int, long double, long double, unsigned long *);
void    csurvey2eq(long double, long double, long double *, long double *);
void    eq2csurvey(long double, long double, long double *, long double *);
void    superpix(int, unsigned long, int, unsigned long *);
void    subpix(int, unsigned long, unsigned long *, unsigned long *, unsigned long *, unsigned long *);
void    pix_bound(int, unsigned long, long double *, long double *, long double *, long double *);
long double  pix_area(int, unsigned long);
void    pix2xyz(int, unsigned long, long double *, long double *, long double *);
void    area_index(int, long double, long double, long double, long double, unsigned long *, unsigned long *, unsigned long *, unsigned long *);
void    area_index_stripe(int, int, unsigned long *, unsigned long *, unsigned long *, unsigned long *);

long double	drandom(void);

#ifdef	GCC
int	cmlim_polys(int npoly, polygon *[npoly], long double, vec);
int	drangle_polys(int npoly, polygon *[npoly], long double, vec, int nth, long double [nth], long double [nth]);
#else
int	cmlim_polys(int npoly, polygon *[/*npoly*/], long double, vec);
int	drangle_polys(int npoly, polygon *[/*npoly*/], long double, vec, int nth, long double [/*nth*/], long double [/*nth*/]);
#endif

void	cmlimpolys_(long double *, vec);
#ifdef	GCC
void	dranglepolys_(long double *, vec, int *nth, long double [*nth], long double [*nth]);
#else
void	dranglepolys_(long double *, vec, int *nth, long double [/**nth*/], long double [/**nth*/]);
#endif

#ifdef	GCC
void	dump_poly(int npoly, polygon *[npoly]);
#else
void	dump_poly(int, polygon *[/*npoly*/]);
#endif

void	fframe_(int *, long double *, long double *, int *, long double *, long double *);

void	findtop(long double [], int, int [], int);
void	findbot(long double [], int, int [], int);
void	findtpa(long double [], int, int [], int);
void	findbta(long double [], int, int [], int);
void	finitop(int [], int, int [], int);
void	finibot(int [], int, int [], int);
void	finitpa(int [], int, int [], int);
void	finibta(int [], int, int [], int);

void	findtop_(long double [], int *, int [], int *);
void	findbot_(long double [], int *, int [], int *);
void	findtpa_(long double [], int *, int [], int *);
void	findbta_(long double [], int *, int [], int *);
void	finitop_(int [], int *, int [], int *);
void	finibot_(int [], int *, int [], int *);
void	finitpa_(int [], int *, int [], int *);
void	finibta_(int [], int *, int [], int *);

polygon *get_pixel(int,char);
int     get_child_pixels(int, int [], char);
int     get_parent_pixels(int, int [], char);
int     get_res(int,char);

void    healpix_ang2pix_nest(int, long double, long double, int *);
polygon *get_healpix_poly(int, int);
int     get_nside(int);
void    healpix_verts(int, int, vec, long double []);
void    pix2vec_nest__(int *, int *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *);
long double  cmrpirpj(vec, vec);

int	garea(polygon *, long double *, int, long double *);
int	gcmlim(polygon *, long double *, vec, long double *, long double *);
int	gphbv(polygon *, int, int, long double *, long double [2], long double [2]);
int	gphi(polygon *, long double *, vec, long double, long double *);
int	gptin(polygon *, vec);
#ifdef	GCC
int	gspher(polygon *, int lmax, long double *, long double *, long double [2], long double [2], harmonic [NW]);
int	gsphera(long double, long double, long double, long double, int lmax, long double *, long double [2], long double [2], harmonic [NW]);
int	gsphr(polygon *, int lmax, long double *, harmonic [NW]);
int	gsphra(long double, long double, long double, long double, int lmax, harmonic [NW]);
#else
int	gspher(polygon *, int lmax, long double *, long double *, long double [2], long double [2], harmonic [/*NW*/]);
int	gsphera(long double, long double, long double, long double, int lmax, long double *, long double [2], long double [2], harmonic [/*NW*/]);
int	gsphr(polygon *, int lmax, long double *, harmonic [/*NW*/]);
int	gsphra(long double, long double, long double, long double, int lmax, harmonic [/*NW*/]);
#endif
int	gverts(polygon *, int, long double *, int, int, int *, vec **, long double **, int **, int **, int *, int *, int **);
#ifdef	GCC
int	gvert(polygon *poly, int, long double *, int nvmax, int, int nve, int *, vec [nvmax * nve], long double [nvmax], int [nvmax], int [poly->np], int *, int *, int [nvmax]);
#else
int	gvert(polygon *poly, int, long double *, int nvmax, int, int nve, int *, vec [/*nvmax * nve*/], long double [/*nvmax*/], int [/*nvmax*/], int [/*poly->np*/], int *, int *, int [/*nvmax*/]);
#endif
int	gvlims(polygon *, int, long double *, vec, int *, vec **, vec **, long double **, long double **, long double **, long double **, int **, int **, int *, int *, int **);
#ifdef	GCC
int	gvlim(polygon *poly, int, long double *, vec, int nvmax, int *, vec [nvmax], vec [nvmax], long double [nvmax], long double [nvmax], long double [poly->np], long double [poly->np], int [nvmax], int [poly->np], int *, int *, int [nvmax]);
#else
int	gvlim(polygon *poly, int, long double *, vec, int nvmax, int *, vec [/*nvmax*/], vec [/*nvmax*/], long double [/*nvmax*/], long double [/*nvmax*/], long double [/*poly->np*/], long double [/*poly->np*/], int [/*nvmax*/], int [/*poly->np*/], int *, int *, int [/*nvmax*/]);
#endif
int	gvphi(polygon *, vec, long double, vec, long double *, long double *, vec);

void	garea_(long double *, vec [], long double [], int *, long double *, int *, long double *, int *, logical *);
void	gaxisi_(vec, vec, vec);
void	gcmlim_(vec [], long double [], int *, vec, long double *, long double *, long double *, long double *, int *);
void	gphbv_(long double [2], long double [2], vec [], long double [], int *, int *, int *, int *, long double *, long double *, int *);
void	gphi_(long double *, vec [], long double [], int *, vec, long double *, long double *, long double *, int *);
logical	gptin_(vec [], long double [], int *, vec);
void	gspher_(long double *, long double [2], long double [2], harmonic [], int *, int *, int *, vec [], long double [], int *, int *, int *, int *, long double *, long double *, int *, long double *, logical *);
void	gsphera_(long double *, long double [2], long double [2], harmonic [], int *, int *, int *, int *, long double *, long double *, long double *, long double *, long double *, long double *);
void	gvert_(vec [], long double [], int [], int [], int [], int *, int *, int *, int *, int *, int *, vec [], long double [], int *, int *, long double *, long double *, int *, long double *, int *, logical *);

void	gvlim_(vec [], vec [], long double [], long double [], long double [], long double [], int [], int [], int [], int *, int *, int *, int *, vec [], long double [], int *, long double [], int *, long double *, long double *, int *, long double *, int *, logical *);
void	gvphi_(long double *, vec, vec [], long double [], int *, vec, long double *, vec, long double *, long double *, int *);

#ifdef	GCC
int	harmonize_polys(int npoly, polygon *[npoly], long double, int lmax, harmonic w[NW]);
#else
int	harmonize_polys(int npoly, polygon *poly[/*npoly*/], long double, int lmax, harmonic w[/*NW*/]);
#endif

void	harmonizepolys_(long double *, int *, harmonic []);

void	ikrand_(int *, long double *);
void	ikrandp_(long double *, long double *);
void	ikrandm_(long double *, long double *);

void	msg(char *, ...);

polygon	*new_poly(int);
void	free_poly(polygon *);
int	room_poly(polygon **, int, int, int);
void	memmsg(void);

vertices	*new_vert(int);
void	free_vert(vertices *);

void	parse_args(int, char *[]);

int	parse_fopt(void);

#ifdef	GCC
int	partition_poly(polygon **, int npolys, polygon *[npolys], long double, int, int, int, int, int *);
int	partition_gpoly(polygon *, int npolys, polygon *[npolys], long double, int, int, int, int *);
int	part_poly(polygon *, int npolys, polygon *[npolys], long double, int, int, int, int *, int *);
int     pixel_list(int npoly, polygon *[npoly], int max_pixel, int [max_pixel], int [max_pixel]);
#else
int	partition_poly(polygon **, int npolys, polygon *[/*npolys*/], long double, int, int, int, int, int *);
int	partition_gpoly(polygon *, int npolys, polygon *[/*npolys*/], long double, int, int, int, int *);
int	part_poly(polygon *, int npolys, polygon *[/*npolys*/], long double, int, int, int, int *, int *);
int     pixel_list(int npoly, polygon *[/*npoly*/], int max_pixel, int [/*max_pixel*/], int [/*max_pixel*/]);
#endif
int     pixel_start(int, char);

long double	places(long double, int);
int     poly_cmp(polygon **, polygon **);

#ifdef	GCC
int	poly_id(int npoly, polygon *[npoly], long double, long double, int **, long double **);
void	poly_sort(int npoly, polygon *[npoly], char);
#else
int	poly_id(int npoly, polygon *[/*npoly*/], long double, long double, int **);
void	poly_sort(int npoly, polygon *[/*npoly*/], char);
#endif


int	prune_poly(polygon *, long double);
int	trim_poly(polygon *);
int	touch_poly(polygon *);

int	rdangle(char *, char **, char, long double *);

#ifdef	GCC
int	rdmask(char *, format *, int npolys, polygon *[npolys]);
#else
int	rdmask(char *, format *, int npolys, polygon *[/*npolys*/]);
#endif

void	rdmask_(void);

int	rdspher(char *, int *, harmonic **);

void	scale(long double *, char, char);
void	scale_azel(azel *, char, char);
void	scale_vert(vertices *, char, char);

#ifdef	GCC
int	search(int n, long double [n], long double);
#else
int	search(int n, long double [/*n*/], long double);
#endif

#ifdef	GCC
int	snap_polys(format *fmt, int npoly, polygon *poly[npoly], int, long double, long double, long double, long double, long double, int, char *);
#else
int	snap_polys(format *fmt, int npoly, polygon *poly[/*npoly*/], int, long double, long double, long double, long double, long double, int, char *);
#endif
int	snap_poly(polygon *, polygon *, long double, long double);
int	snap_polyth(polygon *, polygon *, long double, long double, long double);

int	split_poly(polygon **, polygon *, polygon **, long double, char);
#ifdef	GCC
int	fragment_poly(polygon **, polygon *, int, int npolys, polygon *[npolys], long double, char);
#else
int	fragment_poly(polygon **, polygon *, int, int npolys, polygon *[/*npolys*/], long double, char);
#endif

int	strcmpl(const char *, const char *);
int	strncmpl(const char *, const char *, size_t);

int	strdict(char *, char *[]);
int	strdictl(char *, char *[]);

#ifdef	GCC
int	vmid(polygon *, long double, int nv, int nve, vec [nv * nve], int [nv], int [nv], int *, vec **);
int	vmidc(polygon *, int nv, int nve, vec [nv * nve], int [nv], int [nv], int *, vec **);
#else
int	vmid(polygon *, long double, int nv, int nve, vec [/*nv * nve*/], int [/*nv*/], int [/*nv*/], int *, vec **);
int	vmidc(polygon *, int nv, int nve, vec [/*nv * nve*/], int [/*nv*/], int [/*nv*/], int *, vec **);
#endif

long double	weight_fn(long double, long double, char *);
long double	rdweight(char *);

long double	twoqz_(long double *, long double *, int *);
long double	twodf100k_(long double *, long double *);
long double	twodf230k_(long double *, long double *);

int     which_pixel(long double, long double, int, char);

#ifdef	GCC
void	wrangle(long double, char, int, size_t str_len, char [str_len]);
#else
void	wrangle(long double, char, int, size_t str_len, char [/*str_len*/]);
#endif

#ifdef	GCC
long double	wrho(long double, long double, int lmax, int, harmonic w[NW], long double, long double);
#else
long double	wrho(long double, long double, int lmax, int, harmonic w[/*NW*/], long double, long double);
#endif

long double	wrho_(long double *, long double *, harmonic *, int *, int *, int *, int *, long double *, long double *);

#ifdef	GCC
int	wrmask(char *, format *, int npolys, polygon *[npolys]);
int	wr_circ(char *, format *, int npolys, polygon *[npolys], int);
int	wr_edge(char *, format *, int npolys, polygon *[npolys], int);
int	wr_rect(char *, format *, int npolys, polygon *[npolys], int);
int	wr_poly(char *, format *, int npolys, polygon *[npolys], int);
int	wr_Reg(char *, format *, int npolys, polygon *[npolys], int);
int	wr_area(char *, format *, int npolys, polygon *[npolys], int);
int	wr_id(char *, int npolys, polygon *[npolys], int);
int	wr_midpoint(char *, format *, int npolys, polygon *[npolys], int);
int	wr_weight(char *, format *, int npolys, polygon *[npolys], int);
int     wr_healpix_weight(char *, format *, int numweight, long double [numweight]);
int	wr_list(char *, format *, int npolys, polygon *[npolys], int);
int	discard_poly(int npolys, polygon *[npolys]);
#else
int	wrmask(char *, format *, int npolys, polygon *[/*npolys*/]);
int	wr_circ(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_edge(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_rect(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_poly(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_Reg(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_area(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_id(char *, int npolys, polygon *[/*npolys*/], int);
int	wr_midpoint(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	wr_weight(char *, format *, int npolys, polygon *[/*npolys*/], int);
int     wr_healpix_weight(char *, format *, int numweight, long double [/*numweight*/]);
int	wr_list(char *, format *, int npolys, polygon *[/*npolys*/], int);
int	discard_poly(int npolys, polygon *[/*npolys*/]);
#endif

int	wrrrcoeffs(char *, long double, long double [2], long double [2]);

#ifdef	GCC
int	wrspher(char *, int lmax, harmonic [NW]);
#else
int	wrspher(char *, int lmax, harmonic [/*NW*/]);
#endif

#endif	/* MANGLEFN_H */
