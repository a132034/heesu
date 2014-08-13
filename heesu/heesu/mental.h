#pragma once
#include	<cmath>
#include	<memory>
#include	<fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <vector>
#include <time.h>
#include <string>
#include <memory.h>
#include <algorithm>
#include <functional>      // For greater<int>()
#include <iostream>
#define QX_DEF_IS_PPM 1
#define QX_DEF_IS_PGM 0

#define QX_DEF_MST_KI_MIN_DISTANT			1
#define QX_DEF_MST_KI_4NR_NEIGHBOR			4
#define QX_DEF_MST_KI_8NR_NEIGHBOR			8
#define QX_DEF_MST_KI_PARENT_DEFAULT		-1
#define QX_DEF_MST_KI_SIGMA_RANGE			0.1

#define CONST const

typedef unsigned long       DWORD;
typedef int                 BOOL;
typedef double              DOUBLE;
typedef unsigned char       BYTE, byte;
typedef unsigned short      WORD;
typedef float               FLOAT;
typedef void                VOID;
typedef char                CHAR;
typedef long				LONG;
typedef const char*         LPCSTR;


typedef int                 INT;
typedef unsigned int        UINT;
typedef unsigned int        *PUINT;
using namespace std;  // std c++ libs implemented in std

#define QX_DEF_PADDING					10
#define QX_DEF_THRESHOLD_ZERO			1e-6
#define QX_DEF_PI_DOUBLE				3.14159265359
#define QX_DEF_FLOAT_MAX				1.175494351e+38F
#define QX_DEF_DOUBLE_MAX				1.7E+308
#define QX_DEF_FLOAT_RELATIVE_ACCURACY	2.2204e-016
#define QX_DEF_INI_MAX					2147483647
#define QX_DEF_SHORT_MAX				65535
#define QX_DEF_CHAR_MAX					255
#define	QX_DEF_SEED						42
#define QX_DEF_THRESHOLD_ZERO			1e-6
#define QX_DEF_THRESHOLD_ZERO_DOUBLE	1e-16
#define QX_DEF_ENTER					10
#define QX_DEF_BLANK					32
#define QX_DEF_STRING_LENGTH			300

#define PI acos(-1.0)
#define FORWARD 1
#define BACKWARD -1
#define PLANAR 1
#define INTERLEAVED 2
void ctmf(const byte* src, byte* dst, int width, int height, int src_step_row, int dst_step_row, int r, int channels, unsigned long memsize);
/*********************************** HEADER 등. *********************************/
/* BMP FILE & INFORMATION HEADER*/
typedef struct BMPFILEHEADER
{
	WORD	TypeId;			// "BM" - begin of bmp file
	INT		FileSize;		// file size
	WORD	Reserved1;		// for application program
	WORD	Reserved2;		// for application program
	INT		Offset;			// distance to image data
}BmpFileHeader;

typedef struct BMPINFOHEADER
{
	INT		BmpInfoHeaderSize;			//size of BITMAPINFOHEADER struct
	INT		Width;						//image column
	INT		Height;						//image row
	WORD	ColorPlane;					//color plane(always 1)
	WORD	BitCount;					//number of bits per pixel - typical:1, 4, 8, 16, 24, 32
	INT		CompressionFlag;			//compression method - (0:BI_RGB none, 1:BI_RLE8 for 8bit, 2:BIRLE4 for 4bit, 3:BI_BITFIELDS bit 16,32 bit, 4:BI_JPEG for jpg, 5:BI_PNG for png)
	INT		RawDataSize;				//size of raw bitmap data(image size)
	INT		WidthResolution;			//column resolution of image (pixel per meter)
	INT		HeightResolusion;			//row resolution of image (pixel per meter)
	INT		NumOfColor;					//number of colors
	INT		NumOfImportantColor;		//number of important color used(0 when every color is important)
}BmpInforHearder;
/* END BMP FILE & INFORMATION HEADER */

/*********************************** 신호처리 - 영상처리 기술 *********************************/
/* 1D DISCRETE COSINE TRANSFORM - INPUT , OUTPUT , LENGTH OF SIGNAL */
template <class T>
VOID DCT_1D(T *in, T *out, CONST INT N);
/* 1D INVERSE DISCRETE COSINE TRANSFORM - INPUT , OUTPUT , LENGTH OF SIGNAL */
template <class T>
VOID IDCT_1D(T *in, T *out, CONST INT N);
/* 1D (INVERSE) DISCRETE FOURIER TRANSFORM - INPUT , OUTPUT , LENGTH OF SIGNAL , INVERSE(1 - FORWARD, 2 - BACKWARD) */
template <class T>
VOID DFT_1D(T* real, T* imagine, CONST INT Length, CONST INT Forward = FORWARD);

/* PSNR CALCULATE - ORIGINAL IMAGE DATA, TARGET IMAGE DATA, SIZE OF IMAGE, BIT PER PIXEL(DEFAULT : 8BIT=256) */
template<typename T>
DOUBLE PSNR(T* origin, T* target, CONST INT size, CONST INT bitPerPixel = 8);

/* CLIPPING - ARRAY DATA, MIN VALUE, MAX VALUE (DEFAULT : 0 <= data <= 255) */
template<typename T, typename C>
VOID CLIPPING(T* data, CONST INT size, CONST C min = 0.0, CONST C max = 255.0);


/* REVERSE IMAGE - ORIGINAL IMAGE, SIZE OF IMAGE, BIT PER PIXEL(DEFAULT : 8BIT=256) */
template <typename T>
VOID REVERSE_IMAGE(T* data, CONST INT size, CONST INT bitPerPixel = 8);
template <typename T>
VOID REVERSE_IMAGE(T* data, CONST INT row, CONST INT col, CONST INT bitPerPixel = 8);

/* READ & WRITE IMAGE FILE( FORMAT TO RAW DATA ) - READ IMAGE - ROW(HEIGHT), COL(WIDTH), OUTPUT DATA(RAW), IMAGE FILE NAME - WRITE IMAGE - ROW(HEIGHT), COL(WIDTH), IMAGE RAW DATA, IMAGE FILE NAME < color image => rgb plane > */
/* PGM FILE FORMAT */
BYTE* ReadPgm(INT *row, INT *col, LPCSTR  filename);
VOID WritePgm(INT row, INT col, BYTE* img, LPCSTR  filename); //완

/* PBM FILE FORMAT*/
BYTE* ReadPbm(INT *row, INT *col, LPCSTR  filename);
VOID WritePbm(INT row, INT col, BYTE* img, LPCSTR  filename);

/* PPM FILE FORMAT */
BYTE* ReadPpm(INT *row, INT *col, LPCSTR  filename, INT order = 1);
VOID WritePpm(INT row, INT col, BYTE* img, LPCSTR  filename, INT order = 1);

/* BMP FILE FORMAT */
BYTE* ReadBmp(INT *row, INT *col, LPCSTR  filename);
VOID WriteBmp(INT row, INT col, BYTE* img, LPCSTR  filename);

/* END READ & WRITE IMAGE FILE - READ IMAGE - ROW(HEIGHT), COL(WIDTH), OUTPUT DATA(RAW), IMAGE FILE NAME - WRITE IMAGE - ROW(HEIGHT), COL(WIDTH), IMAGE RAW DATA, IMAGE FILE NAME*/

/* PSNR - GET PSNR - */
template <class T>
DOUBLE gPSNR(T* origin, T* target, CONST INT length, CONST INT max = 255);
template <class T>
DOUBLE gPSNR(T* origin, T* target, CONST INT height, CONST INT width, CONST INT max = 255, CONST INT boundary = 0);


/* 개발 준비중인 함수들 */
template<typename T>
DOUBLE SNR(T* origin, T* target, CONST INT size);

BOOL ReadPng(INT *row, INT *col, BYTE* img, LPCSTR  filename);//http://noteroom.tistory.com/157
BOOL WritePng(INT row, INT col, BYTE* img, LPCSTR  filename);//http://www.fastgraph.com/help/jpeg_header_format.html
BOOL ReadJpg(INT *row, INT *col, BYTE* img, LPCSTR  filename);
BOOL WriteJpg(INT row, INT col, BYTE* img, LPCSTR  filename);
BOOL ReadRaw(INT *row, INT *col, BYTE* img, LPCSTR  filename);//있을거고 - 이건 그냥 가져오면 됨 그냥뺄까..
BOOL WriteRaw(INT row, INT col, BYTE* img, LPCSTR  filename);



inline float qx_max_f3(float*a){ return(std::max(std::max(a[0], a[1]), a[2])); }
inline float qx_min_f3(float*a){ return(std::min(std::min(a[0], a[1]), a[2])); }
inline double qx_div(double x, double y){ return((y != 0) ? (x / y) : 0); }
/*Box filter*/
void boxcar_sliding_window_x(double *out, double *in, int h, int w, int radius);
void boxcar_sliding_window_y(double *out, double *in, int h, int w, int radius);
void boxcar_sliding_window(double **out, double **in, double **temp, int h, int w, int radius);
void boxcar_sliding_window(float**out, float**in, float**temp, int h, int w, int radius);
void boxcar_sliding_window(byte**out, byte**in, byte**temp, int h, int w, int radius);
/*Gaussian filter*/
int gaussian_recursive(double **image, double **temp, double sigma, int order, int h, int w);
void gaussian_recursive_x(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn);
void gaussian_recursive_y(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn);
int gaussian_recursive(float **image, float **temp, float sigma, int order, int h, int w);
/*basic functions*/

inline void qx_image_dot_product(double*out, float*a, float*b, int len){ for (int i = 0; i<len; i++)*out++ = double(*a++)*double(*b++); }
inline void qx_image_dot_product(double*out, float*a, byte*b, int len){ for (int i = 0; i<len; i++)*out++ = double(*a++)*double(*b++); }
inline void qx_image_dot_product(double*out, double*a, double*b, int len){ for (int i = 0; i<len; i++)*out++ = (*a++)*(*b++); }
//inline float min(float a,float b){if(a<b) return(a); else return(b);}
//inline float max(float a,float b){if(a>b) return(a); else return(b);}
inline int qx_sum_u3(byte *a) { return(a[0] + a[1] + a[2]); }
inline double qx_sum_d3(double*a){ return(a[0] + a[1] + a[2]); }
inline byte qx_min_u3(byte *a){ return(std::min(std::min(a[0], a[1]), a[2])); }
inline byte qx_max_u3(byte *a){ return(std::max(std::max(a[0], a[1]), a[2])); }
inline byte qx_max_u3(byte r, byte g, byte b){ return(std::max(std::max(r, g), b)); }
inline void image_zero(float **in, int h, int w, float zero = 0){ memset(in[0], zero, sizeof(float)*h*w); }
inline void image_zero(double **in, int h, int w, double zero = 0){ memset(in[0], zero, sizeof(double)*h*w); }
inline void image_zero(byte**in, int h, int w, byte zero = 0){ memset(in[0], zero, sizeof(byte)*h*w); }
inline void image_zero(double ***in, int h, int w, int d, double zero = 0){ memset(in[0][0], zero, sizeof(double)*h*w*d); }
inline byte rgb_2_gray(byte*in){ return(byte(0.299*in[0] + 0.587*in[1] + 0.114*in[2] + 0.5)); }
inline int qx_square_difference_u3(byte *a, byte *b){ int d1, d2, d3; d1 = (*a++) - (*b++); d2 = (*a++) - (*b++);	d3 = (*a++) - (*b++); return(int(d1*d1 + d2*d2 + d3*d3)); }
void qx_specular_free_image(byte ***image_specular_free, byte ***image_normalized, float **diffuse_chromaticity_max, int h, int w);

inline void qx_sort_increase_using_histogram(int*id, byte*image, int len)
{
	int histogram[QX_DEF_CHAR_MAX + 1];
	int nr_bin = QX_DEF_CHAR_MAX + 1;
	memset(histogram, 0, sizeof(int)*nr_bin);
	for (int i = 0; i<len; i++)
	{
		histogram[image[i]]++;
	}
	int nr_hitted_prev = histogram[0];
	histogram[0] = 0;
	for (int k = 1; k<nr_bin; k++)
	{
		int nr_hitted = histogram[k];
		histogram[k] = nr_hitted_prev + histogram[k - 1];
		nr_hitted_prev = nr_hitted;
	}
	for (int i = 0; i<len; i++)
	{
		byte dist = image[i];
		int index = histogram[dist]++;
		id[index] = i;
	}
}
inline double *get_color_weighted_table(double sigma_range, int len)
{
	double *table_color, *color_table_x; int y;
	table_color = new double[len];
	color_table_x = &table_color[0];
	for (y = 0; y<len; y++) (*color_table_x++) = exp(-double(y*y) / (2 * sigma_range*sigma_range));
	return(table_color);
}
inline void color_weighted_table_update(double *table_color, double dist_color, int len)
{
	double *color_table_x; int y;
	color_table_x = &table_color[0];
	for (y = 0; y<len; y++) (*color_table_x++) = exp(-double(y*y) / (2 * dist_color*dist_color));
}

inline void vec_min_val(float &min_val, float *in, int len)
{
	min_val = in[0];
	for (int i = 1; i<len; i++) if (in[i]<min_val) min_val = in[i];
}
inline void vec_min_val(byte &min_val, byte *in, int len)
{
	min_val = in[0];
	for (int i = 1; i<len; i++) if (in[i]<min_val) min_val = in[i];
}
inline void vec_max_val(float &max_val, float *in, int len)
{
	max_val = in[0];
	for (int i = 1; i<len; i++) if (in[i]>max_val) max_val = in[i];
}
inline void vec_max_val(byte &max_val, byte *in, int len)
{
	max_val = in[0];
	for (int i = 1; i<len; i++) if (in[i]>max_val) max_val = in[i];
}
inline void down_sample_1(byte **out, byte **in, int h, int w, int scale_exp)
{
	int y, x; int ho, wo; byte *out_y, *in_x;
	ho = (h >> scale_exp); wo = (w >> scale_exp);
	for (y = 0; y<ho; y++)
	{
		out_y = &out[y][0]; in_x = in[y << scale_exp];
		for (x = 0; x<wo; x++) *out_y++ = in_x[x << scale_exp];
	}
}
inline void down_sample_1(float**out, float**in, int h, int w, int scale_exp)
{
	int y, x; int ho, wo; float *out_y, *in_x;
	ho = (h >> scale_exp); wo = (w >> scale_exp);
	for (y = 0; y<ho; y++)
	{
		out_y = &out[y][0]; in_x = in[y << scale_exp];
		for (x = 0; x<wo; x++) *out_y++ = in_x[x << scale_exp];
	}
}
inline double qx_linear_interpolate_xy(double **image, double x, double y, int h, int w)
{
	int x0, xt, y0, yt; double dx, dy, dx1, dy1, d00, d0t, dt0, dtt;
	x0 = int(x); xt = min(x0 + 1, w - 1); y0 = int(y); yt = min(y0 + 1, h - 1);
	dx = x - x0; dy = y - y0; dx1 = 1 - dx; dy1 = 1 - dy; d00 = dx1*dy1; d0t = dx*dy1; dt0 = dx1*dy; dtt = dx*dy;
	return(d00*image[y0][x0] + d0t*image[y0][xt] + dt0*image[yt][x0] + dtt*image[yt][xt]);
}
/*memory*/
inline double *** qx_allocd_3(int n, int r, int c, int padding = QX_DEF_PADDING)
{
	double *a, **p, ***pp;
	int rc = r*c;
	int i, j;
	a = (double*)malloc(sizeof(double)*(n*rc + padding));
	if (a == NULL) { printf("qx_allocd_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (double**)malloc(sizeof(double*)*n*r);
	pp = (double***)malloc(sizeof(double**)*n);
	for (i = 0; i<n; i++)
	for (j = 0; j<r; j++)
		p[i*r + j] = &a[i*rc + j*c];
	for (i = 0; i<n; i++)
		pp[i] = &p[i*r];
	return(pp);
}
inline void qx_freed_3(double ***p)
{
	if (p != NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline byte** qx_allocu(int r, int c, int padding = QX_DEF_PADDING)
{
	byte *a, **p;
	a = (byte*)malloc(sizeof(byte)*(r*c + padding));
	if (a == NULL) { printf("qx_allocu() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (byte**)malloc(sizeof(byte*)*r);
	for (int i = 0; i<r; i++) p[i] = &a[i*c];
	return(p);
}
inline void qx_freeu(byte **p)
{
	if (p != NULL)
	{
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline byte *** qx_allocu_3(int n, int r, int c, int padding = QX_DEF_PADDING)
{
	byte *a, **p, ***pp;
	int rc = r*c;
	int i, j;
	a = (byte*)malloc(sizeof(byte)*(n*rc + padding));
	if (a == NULL) { printf("qx_allocu_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (byte**)malloc(sizeof(byte*)*n*r);
	pp = (byte***)malloc(sizeof(byte**)*n);
	for (i = 0; i<n; i++)
	for (j = 0; j<r; j++)
		p[i*r + j] = &a[i*rc + j*c];
	for (i = 0; i<n; i++)
		pp[i] = &p[i*r];
	return(pp);
}
inline void qx_freeu_3(byte ***p)
{
	if (p != NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline void qx_freeu_1(byte*p)
{
	if (p != NULL)
	{
		delete[] p;
		p = NULL;
	}
}
inline float** qx_allocf(int r, int c, int padding = QX_DEF_PADDING)
{
	float *a, **p;
	a = (float*)malloc(sizeof(float)*(r*c + padding));
	if (a == NULL) { printf("qx_allocf() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (float**)malloc(sizeof(float*)*r);
	for (int i = 0; i<r; i++) p[i] = &a[i*c];
	return(p);
}
inline void qx_freef(float **p)
{
	if (p != NULL)
	{
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline float *** qx_allocf_3(int n, int r, int c, int padding = QX_DEF_PADDING)
{
	float *a, **p, ***pp;
	int rc = r*c;
	int i, j;
	a = (float*)malloc(sizeof(float)*(n*rc + padding));
	if (a == NULL) { printf("qx_allocf_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (float**)malloc(sizeof(float*)*n*r);
	pp = (float***)malloc(sizeof(float**)*n);
	for (i = 0; i<n; i++)
	for (j = 0; j<r; j++)
		p[i*r + j] = &a[i*rc + j*c];
	for (i = 0; i<n; i++)
		pp[i] = &p[i*r];
	return(pp);
}
inline void qx_freef_3(float ***p)
{
	if (p != NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline int** qx_alloci(int r, int c, int padding = QX_DEF_PADDING)
{
	int *a, **p;
	a = (int*)malloc(sizeof(int)*(r*c + padding));
	if (a == NULL) { printf("qx_alloci() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (int**)malloc(sizeof(int*)*r);
	for (int i = 0; i<r; i++) p[i] = &a[i*c];
	return(p);
}
inline void qx_freei(int **p)
{
	if (p != NULL)
	{
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline void qx_freei_1(int*p)
{
	if (p != NULL)
	{
		delete[] p;
		p = NULL;
	}
}
inline double** qx_allocd(int r, int c, int padding = QX_DEF_PADDING)
{
	double *a, **p;
	a = (double*)malloc(sizeof(double)*(r*c + padding));
	if (a == NULL) { printf("qx_allocd() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (double**)malloc(sizeof(double*)*r);
	for (int i = 0; i<r; i++) p[i] = &a[i*c];
	return(p);
}
inline void qx_freed(double **p)
{
	if (p != NULL)
	{
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline byte**** qx_allocu_4(int t, int n, int r, int c, int padding = QX_DEF_PADDING)
{
	byte *a, **p, ***pp, ****ppp;
	int nrc = n*r*c, nr = n*r, rc = r*c;
	int i, j, k;
	a = (byte*)malloc(sizeof(byte)*(t*nrc + padding));
	if (a == NULL) { printf("qx_allocu_4() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (byte**)malloc(sizeof(byte*)*t*nr);
	pp = (byte***)malloc(sizeof(byte**)*t*n);
	ppp = (byte****)malloc(sizeof(byte***)*t);
	for (k = 0; k<t; k++)
	for (i = 0; i<n; i++)
	for (j = 0; j<r; j++)
		p[k*nr + i*r + j] = &a[k*nrc + i*rc + j*c];
	for (k = 0; k<t; k++)
	for (i = 0; i<n; i++)
		pp[k*n + i] = &p[k*nr + i*r];
	for (k = 0; k<t; k++)
		ppp[k] = &pp[k*n];
	return(ppp);
}
inline void qx_freeu_4(byte ****p)
{
	if (p != NULL)
	{
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
inline double**** qx_allocd_4(int t, int n, int r, int c, int padding = QX_DEF_PADDING)
{
	double *a, **p, ***pp, ****ppp;
	int nrc = n*r*c, nr = n*r, rc = r*c;
	int i, j, k;
	a = (double*)malloc(sizeof(double)*(t*nrc + padding));
	if (a == NULL) { printf("qx_allocd_4() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p = (double**)malloc(sizeof(double*)*t*nr);
	pp = (double***)malloc(sizeof(double**)*t*n);
	ppp = (double****)malloc(sizeof(double***)*t);
	for (k = 0; k<t; k++)
	for (i = 0; i<n; i++)
	for (j = 0; j<r; j++)
		p[k*nr + i*r + j] = &a[k*nrc + i*rc + j*c];
	for (k = 0; k<t; k++)
	for (i = 0; i<n; i++)
		pp[k*n + i] = &p[k*nr + i*r];
	for (k = 0; k<t; k++)
		ppp[k] = &pp[k*n];
	return(ppp);
}
inline void qx_freed_4(double ****p)
{

	if (p != NULL)
	{
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
void qx_stereo_flip_corr_vol(double***corr_vol_right, double***corr_vol, int h, int w, int nr_plane);
inline void qx_memcpy_u3(byte a[3], byte b[3]){ *a++ = *b++; *a++ = *b++; *a++ = *b++; }
inline void image_copy(double***out, double***in, int h, int w, int d){ memcpy(out[0][0], in[0][0], sizeof(double)*h*w*d); }
inline void image_copy(byte**out, byte**in, int h, int w){ memcpy(out[0], in[0], sizeof(byte)*h*w); }
void depth_best_cost(byte**depth, double***evidence, int h, int w, int nr_planes);
void vec_min_pos(int &min_pos, double *in, int len);
void qx_detect_occlusion_left_right(byte**mask_left, byte**depth_left, byte**depth_right, int h, int w, int nr_plane);
int file_open_ascii(char *file_path, int *out, int len);


class qx_queue_mst
{
public:
	qx_queue_mst();
	~qx_queue_mst();
	void clean();
	int	init(int len);
	void reinit();
	void push(int x);
	int pull();//return x: the first item in the queue;
	int *queue;
	int	length;
	int	first;
	int	last;
};
class qx_mst_kruskals_image
{
public:
	qx_mst_kruskals_image();
	~qx_mst_kruskals_image();
	void clean();
	int init(int h, int w, int nr_channel, int nr_neighbor = QX_DEF_MST_KI_4NR_NEIGHBOR);
	int mst(byte*image, bool print_edges = false);
public://output
	void update_table(double sigma_range);
	int*get_rank(){ return(m_rank); }
	int*get_parent(){ return(m_parent); }
	int*get_nr_child(){ return(m_nr_child); }
	int**get_children(){ return(m_children); }
	byte*get_weight(){ return(m_weight); }
	int*get_node_id(){ return(m_node_id_from_parent_to_child); }
private:
	int m_h, m_w, m_nr_channel, m_nr_neighbor;
	int m_nr_edge, m_nr_vertices, m_parent_default, m_tree_size;
	int**m_edge, *m_id_edge;//compute edges
	byte*m_distance;
	byte*m_image;

	int*m_parent;//obtain tree edges
	int*m_nr_child, m_max_nr_child, **m_children;
	byte*m_weight;
	int**m_connected;
	byte**m_connected_distance;
	int*m_nr_connected;

	qx_queue_mst m_queue;//build tree
	int*m_node_id_from_parent_to_child;
	int*m_rank;
private:
	void init_mst();
	int findset(int x);
	void kruskal();
	int build_tree();
	//void reorder_edges(vector< pair< int, QX_PAIR > >&min_spanning_tree);
	//void print(vector< pair< int, QX_PAIR > >&min_spanning_tree);
};
void test_qx_mst_kruskals_image();



/*basic function*/
/*get the size of the target image, and the properties of the image, and use the following
4 functions to read the image data*/
byte* loadimage(char* file_name, int &h, int &w, bool &is_ppm);
int loadimage(byte* out, char* file_name, int h, int w);
int loadimage(float* out, byte *out_u, char* file_name, int h, int w);
void get_ascii_pgm(byte *image, FILE* file_in, int h, int w);
void get_ascii_ppm(byte *image, FILE* file_in, int h, int w);
byte* get_binary_ppm(FILE* file_in, int h, int w);/*If the image is a ppm image in binary mode*/
byte* get_binary_pgm(FILE* file_in, int h, int w);/*If the image is a pgm image in binary mode*/
byte* get_ascii_ppm(FILE* file_in, int h, int w);/*If the image is a ppm image in ascii mode*/
byte* get_ascii_pgm(FILE* file_in, int h, int w);/*If the image is a pgm image in ascii mode*/
/*both saveimage_ppm and save_image_pgm use the following Four functions to save an image*/
void saveimage_ppm(char *file_name, byte *image, int h, int w, bool is_binary = true);
void saveimage_pgm(char *file_name, byte *image, int h, int w, bool is_binary = true);
void write_binary_ppm(char* file_name, byte* image, int h, int w);/*If the image is a ppm image in binary mode*/
void write_binary_pgm(char* file_name, byte* image, int h, int w);/*If the image is a pgm image in binary mode*/
void write_ascii_ppm(char* file_name, byte *image, int h, int w);/*If the image is a ppm image in ascii mode*/
void write_ascii_pgm(char* file_name, byte *image, int h, int w);/*If the image is a pgm image in ascii mode*/

/*extended function*/
/*To get the size of the image*/
void qx_image_size(char* file_name, int &h, int &w, int *nr_channel = NULL);
/*load image (with height h and width w) as a float matrix
(if image is color image, do color to gray conversion)*/
float** loadimage_pgm(char *file_name, int &h, int &w);
/*load ppm image (with height h and width w) as a float cubic (size: h x w x 3)*/
float*** loadimage_ppm(char *file_name, int &h, int &w);
/*load ppm image (with height h and width w) as a byte cubic (size: h x w x 3)*/
byte*** loadimage_ppm_u(char *file_name, int &h, int &w);
/*load pgm image as a float matrix, and save the height as h, and width as w*/
int** loadimage_pgm_i(char *file_name, int &h, int &w);
/*load image (with height h and width w) as a unsigned matrix
(if image is color image, do color to gray conversion)*/
byte** loadimage_pgm_u(char *file_name, int &h, int &w);
/*The following Four functions are used to save an image of different data types*/
void saveimage_pgm(char *file_name, float **image, int h, int w, int scale = 1);
void saveimage_pgm(char *file_name, byte **image, int h, int w, int scale = 1);
void saveimage_pgm_ascii(char *file_name, float **image, int h, int w, int scale = 1);
void saveimage_pgm_ascii(char *file_name, int **image, int h, int w, int scale = 1);
void saveimage_pgm(char *file_name, int **image, int h, int w, int scale = 1);
void saveimage_ppm(char *file_name, float ***image, int h, int w, int scale = 1);
void saveimage_ppm(char *file_name, double ***image, int h, int w, int scale = 1);
void saveimage_ppm(char *file_name, byte ***image, int h, int w, int scale = 1);


/*load image*/
int qx_loadimage(char*file_name, byte*image, int h, int w, int *nr_channel = NULL);
int qx_loadimage(char*file_name, short*image, int h, int w, int *nr_channel = NULL);
int qx_loadimage(char*file_name, float*image, int h, int w, int *nr_channel = NULL);
/*save image*/
void qx_saveimage(char*file_name, float*image, int h, int w, int channel);
void qx_saveimage(char*file_name, double*image_d, int h, int w, int channel);
void qx_saveimage(char*file_name, short*image, int h, int w, int channel);
void qx_saveimage(char*file_name, byte*image, int h, int w, int channel);
float *loadimage(char* file_name, int &h, int &w, int *is_ppm = NULL);
int loadimage(char* file_name, float *image, int &h, int &w, int *nr_channel = NULL);
float *get_binary(FILE* file_in, int h, int w, int channel);
void saveimage(char* file_name, float *image, int h, int w, int channel = 1);
void saveimage(char* file_name, double *image_d, int h, int w, int channel);
float ***loadimage_ftif(char *file_name, int &h, int &w, int &nr_channel);



class qx_tree_filter
{
public:
	qx_tree_filter();
	~qx_tree_filter();
	void clean();
	int init(int h, int w, int nr_channel, double sigma_range = QX_DEF_MST_KI_SIGMA_RANGE, int nr_neighbor = QX_DEF_MST_KI_4NR_NEIGHBOR);
	int build_tree(byte*texture);
	//void init_tree_value(byte*image);
	template<typename T>void init_tree_value(T*image, bool compute_weight);
	template<typename T>void combine_tree(T*image_filtered);
	int filter(double*cost, double*cost_backup, int nr_plane);
	int*get_rank(){ return(m_mst_rank); };
	void update_table(double sigma_range);
private:
	qx_mst_kruskals_image m_mst;
	int m_h, m_w, m_nr_channel; int m_nr_pixel;
	int*m_mst_parent;
	int*m_mst_nr_child;
	int**m_mst_children;//[QX_DEF_MST_NODE_MAX_NR_CHILDREN];
	int*m_mst_rank;
	byte*m_mst_weight;//cost between this node and its parent

	double m_table[QX_DEF_CHAR_MAX + 1];
	int*m_node_id;
private:
	void filter_main(bool compute_weight);
};
void test_qx_tree_filter();


