#include	"mental.h"
#include <Windows.h>//The win32 API library 
INT Read_Char_For_P(std::ifstream& file);
INT Read_Int_For_P(std::ifstream& file);
#pragma warning(disable:4996)

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Type declarations */
#ifdef _MSC_VER
#include <basetsd.h>
typedef UINT8 uint8_t;
typedef UINT16 uint16_t;
typedef UINT32 uint32_t;
#pragma warning( disable: 4799 )
#else
#include <stdint.h>
#endif

/* Intrinsic declarations */
#if defined(__SSE2__) || defined(__MMX__)
#if defined(__SSE2__)
#include <emmintrin.h>
#elif defined(__MMX__)
#include <mmintrin.h>
#endif
#if defined(__GNUC__)
#include <mm_malloc.h>
#elif defined(_MSC_VER)
#include <malloc.h>
#endif
#elif defined(__ALTIVEC__)
#include <altivec.h>
#endif

/* Compiler peculiarities */
#if defined(__GNUC__)
#include <stdint.h>
#define inline __inline__
#define align(x) __attribute__ ((aligned (x)))
#elif defined(_MSC_VER)
#define inline __inline
#define align(x) __declspec(align(x))
#else
#define inline
#define align(x)
#endif

#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#endif

/**
* This structure represents a two-tier histogram. The first tier (known as the
* "coarse" level) is 4 bit wide and the second tier (known as the "fine" level)
* is 8 bit wide. Pixels inserted in the fine level also get inserted into the
* coarse bucket designated by the 4 MSBs of the fine bucket value.
*
* The structure is aligned on 16 bytes, which is a prerequisite for SIMD
* instructions. Each bucket is 16 bit wide, which means that extra care must be
* taken to prevent overflow.
*/
typedef struct align(16)
{
	uint16_t coarse[16];
	uint16_t fine[16][16];
} Histogram;

/**
* HOP is short for Histogram OPeration. This macro makes an operation \a op on
* histogram \a h for pixel value \a x. It takes care of handling both levels.
*/
#define HOP(h,x,op) \
	h.coarse[x >> 4] op; \
	*((uint16_t*)h.fine + x) op;

#define COP(c,j,x,op) \
	h_coarse[16 * (n*c + j) + (x >> 4)] op; \
	h_fine[16 * (n*(16 * c + (x >> 4)) + j) + (x & 0xF)] op;

/**
* Adds histograms \a x and \a y and stores the result in \a y. Makes use of
* SSE2, MMX or Altivec, if available.
*/
#if defined(__SSE2__)
static inline void histogram_add(const uint16_t x[16], uint16_t y[16])
{
	*(__m128i*) &y[0] = _mm_add_epi16(*(__m128i*) &y[0], *(__m128i*) &x[0]);
	*(__m128i*) &y[8] = _mm_add_epi16(*(__m128i*) &y[8], *(__m128i*) &x[8]);
}
#elif defined(__MMX__)
static inline void histogram_add(const uint16_t x[16], uint16_t y[16])
{
	*(__m64*) &y[0] = _mm_add_pi16(*(__m64*) &y[0], *(__m64*) &x[0]);
	*(__m64*) &y[4] = _mm_add_pi16(*(__m64*) &y[4], *(__m64*) &x[4]);
	*(__m64*) &y[8] = _mm_add_pi16(*(__m64*) &y[8], *(__m64*) &x[8]);
	*(__m64*) &y[12] = _mm_add_pi16(*(__m64*) &y[12], *(__m64*) &x[12]);
}
#elif defined(__ALTIVEC__)
static inline void histogram_add(const uint16_t x[16], uint16_t y[16])
{
	*(vector unsigned short*) &y[0] = vec_add(*(vector unsigned short*) &y[0], *(vector unsigned short*) &x[0]);
	*(vector unsigned short*) &y[8] = vec_add(*(vector unsigned short*) &y[8], *(vector unsigned short*) &x[8]);
}
#else
static inline void histogram_add(const uint16_t x[16], uint16_t y[16])
{
	int i;
	for (i = 0; i < 16; ++i) {
		y[i] += x[i];
	}
}
#endif

/**
* Subtracts histogram \a x from \a y and stores the result in \a y. Makes use
* of SSE2, MMX or Altivec, if available.
*/
#if defined(__SSE2__)
static inline void histogram_sub(const uint16_t x[16], uint16_t y[16])
{
	*(__m128i*) &y[0] = _mm_sub_epi16(*(__m128i*) &y[0], *(__m128i*) &x[0]);
	*(__m128i*) &y[8] = _mm_sub_epi16(*(__m128i*) &y[8], *(__m128i*) &x[8]);
}
#elif defined(__MMX__)
static inline void histogram_sub(const uint16_t x[16], uint16_t y[16])
{
	*(__m64*) &y[0] = _mm_sub_pi16(*(__m64*) &y[0], *(__m64*) &x[0]);
	*(__m64*) &y[4] = _mm_sub_pi16(*(__m64*) &y[4], *(__m64*) &x[4]);
	*(__m64*) &y[8] = _mm_sub_pi16(*(__m64*) &y[8], *(__m64*) &x[8]);
	*(__m64*) &y[12] = _mm_sub_pi16(*(__m64*) &y[12], *(__m64*) &x[12]);
}
#elif defined(__ALTIVEC__)
static inline void histogram_sub(const uint16_t x[16], uint16_t y[16])
{
	*(vector unsigned short*) &y[0] = vec_sub(*(vector unsigned short*) &y[0], *(vector unsigned short*) &x[0]);
	*(vector unsigned short*) &y[8] = vec_sub(*(vector unsigned short*) &y[8], *(vector unsigned short*) &x[8]);
}
#else
static inline void histogram_sub(const uint16_t x[16], uint16_t y[16])
{
	int i;
	for (i = 0; i < 16; ++i) {
		y[i] -= x[i];
	}
}
#endif

static inline void histogram_muladd(const uint16_t a, const uint16_t x[16],
	uint16_t y[16])
{
	int i;
	for (i = 0; i < 16; ++i) {
		y[i] += a * x[i];
	}
}

static void ctmf_helper(
	const byte* const src, byte* const dst,
	const int width, const int height,
	const int src_step, const int dst_step,
	const int r, const int cn,
	const int pad_left, const int pad_right
	)
{
	const int m = height, n = width;
	int i, j, k, c;
	const byte *p, *q;

	Histogram H[4];
	uint16_t *h_coarse, *h_fine, luc[4][16];

	assert(src);
	assert(dst);
	assert(r >= 0);
	assert(width >= 2 * r + 1);
	assert(height >= 2 * r + 1);
	assert(src_step != 0);
	assert(dst_step != 0);

	/* SSE2 and MMX need aligned memory, provided by _mm_malloc(). */
#if defined(__SSE2__) || defined(__MMX__)
	h_coarse = (uint16_t*)_mm_malloc(1 * 16 * n * cn * sizeof(uint16_t), 16);
	h_fine = (uint16_t*)_mm_malloc(16 * 16 * n * cn * sizeof(uint16_t), 16);
	memset(h_coarse, 0, 1 * 16 * n * cn * sizeof(uint16_t));
	memset(h_fine, 0, 16 * 16 * n * cn * sizeof(uint16_t));
#else
	h_coarse = (uint16_t*)calloc(1 * 16 * n * cn, sizeof(uint16_t));
	h_fine = (uint16_t*)calloc(16 * 16 * n * cn, sizeof(uint16_t));
#endif

	/* First row initialization */
	for (j = 0; j < n; ++j) {
		for (c = 0; c < cn; ++c) {
			COP(c, j, src[cn*j + c], += r + 1);
		}
	}
	for (i = 0; i < r; ++i) {
		for (j = 0; j < n; ++j) {
			for (c = 0; c < cn; ++c) {
				COP(c, j, src[src_step*i + cn*j + c], ++);
			}
		}
	}

	for (i = 0; i < m; ++i) {

		/* Update column histograms for entire row. */
		p = src + src_step * MAX(0, i - r - 1);
		q = p + cn * n;
		for (j = 0; p != q; ++j) {
			for (c = 0; c < cn; ++c, ++p) {
				COP(c, j, *p, --);
			}
		}

		p = src + src_step * MIN(m - 1, i + r);
		q = p + cn * n;
		for (j = 0; p != q; ++j) {
			for (c = 0; c < cn; ++c, ++p) {
				COP(c, j, *p, ++);
			}
		}

		/* First column initialization */
		memset(H, 0, cn*sizeof(H[0]));
		memset(luc, 0, cn*sizeof(luc[0]));
		if (pad_left) {
			for (c = 0; c < cn; ++c) {
				histogram_muladd(r, &h_coarse[16 * n*c], H[c].coarse);
			}
		}
		for (j = 0; j < (pad_left ? r : 2 * r); ++j) {
			for (c = 0; c < cn; ++c) {
				histogram_add(&h_coarse[16 * (n*c + j)], H[c].coarse);
			}
		}
		for (c = 0; c < cn; ++c) {
			for (k = 0; k < 16; ++k) {
				histogram_muladd(2 * r + 1, &h_fine[16 * n*(16 * c + k)], &H[c].fine[k][0]);
			}
		}

		for (j = pad_left ? 0 : r; j < (pad_right ? n : n - r); ++j) {
			for (c = 0; c < cn; ++c) {
				const uint16_t t = 2 * r*r + 2 * r;
				uint16_t sum = 0, *segment;
				int b;

				histogram_add(&h_coarse[16 * (n*c + MIN(j + r, n - 1))], H[c].coarse);

				/* Find median at coarse level */
				for (k = 0; k < 16; ++k) {
					sum += H[c].coarse[k];
					if (sum > t) {
						sum -= H[c].coarse[k];
						break;
					}
				}
				assert(k < 16);

				/* Update corresponding histogram segment */
				if (luc[c][k] <= j - r) {
					memset(&H[c].fine[k], 0, 16 * sizeof(uint16_t));
					for (luc[c][k] = j - r; luc[c][k] < MIN(j + r + 1, n); ++luc[c][k]) {
						histogram_add(&h_fine[16 * (n*(16 * c + k) + luc[c][k])], H[c].fine[k]);
					}
					if (luc[c][k] < j + r + 1) {
						histogram_muladd(j + r + 1 - n, &h_fine[16 * (n*(16 * c + k) + (n - 1))], &H[c].fine[k][0]);
						luc[c][k] = j + r + 1;
					}
				}
				else {
					for (; luc[c][k] < j + r + 1; ++luc[c][k]) {
						histogram_sub(&h_fine[16 * (n*(16 * c + k) + MAX(luc[c][k] - 2 * r - 1, 0))], H[c].fine[k]);
						histogram_add(&h_fine[16 * (n*(16 * c + k) + MIN(luc[c][k], n - 1))], H[c].fine[k]);
					}
				}

				histogram_sub(&h_coarse[16 * (n*c + MAX(j - r, 0))], H[c].coarse);

				/* Find median in segment */
				segment = H[c].fine[k];
				for (b = 0; b < 16; ++b) {
					sum += segment[b];
					if (sum > t) {
						dst[dst_step*i + cn*j + c] = 16 * k + b;
						break;
					}
				}
				assert(b < 16);
			}
		}
	}

#if defined(__SSE2__) || defined(__MMX__)
	_mm_empty();
	_mm_free(h_coarse);
	_mm_free(h_fine);
#else
	free(h_coarse);
	free(h_fine);
#endif
}

/**
* \brief Constant-time median filtering
*
* This function does a median filtering of an 8-bit image. The source image is
* processed as if it was padded with zeros. The median kernel is square with
* odd dimensions. Images of arbitrary size may be processed.
*
* To process multi-channel images, you must call this function multiple times,
* changing the source and destination adresses and steps such that each channel
* is processed as an independent single-channel image.
*
* Processing images of arbitrary bit depth is not supported.
*
* The computing time is O(1) per pixel, independent of the radius of the
* filter. The algorithm's initialization is O(r*width), but it is negligible.
* Memory usage is simple: it will be as big as the cache size, or smaller if
* the image is small. For efficiency, the histograms' bins are 16-bit wide.
* This may become too small and lead to overflow as \a r increases.
*
* \param src           Source image data.
* \param dst           Destination image data. Must be preallocated.
* \param width         Image width, in pixels.
* \param height        Image height, in pixels.
* \param src_step      Distance between adjacent pixels on the same column in
*                      the source image, in bytes.
* \param dst_step      Distance between adjacent pixels on the same column in
*                      the destination image, in bytes.
* \param r             Median filter radius. The kernel will be a 2*r+1 by
*                      2*r+1 square.
* \param cn            Number of channels. For example, a grayscale image would
*                      have cn=1 while an RGB image would have cn=3.
* \param memsize       Maximum amount of memory to use, in bytes. Set this to
*                      the size of the L2 cache, then vary it slightly and
*                      measure the processing time to find the optimal value.
*                      For example, a 512 kB L2 cache would have
*                      memsize=512*1024 initially.
*/
void ctmf(
	const byte* const src, byte* const dst,
	const int width, const int height,
	const int src_step, const int dst_step,
	const int r, const int cn, const long unsigned int memsize
	)
{
	/*
	* Processing the image in vertical stripes is an optimization made
	* necessary by the limited size of the CPU cache. Each histogram is 544
	* bytes big and therefore I can fit a limited number of them in the cache.
	* That number may sometimes be smaller than the image width, which would be
	* the number of histograms I would need without stripes.
	*
	* I need to keep histograms in the cache so that they are available
	* quickly when processing a new row. Each row needs access to the previous
	* row's histograms. If there are too many histograms to fit in the cache,
	* thrashing to RAM happens.
	*
	* To solve this problem, I figure out the maximum number of histograms
	* that can fit in cache. From this is determined the number of stripes in
	* an image. The formulas below make the stripes all the same size and use
	* as few stripes as possible.
	*
	* Note that each stripe causes an overlap on the neighboring stripes, as
	* when mowing the lawn. That overlap is proportional to r. When the overlap
	* is a significant size in comparison with the stripe size, then we are not
	* O(1) anymore, but O(r). In fact, we have been O(r) all along, but the
	* initialization term was neglected, as it has been (and rightly so) in B.
	* Weiss, "Fast Median and Bilateral Filtering", SIGGRAPH, 2006. Processing
	* by stripes only makes that initialization term bigger.
	*
	* Also, note that the leftmost and rightmost stripes don't need overlap.
	* A flag is passed to ctmf_helper() so that it treats these cases as if the
	* image was zero-padded.
	*/
	int stripes = (int)ceil((double)(width - 2 * r) / (memsize / sizeof(Histogram)-2 * r));
	int stripe_size = (int)ceil((double)(width + stripes * 2 * r - 2 * r) / stripes);

	int i;

	for (i = 0; i < width; i += stripe_size - 2 * r) {
		int stripe = stripe_size;
		/* Make sure that the filter kernel fits into one stripe. */
		if (i + stripe_size - 2 * r >= width || width - (i + stripe_size - 2 * r) < 2 * r + 1) {
			stripe = width - i;
		}

		ctmf_helper(src + cn*i, dst + cn*i, stripe, height, src_step, dst_step, r, cn,
			i == 0, stripe == width - i);

		if (stripe == width - i) {
			break;
		}
	}
}

#define LEN_MAX	256
template <class T>
VOID DCT_1D(T *in, T *out, CONST INT N)
{
	register INT k, n;
	T dct = 0;
	T temp = 0;
	T omega = 0;
	for (k = 1; k <= N; ++k) // 결과 Loop
	{
		dct = 0;
		for (n = 1; n <= N; ++n) // sigma
		{
			if (k == 1)
				omega = 1.0 / sqrt((T)N);
			else // 2 <= k <= N
			{
				omega = 2.0 / N;
				omega = sqrt(omega);
			}

			temp = (PI * (2.0 * n - 1) * (k - 1)) / (2 * N);
			dct += omega * in[n - 1] * cos(temp);
		}
		out[k - 1] = dct;
	}

}

template <class T>
VOID IDCT_1D(T *in, T *out, CONST INT N)
{
	register INT k, n;
	T dct = 0;
	T temp = 0;
	T omega = 0;
	for (n = 1; n <= N; ++n) // 결과 Loop
	{
		dct = 0;
		for (k = 1; k <= N; ++k) // sigma
		{
			if (k == 1)
				omega = 1.0 / sqrt((T)N);
			else // 2 <= k <= N
			{
				omega = 2.0 / N;
				omega = sqrt(omega);
			}

			temp = (PI * (2.0 * n - 1) * (k - 1)) / (2 * N);
			dct += omega * in[k - 1] * cos(temp);
		}
		out[n - 1] = dct;
	}


}

template <class T>
VOID DFT_1D(T* real, T* imagine, CONST INT Length, CONST INT Forward)
{
	register INT i, j, k, n;
	INT N = Length;

	T* tr = new T[sizeof(T)*N]; // temp re 
	T* ti = new T[sizeof(T)*N]; // temp im

	memcpy(tr, real, sizeof(T)*N);
	memcpy(ti, imagine, sizeof(T)*N);

	T re, im, temp;
	T cosine, sine;

	for (k = 1; k <= N; ++k)
	{
		re = 0.0;
		im = 0.0;
		for (j = 1; j <= N; ++j)
		{
			temp = 2 * Forward*PI*(j - 1)*(k - 1) / (N); //-j
			cosine = cos(temp);
			sine = sin(temp);
			re += (tr[j - 1] * cosine + ti[j - 1] * sine);
			im += (ti[j - 1] * cosine - tr[j - 1] * sine);
		}
		real[k - 1] = re;
		imagine[k - 1] = im;
	}

	if (Forward == -1) // IDFT
	{
		for (i = 0; i < N; i++)
		{
			real[i] /= (T)N;
			imagine[i] /= (T)N;
		}
	}

	delete[] tr;
	delete[] ti;
}

template <class T>
DOUBLE PSNR(T* origin, T* target, CONST INT size, INT bitPerPixel)
{
	DOUBLE psnr = 0.0;
	INT Peak = pow(2, bitPerPixel);
	DOUBLE mse = MSE(origin, target, size);

	psnr = 10.0 * log10((Peak*Peak) / mse);

	return psnr;
}

template <class T>
DOUBLE MSE(T* origin, T* target, CONST INT size)
{
	register INT i;
	DOUBLE mse = 0.0;

	for (i = 0; i < size; ++i)
		mse += (origin[i] - target[i]) * (origin[i] - target[i]);

	mse /= size;

	return mse;
}

template<typename T, typename C>
VOID CLIPPING(T* data, CONST INT size, CONST C min, CONST C max)
{
	register INT i;
	for (i = 0; i < size; ++i)
	{
		if (data[i] > max)
			data[i] = max;
		else if (data[i] < min)
			data[i] = min;
	}
}


template <typename T>
VOID REVERSE_IMAGE(T* data, CONST INT size, CONST INT bitPerPixel)
{
	register INT i;
	INT max = pow(2, bitPerPixel);
	for (i = 0; i < size; ++i)
		data[i] = max - data[i];

}
template <typename T>
VOID REVERSE_IMAGE(T* data, CONST INT row, CONST INT col, CONST INT bitPerPixel)
{
	INT size = row*col;
	register INT i;
	INT max = pow(2, bitPerPixel);
	for (i = 0; i < size; ++i)
		data[i] = max - data[i];
}

VOID WritePgm(INT row, INT col, BYTE* img, LPCSTR  filename)
{
	std::ofstream file;
	char* header = "P5\n# Created by Vision Image Processing(VIP) Lab, Sangmyung Univ. Korea \n";

	file.open(filename, std::ios::out | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in WritePgm function" << std::endl;
		exit(1);
	}

	//header information writing phase 'P5' is for pgm file format
	file << header << col << ' ' << row << std::endl << "255" << std::endl;
	file.write((char*)img, row*col);

	file.close();

}

BYTE* ReadPgm(INT *row, INT *col, LPCSTR  filename)
{
	std::ifstream file;
	BYTE *buffer, *handle;
	INT size;

	file.open(filename, std::ios::in | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in ReadPgm function" << std::endl;
		exit(1);
	}

	//header information reading phase 'P5' is for pgm file format
	if ((file.get() != 'P') || (file.get() != '5'))
	{
		std::cerr << "error : image format is not a pgm" << std::endl;
		exit(1);
	}

	//read col & row
	*col = Read_Int_For_P(file);
	*row = Read_Int_For_P(file);
	Read_Int_For_P(file);

	size = (*row)*(*col);

	buffer = new byte[size];
	handle = buffer;

	if (buffer == NULL)
	{
		std::cerr << "error: out of memory in ReadPgm function" << std::endl;
		exit(1);
	}

	// read image from file
	file.read((char*)buffer, size);
	file.close();
	return handle;
}


/* FOR PGM, PPM, AND PBM FILE FORMAT */
INT Read_Char_For_P(std::ifstream& file)
{
	INT c;

	c = file.get();
	if (c == '#')
	{
		do{
			c = file.get();
		} while (c != '\n' && c != file.eof());
	}

	return c;
}
INT Read_Int_For_P(std::ifstream& file)
{
	INT c, val;

	do{
		c = Read_Char_For_P(file);
	} while (c == ' ' || c == '\t' || c == '\n');

	val = c - '0';

	while ((c = Read_Char_For_P(file)) >= '0' && c <= '9')
	{
		val *= 10;
		val += c - '0';
	}

	return val;
}
/* END FOR PGM, PPM, AND PBM FILE FORMAT */

VOID WriteBmp(INT row, INT col, BYTE* img, LPCSTR  filename)
{
	register INT i, j;
	std::ofstream file;
	INT m_col;
	INT tempInt;
	WORD tempWord;
	WORD check = 19778; // "BM"
	BYTE *rawdata;

	file.open(filename, std::ios::out | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in WriteBmp function" << std::endl;
		exit(1);
	}

	//type "BM" write
	file.write((char*)&check, 2);

	if ((col * 3) % 4 == 0)
	{
		m_col = col * 3;
	}
	else
	{
		m_col = 4 * ((col * 3) / 4) + 4;
	}

	tempInt = 54 + m_col*row; file.write((char*)&tempInt, 4);		//file size
	tempInt = 0; file.write((char*)&tempInt, 4);					//reserved1
	tempInt = 54; file.write((char*)&tempInt, 4);					//reserved2
	tempInt = 40; file.write((char*)&tempInt, 4);					//offset

	//col & row write
	file.write((char*)&col, 4);
	file.write((char*)&row, 4);

	tempWord = 1; file.write((char*)&tempWord, 2);					//color plane

	// num of bit(24 bit)
	tempWord = 24; file.write((char*)tempWord, 2);

	tempInt = 0; file.write((char*)&tempInt, 4);					//compresstion flag
	tempInt = m_col*row; file.write((char*)&tempInt, 4);			//rawdata size
	tempInt = 3937; file.write((char*)&tempInt, 4);					//width resolution
	tempInt = 3937; file.write((char*)&tempInt, 4);					//height resolution

	// load color num
	tempInt = 0; file.write((char*)&tempInt, 4);

	tempInt = 0; file.write((char*)&tempInt, 4);					//important color

	//load file
	rawdata = new BYTE[row*col * 3];

	//raw data(rgb plane) to bmp format
	for (i = 0; i < row; ++i) {
		for (j = 0; j < col; ++j) {
			rawdata[3 * ((row - i - 1)*col + j) + 2] = img[i*col + j];
			rawdata[3 * ((row - i - 1)*col + j) + 1] = img[i*col + j + col*row];
			rawdata[3 * ((row - i - 1)*col + j)] = img[i*col + j + col*row * 2];
		}
	}

	//rgb0 rgb0 ... 
	for (i = 0; i < row; ++i) {
		file.write((char*)(rawdata + i * col * 3), col * 3);
		if (m_col != (col * 3)) {
			for (j = 0; j < (m_col - col * 3); j++)	file.put(0);
		}
	}

	delete[]rawdata;
	file.close();
}

BYTE* ReadBmp(INT *row, INT *col, LPCSTR  filename)
{
	register INT i, j;
	std::ifstream file;
	BmpInforHearder Bmp;
	WORD check;
	INT tempInt;
	INT Size;
	INT tempImg;
	BYTE *rawdata, *reverse;
	BYTE tempByte;

	file.open(filename, std::ios::in | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in ReadBmp function" << std::endl;
		exit(1);
	}

	//check 'BM' type
	file.read((char*)&check, 2);
	if (check != 19778)
	{
		std::cerr << "error : image format is not a bmp" << std::endl;
		exit(1);
	}
	/* // alter
	if ((file.get() != 'B') || (file.get() != 'M'))
	{
	cerr << "errer : image format is not a bmp" << endl;
	exit(1);
	}*/

	file.read((char*)&tempInt, 4);					// file size
	file.read((char*)&tempInt, 4);					// reserved1
	file.read((char*)&tempInt, 4);					// reserved2
	file.read((char*)&tempInt, 4);					// offset

	//width & height
	file.read((char*)&Bmp.Width, 4);
	file.read((char*)&Bmp.Height, 4);
	*col = Bmp.Width;
	*row = Bmp.Height;

	file.read((char*)&tempInt, 2);					// colorplane is always 1

	//bit per pixel
	file.read((char*)&Bmp.BitCount, 2);

	file.read((char*)&tempInt, 4);					// compresstion flag
	file.read((char*)&tempInt, 4);					// raw data size
	file.read((char*)&tempInt, 4);					// width resolution
	file.read((char*)&tempInt, 4);					// height resolution

	//color number
	file.read((char*)&Bmp.NumOfColor, 4);
	Bmp.NumOfColor = Bmp.NumOfColor;

	file.read((char*)&tempInt, 4); // important color

	Size = Bmp.Width * Bmp.Height * 3;
	rawdata = new BYTE[Size];
	reverse = new BYTE[Size];

	if (Bmp.NumOfColor != 0)
	{
		std::cout << "not 24 bit" << std::endl;
		exit(1);
	}
	else
	{
		for (i = 0; i < Size * 3; i += 3)
		{
			if (((Bmp.Width % 4) != 0) && ((Bmp.Height % 3) == 0) && (i != 0))
			{
				for (j = 0; j < (4 - (Bmp.Width * 3) % 4); ++j)
				{
					file.read((char*)&tempByte, 1);
				}
			}
			tempImg = ((Bmp.Height - ((i / 3) / Bmp.Width) - 1)*Bmp.Width + (i / 3) % Bmp.Width) * 3;

			file.read((char*)(rawdata + tempImg + 2), 1);
			file.read((char*)(rawdata + tempImg + 1), 1);
			file.read((char*)(rawdata + tempImg), 1);
		}
	}
	file.close();

	memcpy(reverse, rawdata, Size);
	//std::memcpy(reverse, rawdata, Size);

	for (i = 0; i < Bmp.Height; ++i)
	{
		for (j = 0; j < Bmp.Width; ++j)
		{
			rawdata[i*Bmp.Width + j] = reverse[3 * (i*Bmp.Width + j)];
			rawdata[i*Bmp.Width + j + Bmp.Width] = reverse[3 * (i*Bmp.Width + j) + 1];
			rawdata[i*Bmp.Width + j + Bmp.Width * 2] = reverse[3 * (i*Bmp.Width + j) + 2];
		}
	}

	delete[] reverse;
	return rawdata;
}

VOID WritePpm(INT row, INT col, BYTE* img, LPCSTR  filename, INT order)
{
	std::ofstream file;
	register INT i;
	INT size;
	BYTE *buffer = NULL, *rgb, *handle = NULL;
	INT quarter = 4, round, qSize;
	char* header = "P6\n# Created by Vision Image Processing(VIP) Lab, Sangmyung Univ. Korea \n";

	file.open(filename, std::ios::out | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in WritePgm function" << std::endl;
		exit(1);
	}

	// planar
	else if (order == PLANAR)
	{
		size = row * col;
		qSize = size / quarter;
		round = size - quarter * qSize;

		rgb = img;
		buffer = new BYTE[size * 3];
		handle = buffer;

		if (buffer == NULL)
		{
			std::cerr << "error: out of memory in WritePpm function" << std::endl;
			exit(1);
		}

		for (i = 0; i < qSize; ++i)
		{
			handle[0] = *(rgb);
			handle[1] = *(rgb + size);
			handle[2] = *(rgb + 2 * size);
			rgb++;

			handle[3] = *(rgb);
			handle[4] = *(rgb + size);
			handle[5] = *(rgb + 2 * size);
			rgb++;

			handle[6] = *(rgb);
			handle[7] = *(rgb + size);
			handle[8] = *(rgb + 2 * size);
			rgb++;

			handle[9] = *(rgb);
			handle[10] = *(rgb + size);
			handle[11] = *(rgb + 2 * size);
			rgb++;

			handle += 12;
		}

		for (i = 0; i < round; ++i)
		{
			*handle = *(rgb); handle++;
			*handle = *(rgb + size); handle++;
			*handle = *(rgb + 2 * size); handle++;
			rgb++;
		}
	}

	//header information writing phase 'P6' is for ppm file format
	file << header << col << ' ' << row << std::endl << "255" << std::endl;

	if (order == INTERLEAVED)
	{
		file.write((char*)img, row*col * 3);
	}
	else if (order == PLANAR)
	{
		file.write((char*)buffer, row*col * 3);
	}
	file.close();
}
BYTE* ReadPpm(INT *row, INT *col, LPCSTR  filename, INT order)
{
	std::ifstream file;
	register INT i;
	INT size;
	BYTE *buffer = NULL, *rgb, *handle;

	file.open(filename, std::ios::in | std::ios::binary);
	if (file.fail())
	{
		std::cout << "file open error in ReadPpm function" << std::endl;
		exit(1);
	}
	//header information reading phase 'P5' is for pgm file format
	if ((file.get() != 'P') || (file.get() != '6'))
	{
		std::cerr << "error : image format is not a ppm" << std::endl;
		exit(1);
	}
	//read col & row
	*col = Read_Int_For_P(file);
	*row = Read_Int_For_P(file);
	Read_Int_For_P(file);

	size = (*row)*(*col);

	rgb = new byte[size * 3];
	handle = rgb;

	buffer = new byte[size * 3];
	if (buffer == NULL)
	{
		std::cerr << "error: out of memory in ReadPpm function" << std::endl;
		exit(1);
	}

	// read image from file
	file.read((char*)buffer, size * 3);
	std::cout << "444444" << std::endl;
	if (order == PLANAR)
	{
		// interleaved RGB order into planar RGB order
		for (i = 0; i < size; ++i)
		{
			rgb[0] = buffer[0];				//R
			rgb[size] = buffer[1];			//G
			rgb[2 * size] = buffer[2];		//B
			rgb++;
			buffer += 3;
		}

		buffer -= size;
		buffer -= size;
		buffer -= size;

		delete[]buffer;
	}

	file.close();
	return handle;
}




template <class T>
DOUBLE gPSNR(T* origin, T* target, CONST INT length, CONST INT max)
{
	register INT i;
	DOUBLE PSNR = 0.0;
	DOUBLE MSE = 0.0;
	DOUBLE SUM = 0.0;

	for (i = 0; i < length; ++i)
	{
		SUM += (origin[i] - target[i]) * (origin[i] - target[i]);
	}

	MSE = SUM / length;

	PSNR = 20.0 * log10(max / sqrt(MSE));

	return PSNR;
}

template <class T>
DOUBLE gPSNR(T* origin, T* target, CONST INT height, CONST INT width, CONST INT max, CONST INT boundary)
{
	register INT i, j;
	DOUBLE PSNR = 0.0;
	DOUBLE MSE = 0.0;
	DOUBLE SUM = 0.0;
	LONG LENGTH = 0;

	for (i = boundary; i < height - boundary; ++i)
	{
		for (j = boundary; j < width - boundary; ++j)
		{
			SUM += (origin[i*width + j] - target[i*width + j])*(origin[i*width + j] - target[i*width + j]);
			LENGTH++;
		}
	}

	MSE = SUM / LENGTH;

	PSNR = 20.0 * log10(max / sqrt(MSE));

	return PSNR;
}

void boxcar_sliding_window_x(double *out, double *in, int h, int w, int radius)
{
	double scale = 1.0f / (2 * radius + 1);
	for (int y = 0; y < h; y++) {
		double t;
		// do left edge
		t = in[y*w] * radius;
		for (int x = 0; x < radius + 1; x++) {
			t += in[y*w + x];
		}
		out[y*w] = t * scale;
		for (int x = 1; x < radius + 1; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[y*w];
			out[c] = t * scale;
		}
		// main loop
		for (int x = radius + 1; x < w - radius; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}
		// do right edge
		for (int x = w - radius; x < w; x++) {
			int c = y*w + x;
			t += in[(y*w) + w - 1];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}

	}
}
void boxcar_sliding_window_y(double *out, double *in, int h, int w, int radius)
{
	double scale = 1.0f / (2 * radius + 1);
	for (int x = 0; x < w; x++)
	{
		double t;
		// do left edge
		t = in[x] * radius;
		for (int y = 0; y < radius + 1; y++) {
			t += in[y*w + x];
		}
		out[x] = t * scale;
		for (int y = 1; y < radius + 1; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[x];
			out[c] = t * scale;
		}
		// main loop
		for (int y = radius + 1; y < h - radius; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
		// do right edge
		for (int y = h - radius; y < h; y++) {
			int c = y*w + x;
			t += in[(h - 1)*w + x];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
	}
}
void boxcar_sliding_window(double **out, double **in, double **temp, int h, int w, int radius)
{
	boxcar_sliding_window_x(temp[0], in[0], h, w, radius);
	boxcar_sliding_window_y(out[0], temp[0], h, w, radius);
}

void boxcar_sliding_window_x(float *out, float *in, int h, int w, int radius)
{
	float scale = 1.0f / (2 * radius + 1);
	for (int y = 0; y < h; y++) {
		float t;
		// do left edge
		t = in[y*w] * radius;
		for (int x = 0; x < radius + 1; x++) {
			t += in[y*w + x];
		}
		out[y*w] = t * scale;
		for (int x = 1; x < radius + 1; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[y*w];
			out[c] = t * scale;
		}
		// main loop
		for (int x = radius + 1; x < w - radius; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}
		// do right edge
		for (int x = w - radius; x < w; x++) {
			int c = y*w + x;
			t += in[(y*w) + w - 1];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}

	}
}
void boxcar_sliding_window_y(float *out, float *in, int h, int w, int radius)
{
	float scale = 1.0f / (2 * radius + 1);
	for (int x = 0; x < w; x++)
	{
		float t;
		// do left edge
		t = in[x] * radius;
		for (int y = 0; y < radius + 1; y++) {
			t += in[y*w + x];
		}
		out[x] = t * scale;
		for (int y = 1; y < radius + 1; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[x];
			out[c] = t * scale;
		}
		// main loop
		for (int y = radius + 1; y < h - radius; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
		// do right edge
		for (int y = h - radius; y < h; y++) {
			int c = y*w + x;
			t += in[(h - 1)*w + x];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
	}
}
void boxcar_sliding_window(float**out, float**in, float**temp, int h, int w, int radius)
{
	int min_hw = min(h, w);
	if (radius >= min_hw)
	{
		double dsum = 0;
		float*in_ = in[0];
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) dsum += *in_++;
		dsum /= (h*w);
		float fsum = (float)dsum;
		float*out_ = out[0];
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) *out_++ = fsum;
	}
	else if (radius>0)
	{
		boxcar_sliding_window_x(temp[0], in[0], h, w, radius);
		boxcar_sliding_window_y(out[0], temp[0], h, w, radius);
	}
	else memcpy(out[0], in[0], sizeof(float)*h*w);
}
void boxcar_sliding_window_x(byte*out, byte*in, int h, int w, int radius)
{
	float scale = 1.0f / (2 * radius + 1);
	for (int y = 0; y < h; y++) {
		float t;
		// do left edge
		t = in[y*w] * radius;
		for (int x = 0; x < radius + 1; x++) {
			t += in[y*w + x];
		}
		out[y*w] = byte(t * scale + 0.5);
		for (int x = 1; x < radius + 1; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[y*w];
			out[c] = byte(t * scale + 0.5);
		}
		// main loop
		for (int x = radius + 1; x < w - radius; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[c - radius - 1];
			out[c] = byte(t * scale + 0.5);
		}
		// do right edge
		for (int x = w - radius; x < w; x++) {
			int c = y*w + x;
			t += in[(y*w) + w - 1];
			t -= in[c - radius - 1];
			out[c] = byte(t * scale + 0.5);
		}

	}
}
void boxcar_sliding_window_y(byte*out, byte*in, int h, int w, int radius)
{
	float scale = 1.0f / (2 * radius + 1);
	for (int x = 0; x < w; x++)
	{
		float t;
		// do left edge
		t = in[x] * radius;
		for (int y = 0; y < radius + 1; y++) {
			t += in[y*w + x];
		}
		out[x] = byte(t * scale + 0.5);
		for (int y = 1; y < radius + 1; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[x];
			out[c] = byte(t * scale + 0.5);
		}
		// main loop
		for (int y = radius + 1; y < h - radius; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[c - (radius*w) - w];
			out[c] = byte(t * scale + 0.5);
		}
		// do right edge
		for (int y = h - radius; y < h; y++) {
			int c = y*w + x;
			t += in[(h - 1)*w + x];
			t -= in[c - (radius*w) - w];
			out[c] = byte(t * scale + 0.5);
		}
	}
}
void boxcar_sliding_window(byte**out, byte**in, byte**temp, int h, int w, int radius)
{
	int min_hw = min(h, w);
	if (radius<min_hw)
	{
		boxcar_sliding_window_x(temp[0], in[0], h, w, radius);
		boxcar_sliding_window_y(out[0], temp[0], h, w, radius);
	}
	else
	{
		double dsum = 0;
		byte*in_ = in[0];
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) dsum += *in_++;
		dsum /= (h*w);
		byte usum = byte(dsum + 0.5);
		byte*out_ = out[0];
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) *out_++ = usum;
	}
}
void gaussian_recursive_x(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn)
{
	double xp = 0.0f;  // previous input
	double yp = 0.0f;  // previous output
	double yb = 0.0f;  // previous output by 2
	for (int y = 0; y<h; y++)
	{
		xp = id[y][0]; yb = coefp*xp; yp = yb;
		for (int x = 0; x < w; x++)
		{
			double xc = id[y][x];
			double yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}
	// reverse pass
	// ensures response is symmetrical
	double xn = 0.f;
	double xa = 0.f;
	double yn = 0.f;
	double ya = 0.f;
	for (int y = 0; y<h; y++)
	{
		xn = xa = id[y][w - 1]; yn = coefn*xn; ya = yn;
		for (int x = w - 1; x >= 0; x--) {
			double xc = id[y][x];
			double yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
void gaussian_recursive_y(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn)
{
	double xp = 0.0f;  // previous input
	double yp = 0.0f;  // previous output
	double yb = 0.0f;  // previous output by 2
	for (int x = 0; x < w; x++)
	{
		xp = id[0][x]; yb = coefp*xp; yp = yb;
		for (int y = 0; y<h; y++)
		{
			double xc = id[y][x];
			double yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}


	// reverse pass
	// ensures response is symmetrical
	double xn = 0.f;
	double xa = 0.f;
	double yn = 0.f;
	double ya = 0.f;
	for (int x = 0; x < w; x++)
	{
		xn = xa = id[h - 1][x]; yn = coefn*xn; ya = yn;
		for (int y = h - 1; y >= 0; y--)
		{
			double xc = id[y][x];
			double yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
int gaussian_recursive(double **image, double **temp, double sigma, int order, int h, int w)
{
	const double
		nsigma = sigma < 0.1f ? 0.1f : sigma,
		alpha = 1.695f / nsigma,
		ema = exp(-alpha),
		ema2 = exp(-2 * alpha),
		b1 = -2 * ema,
		b2 = ema2;
	double a0 = 0, a1 = 0, a2 = 0, a3 = 0, coefp = 0, coefn = 0;
	switch (order) {
	case 0: {
				const double k = (1 - ema)*(1 - ema) / (1 + 2 * alpha*ema - ema2);
				a0 = k;
				a1 = k*(alpha - 1)*ema;
				a2 = k*(alpha + 1)*ema;
				a3 = -k*ema2;
	} break;

	case 1: {
				const double k = (1 - ema)*(1 - ema) / ema;
				a0 = k*ema;
				a1 = a3 = 0;
				a2 = -a0;
	} break;

	case 2: {
				const double
					ea = exp(-alpha),
					k = -(ema2 - 1) / (2 * alpha*ema),
					kn = (-2 * (-1 + 3 * ea - 3 * ea*ea + ea*ea*ea) / (3 * ea + 1 + 3 * ea*ea + ea*ea*ea));
				a0 = kn;
				a1 = -kn*(1 + k*alpha)*ema;
				a2 = kn*(1 - k*alpha)*ema;
				a3 = -kn*ema2;
	} break;

	default:
		fprintf(stderr, "gaussianFilter: invalid order parameter!\n");
		return 0;
	}
	coefp = (a0 + a1) / (1 + b1 + b2);
	coefn = (a2 + a3) / (1 + b1 + b2);
	//timer.start();
	gaussian_recursive_x(temp, image, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	//image_display(temp,h,w);
	gaussian_recursive_y(image, temp, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	//timer.fps_display();
	return(0);
}
void gaussian_recursive_x(float **od, float **id, int w, int h, float a0, float a1, float a2, float a3, float b1, float b2, float coefp, float coefn)
{
	//image_display(id,h,w);
	float xp = 0.0f;  // previous input
	float yp = 0.0f;  // previous output
	float yb = 0.0f;  // previous output by 2
	for (int y = 0; y<h; y++)
	{
		xp = id[y][0]; yb = coefp*xp; yp = yb;
		for (int x = 0; x < w; x++)
		{
			float xc = id[y][x];
			float yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}

	//image_display(od,h,w);
	// reverse pass
	// ensures response is symmetrical
	float xn = 0.f;
	float xa = 0.f;
	float yn = 0.f;
	float ya = 0.f;
	for (int y = 0; y<h; y++)
	{
		xn = xa = id[y][w - 1]; yn = coefn*xn; ya = yn;
		for (int x = w - 1; x >= 0; x--) {
			float xc = id[y][x];
			float yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
void gaussian_recursive_y(float **od, float **id, int w, int h, float a0, float a1, float a2, float a3, float b1, float b2, float coefp, float coefn)
{
	float xp = 0.0f;  // previous input
	float yp = 0.0f;  // previous output
	float yb = 0.0f;  // previous output by 2
	for (int x = 0; x < w; x++)
	{
		xp = id[0][x]; yb = coefp*xp; yp = yb;
		for (int y = 0; y<h; y++)
		{
			float xc = id[y][x];
			float yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}


	// reverse pass
	// ensures response is symmetrical
	float xn = 0.f;
	float xa = 0.f;
	float yn = 0.f;
	float ya = 0.f;
	for (int x = 0; x < w; x++)
	{
		xn = xa = id[h - 1][x]; yn = coefn*xn; ya = yn;
		for (int y = h - 1; y >= 0; y--)
		{
			float xc = id[y][x];
			float yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
int gaussian_recursive(float **image, float **temp, float sigma, int order, int h, int w)
{
	const float
		nsigma = sigma < 0.1f ? 0.1f : sigma,
		alpha = 1.695f / nsigma,
		ema = exp(-alpha),
		ema2 = exp(-2 * alpha),
		b1 = -2 * ema,
		b2 = ema2;
	float a0 = 0, a1 = 0, a2 = 0, a3 = 0, coefp = 0, coefn = 0;
	switch (order) {
	case 0: {
				const float k = (1 - ema)*(1 - ema) / (1 + 2 * alpha*ema - ema2);
				a0 = k;
				a1 = k*(alpha - 1)*ema;
				a2 = k*(alpha + 1)*ema;
				a3 = -k*ema2;
	} break;

	case 1: {
				const float k = (1 - ema)*(1 - ema) / ema;
				a0 = k*ema;
				a1 = a3 = 0;
				a2 = -a0;
	} break;

	case 2: {
				const float
					ea = exp(-alpha),
					k = -(ema2 - 1) / (2 * alpha*ema),
					kn = (-2 * (-1 + 3 * ea - 3 * ea*ea + ea*ea*ea) / (3 * ea + 1 + 3 * ea*ea + ea*ea*ea));
				a0 = kn;
				a1 = -kn*(1 + k*alpha)*ema;
				a2 = kn*(1 - k*alpha)*ema;
				a3 = -kn*ema2;
	} break;

	default:
		fprintf(stderr, "gaussianFilter: invalid order parameter!\n");
		return 0;
	}
	coefp = (a0 + a1) / (1 + b1 + b2);
	coefn = (a2 + a3) / (1 + b1 + b2);
	//timer.start();
	gaussian_recursive_x(temp, image, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	gaussian_recursive_y(image, temp, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	//timer.fps_display();
}


void qx_specular_free_image(byte ***image_specular_free, byte ***image_normalized, float **diffuse_chromaticity_max, int h, int w)
{
	int y, x;
	byte *image_specular_free_x, *image_normalized_x; float *diffuse_chromaticity_max_x;
	byte r, g, b; double imax, isum; float rf, gf, bf, c, t0, t1, t2, t3, diffuse, specular;
	//*image_sum_x,*image_max_x,*chromaticity_max_x,
	image_specular_free_x = image_specular_free[0][0];
	image_normalized_x = image_normalized[0][0];
	diffuse_chromaticity_max_x = diffuse_chromaticity_max[0];
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			t1 = 3.f*(*diffuse_chromaticity_max_x++) - 1.f;
			t3 = 1.0f / 3.0f;
			r = (*image_normalized_x++);
			g = (*image_normalized_x++);
			b = (*image_normalized_x++);
			if (t1>0)
			{
				isum = r + g + b;
				if (isum == 0) c = 0;
				else
				{
					imax = max(max(r, g), b);
					c = (float)(imax / isum);
				}
				t0 = t1*c;
				if (fabs(t0)<QX_DEF_THRESHOLD_ZERO)
				{
					*image_specular_free_x++ = r;
					*image_specular_free_x++ = g;
					*image_specular_free_x++ = b;
				}
				else
				{
					t2 = (3.0f*c - 1.f);
					diffuse = float(imax*t2 / t0);
					specular = float(t3*(isum - diffuse));
					rf = r - specular;
					gf = g - specular;
					bf = b - specular;
					if (rf<0.f) rf = 0.f; else if (rf>255.f) rf = 255.f;
					if (gf<0.f) gf = 0.f; else if (gf>255.f) gf = 255.f;
					if (bf<0.f) bf = 0.f; else if (bf>255.f) bf = 255.f;
					*image_specular_free_x++ = byte(rf + 0.5f);
					*image_specular_free_x++ = byte(gf + 0.5f);
					*image_specular_free_x++ = byte(bf + 0.5f);
				}
			}
			else
			{
				*image_specular_free_x++ = r;
				*image_specular_free_x++ = g;
				*image_specular_free_x++ = b;
			}
		}
	}
}

void qx_stereo_flip_corr_vol(double***corr_vol_right, double***corr_vol, int h, int w, int nr_plane)
{
	for (int y = 0; y<h; y++)
	{
		for (int x = 0; x<w - nr_plane; x++) for (int d = 0; d<nr_plane; d++) corr_vol_right[y][x][d] = corr_vol[y][x + d][d];
		for (int x = w - nr_plane; x<w; x++) for (int d = 0; d<nr_plane; d++)
		{
			if ((x + d)<w) corr_vol_right[y][x][d] = corr_vol[y][x + d][d];
			else corr_vol_right[y][x][d] = corr_vol_right[y][x][d - 1];
		}
	}
}
void depth_best_cost(byte**depth, double***evidence, int h, int w, int nr_planes)
{
	for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) { int d; vec_min_pos(d, evidence[y][x], nr_planes); depth[y][x] = d; }
}
void vec_min_pos(int &min_pos, double *in, int len)
{
	double min_val = in[0];
	min_pos = 0;
	for (int i = 1; i<len; i++) if (in[i]<min_val)
	{
		min_val = in[i];
		min_pos = i;
	}
}
void qx_detect_occlusion_left_right(byte**mask_left, byte**depth_left, byte**depth_right, int h, int w, int nr_plane)
{
	memset(mask_left[0], 0, sizeof(char)*h*w);
	for (int y = 0; y<h; y++)
	{
		for (int x = 0; x<w; x++)
		{
			int d = depth_left[y][x];
			int xr = x - d;
			if (xr >= 0)
			{
				if (d == 0 || abs(d - depth_right[y][xr]) >= 1)
				{
					//depth_left[y][x]=min(depth_left[y][x],depth_right[y][xr]);
					mask_left[y][x] = 255;
				}
			}
			else mask_left[y][x] = 255;
		}
	}
}

int file_open_ascii(char *file_path, int *out, int len)
{
	FILE *file_in; char str[65]; int i;
	file_in=fopen(file_path,"r");
	//fopen(&file_in, file_path, "r");
	if (file_in != NULL)
	{
		fseek(file_in, 0, SEEK_SET);
		for (i = 0; i<len; i++)
		{
			//fscanf(file_in,"%s",str); 
			fscanf(file_in, "%s", str, 65);
			out[i] = atoi(str);
		}
		fclose(file_in);
	}
	else
	{
		printf("qx_basic_file: Can not open file: %s\n", file_path);
		getchar();
		exit(-1);
	}
	return(0);
}

inline int qx_mst_yx_2_image_index(int y, int x, int h, int w){ return(y*w + x); }
inline int qx_mst_compute_nr_edge_4neighbor(int h, int w){ if (h <= 2 && w <= 2) return(0); else return((h - 1)*w + (w - 1)*h); }
inline int qx_mst_compute_nr_edge_8neighbor(int h, int w){ if (h <= 2 && w <= 2) return(0); else return((h - 1)*w + (w - 1)*h + (h - 1)*(w - 1) * 2); }
inline void qx_mst_compute_edges_per_pixel(int**edges, byte*distance, byte*image, int nr_channel, int &nr_edge, int y0, int x0, int yt, int xt, int h, int w)
{
	int id0 = edges[nr_edge][0] = qx_mst_yx_2_image_index(y0, x0, h, w);
	int idt = edges[nr_edge][1] = qx_mst_yx_2_image_index(yt, xt, h, w);
	if (nr_channel == 1)
		distance[nr_edge] = abs(image[idt] - image[id0]);
	//edges[nr_edge++][2]=min(10,abs(image[edges[nr_edge][1]]-image[edges[nr_edge][0]])); 
	else if (nr_channel == 3)
	{
		id0 *= nr_channel;
		idt *= nr_channel;
		//distance[nr_edge]=euro_dist_rgb_max(&(image[idt]),&(image[id0]));


		byte r = abs(image[idt++] - image[id0++]);
		byte g = abs(image[idt++] - image[id0++]);
		byte b = abs(image[idt++] - image[id0++]);
		distance[nr_edge] = qx_max_u3(r, g, b);
		//distance[nr_edge]=int((r+g+b)*0.33+0.5);
		//double e=(r+g+b)*0.5;//euro_dist_rgb_mean(&(image[edges[nr_edge][1]*3]),&(image[edges[nr_edge][0]*3]));
		//edges[nr_edge++][2]=min(10,euro_dist_rgb_max(&(image[edges[nr_edge][1]*3]),&(image[edges[nr_edge][0]*3])));
		//if(e<QX_DEF_MST_KI_MIN_DISTANT) distance[nr_edge]=0;
		//else distance[nr_edge]=min(255,int(e+0.5));

	}
	else
	{
		id0 *= nr_channel;
		idt *= nr_channel;
		int cost_max = 0;
		for (int i = 0; i<nr_channel; i++)
		{
			cost_max = max(cost_max, abs(image[idt + i] - image[id0 + i]));
		}
		distance[nr_edge] = cost_max;
	}
	nr_edge++;
}
inline int qx_mst_compute_edges_4neighbor(int**edges, byte*distance, byte*image, int nr_channel, int h, int w)
{
	int y0, yt, x0, xt;
	int nr_edge = 0;
	for (y0 = 0; y0<h; y0++)
	{
		yt = y0;
		for (int x0 = 0; x0<w - 1; x0++)
		{
			xt = x0 + 1;
			qx_mst_compute_edges_per_pixel(edges, distance, image, nr_channel, nr_edge, y0, x0, yt, xt, h, w);
		}
	}
	for (int x0 = 0; x0<w; x0++)
	{
		xt = x0;
		for (y0 = 0; y0<h - 1; y0++)
		{
			yt = y0 + 1;
			qx_mst_compute_edges_per_pixel(edges, distance, image, nr_channel, nr_edge, y0, x0, yt, xt, h, w);
		}
	}
	return(nr_edge);
}

inline int qx_mst_compute_edges_8neighbor(int**edges, byte*distance, byte*image, int nr_channel, int h, int w)
{
	int y0, yt, x0, xt;
	int nr_edge = qx_mst_compute_edges_4neighbor(edges, distance, image, nr_channel, h, w);
	for (y0 = 0; y0<h - 1; y0++)
	{
		yt = y0 + 1;
		for (int x0 = 0; x0<w - 1; x0++)
		{
			xt = x0 + 1;
			qx_mst_compute_edges_per_pixel(edges, distance, image, nr_channel, nr_edge, y0, x0, yt, xt, h, w);
			qx_mst_compute_edges_per_pixel(edges, distance, image, nr_channel, nr_edge, y0, xt, yt, x0, h, w);
		}
	}
	return(nr_edge);
}


qx_mst_kruskals_image::qx_mst_kruskals_image()
{
	m_image = NULL;
	m_edge = NULL;
	m_distance = NULL;
	m_id_edge = NULL;

	m_parent = NULL;
	m_nr_child = NULL;
	m_children = NULL;
	m_weight = NULL;


	m_connected = NULL;
	m_connected_distance = NULL;
	m_nr_connected = NULL;

	m_node_id_from_parent_to_child = NULL;
	m_rank = NULL;

	m_parent_default = QX_DEF_MST_KI_PARENT_DEFAULT;
}
qx_mst_kruskals_image::~qx_mst_kruskals_image()
{
	clean();
}
void qx_mst_kruskals_image::clean()
{
	qx_freeu_1(m_image); m_image = NULL;
	qx_freei(m_edge); m_edge = NULL;
	qx_freeu_1(m_distance); m_distance = NULL;
	qx_freei_1(m_id_edge); m_id_edge = NULL;

	qx_freei_1(m_parent); m_parent = NULL;
	qx_freei_1(m_nr_child); m_nr_child = NULL;
	qx_freeu_1(m_weight); m_weight = NULL;
	qx_freei(m_children); m_children = NULL;

	qx_freei(m_connected); m_connected = NULL;
	qx_freeu(m_connected_distance); m_connected_distance = NULL;
	qx_freei_1(m_nr_connected); m_nr_connected = NULL;

	qx_freei_1(m_node_id_from_parent_to_child); m_node_id_from_parent_to_child = NULL;
	qx_freei_1(m_rank); m_rank = NULL;
}
int qx_mst_kruskals_image::init(int h, int w, int nr_channel, int nr_neighbor)
{
	clean();
	m_h = h; m_w = w; m_nr_channel = nr_channel; m_nr_neighbor = nr_neighbor;
	m_max_nr_child = m_nr_neighbor - 1;
	m_nr_vertices = m_h*m_w;
	if (m_nr_neighbor == QX_DEF_MST_KI_4NR_NEIGHBOR) m_nr_edge = qx_mst_compute_nr_edge_4neighbor(m_h, m_w);
	else m_nr_edge = qx_mst_compute_nr_edge_8neighbor(m_h, m_w);
	m_image = new byte[m_h*m_w*m_nr_channel];
	m_edge = qx_alloci(m_nr_edge, 2);
	m_distance = new byte[m_nr_edge];
	m_id_edge = new int[m_nr_edge];

	m_parent = new int[m_nr_vertices];
	m_nr_child = new int[m_nr_vertices];
	m_weight = new byte[m_nr_vertices];
	m_children = qx_alloci(m_nr_vertices, m_max_nr_child);

	m_connected = qx_alloci(m_nr_vertices, m_nr_neighbor);
	m_connected_distance = qx_allocu(m_nr_vertices, m_nr_neighbor);
	m_nr_connected = new int[m_nr_vertices];

	m_node_id_from_parent_to_child = new int[m_nr_vertices];//build tree
	m_rank = new int[m_nr_vertices];
	m_queue.init(m_nr_vertices);
	return(0);
}
void qx_mst_kruskals_image::init_mst()
{
	int*parent = m_parent;
	for (int i = 0; i<m_nr_vertices; i++) *parent++ = i;
	memset(m_nr_connected, 0, sizeof(int)*m_nr_vertices);
}
int qx_mst_kruskals_image::mst(byte*image, bool print_edges)
{
	
	init_mst();
	//timer.fps_display();
	
	ctmf(image, m_image, m_w, m_h, m_w*m_nr_channel, m_w*m_nr_channel, 1, m_nr_channel, m_h*m_w*m_nr_channel);
	if (m_nr_neighbor == QX_DEF_MST_KI_4NR_NEIGHBOR) qx_mst_compute_edges_4neighbor(m_edge, m_distance, m_image, m_nr_channel, m_h, m_w);//find edges
	else qx_mst_compute_edges_8neighbor(m_edge, m_distance, m_image, m_nr_channel, m_h, m_w);
	//timer.fps_display();
	
	qx_sort_increase_using_histogram(m_id_edge, m_distance, m_nr_edge);//sorting edges in a nondesc
	//timer.fps_display("qx_sort_increase_using_histogram");
	
	kruskal();
	//timer.fps_display("kruskal");
	build_tree();
	return(0);
}

int qx_mst_kruskals_image::findset(int x)
{
	//if(x==m_nr_vertices) printf("[%d - %d]",x,m_nr_vertices);
	int parent = m_parent[x];
	if (x != parent)
	{
		m_parent[x] = findset(parent);
	}
	return m_parent[x];
}
void qx_mst_kruskals_image::kruskal()
{
	m_tree_size = 0;
	for (int j = 0; j<m_nr_edge; j++)
	{
		int i = m_id_edge[j];
		int*edge = m_edge[i];
		int u = edge[0];
		int v = edge[1];
		int pu = findset(u);
		int pv = findset(v);
		if (pu != pv)
		{
			int nr_connected = m_nr_connected[u];
			m_connected[u][nr_connected] = v;
			m_connected_distance[u][nr_connected] = m_distance[i];
			m_nr_connected[u]++;

			nr_connected = m_nr_connected[v];
			m_connected[v][nr_connected] = u;
			m_connected_distance[v][nr_connected] = m_distance[i];
			m_nr_connected[v]++;

			m_tree_size++;
			//printf("( %d, %d ): %d\n", u,v,m_distance[i]);
			m_parent[pu] = m_parent[pv]; // link
		}
	}
	//printf("m_total_weight: %d\n",total_weight);
	//printf("size: [%d - %d]\n",m_tree_size,m_h*m_w);
}
int qx_mst_kruskals_image::build_tree()
{

	int tree_parent = 0;
	int parent = tree_parent;
	memset(m_parent, m_parent_default, sizeof(int)*m_nr_vertices);
	memset(m_nr_child, 0, sizeof(int)*m_nr_vertices);
	memset(m_rank, 0, sizeof(int)*m_nr_vertices);

	m_parent[parent] = parent;
	m_weight[parent] = 0;
	m_node_id_from_parent_to_child[0] = parent;
	int len = 1;
	m_queue.reinit();
	m_queue.push(parent);
	//m_max_rank=0; 
	while (m_queue.length>0)
	{
		parent = m_queue.pull();
		int nr_connected = m_nr_connected[parent];
		for (int i = 0; i<nr_connected; i++)
		{
			int potential_child = m_connected[parent][i];
			if (m_parent[potential_child] == m_parent_default)//&&m_parent[potential_child]!=tree_parent)
			{
				m_queue.push(potential_child);
				m_parent[potential_child] = parent;
				//printf("[%d <- %d]\n",parent,potential_child);
				m_rank[potential_child] = m_rank[parent] + 1;
				byte weight = m_connected_distance[parent][i];
				m_weight[potential_child] = weight;
				//m_nodes[potential_child].value_to_be_filtered=m_image[potential_child];
				m_children[parent][(m_nr_child[parent]++)] = potential_child;

				m_node_id_from_parent_to_child[len] = potential_child;
				len++;
				if (len>m_nr_vertices)
				{
					printf("len>m_nr_vertices!!");
					getchar();
					exit(0);
				}
			}
		}
	}
	//timer.fps_display("build_tree");
	return(0);
}
void test_qx_mst_kruskals_image()
{
	qx_mst_kruskals_image m_mst;
	int h = 3, w = 3;
	//byte***a=loadimage_ppm_u("a.ppm",h,w);
	byte aa[3][3] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	byte**a = qx_allocu(h, w);
	for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) a[y][x] = aa[y][x];
	m_mst.init(h, w, 1);
	m_mst.mst(a[0], true);
	qx_freeu(a);
}


qx_queue_mst::qx_queue_mst()
{
	queue = NULL;
};
qx_queue_mst::~qx_queue_mst()
{
	clean();
};
void qx_queue_mst::clean()
{
	if (queue != NULL) delete[] queue; queue = NULL;
}
int qx_queue_mst::init(int len)
{
	if (len <= 0)
	{
		cout << "The length is: " << len << "!!";
		exit(0);
	}
	queue = new int[len * 2];
	first = 0;
	last = -1;
	length = 0;
	return(0);
}
void qx_queue_mst::reinit()
{
	first = 0;
	last = -1;
	length = 0;
}
void qx_queue_mst::push(int x)
{
	length++;
	queue[++last] = x;
}
int qx_queue_mst::pull()
{
	length--;
	return(queue[first++]);
}


float* get_binary(FILE* file_in, int h, int w, int channel)
{
	float *image = NULL;
	image = (float*)malloc(sizeof(float)*w*h*channel);
	memset(image, 0, sizeof(float)*w*h*channel);
	fread(image, sizeof(float), w*h*channel, file_in);
	return image;
}
void saveimage(char* file_name, float *image, int h, int w, int channel)
{
	FILE* file_out; int i; float maxx;
	file_out=fopen(file_name,"wb");
	//fopen(&file_out, file_name, "wb");
	maxx = image[0];
	for (i = 0; i<h*w*channel; i++) maxx = (maxx>image[i] ? maxx : image[i]);
	if (channel == 1) fprintf(file_out, "P7\n%d %d\n%f\n", w, h, maxx);
	else if (channel == 3) fprintf(file_out, "P8\n%d %d\n%f\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%f\n", w, h, channel, maxx);
	fwrite(image, sizeof(float), w*h*channel, file_out);
	fclose(file_out);
}
void saveimage(char* file_name, double *image_d, int h, int w, int channel)
{
	FILE* file_out; int i; float maxx;
	float *image = new float[h*w*channel];
	for (i = 0; i<h*w*channel; i++) image[i] = (float)image_d[i];
	file_out=fopen(file_name,"wb");
//	fopen(&file_out, file_name, "wb");
	maxx = image[0];
	for (i = 0; i<h*w*channel; i++) maxx = (maxx>image[i] ? maxx : image[i]);
	if (channel == 1) fprintf(file_out, "P7\n%d %d\n%f\n", w, h, maxx);
	else if (channel == 3) fprintf(file_out, "P8\n%d %d\n%f\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%f\n", w, h, channel, maxx);
	fwrite(image, sizeof(float), w*h*channel, file_out);
	fclose(file_out);
	delete[] image;
}
/*P2(PGM ASCII format), P3(is_ppm ASCII format), P5 (PGM Binary format), P6 (is_ppm Binary format), P7 (pgm Binary format)*/
float* loadimage(char* file_name, int &h, int &w, int *is_ppm)
{
	FILE * file_in;
	char line[LEN_MAX];
	int	i;
	float imax;
	float *image = NULL;
	if (is_ppm != NULL) *is_ppm = QX_DEF_IS_PGM;
	file_in=fopen(file_name,"rb");
//	fopen(&file_in, file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		//fscanf(file_in,"%d\n",&i);
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf(line, "%d %d\n",	&w,	&h);
			sscanf(line, "%d %d\n", &w, &h);
			break;
		}
	}
	//fscanf(file_in, "%f\n", &imax);
	fscanf(file_in, "%f\n", &imax);
	switch (i)
	{
	case 7:
		if (is_ppm != NULL) *is_ppm = QX_DEF_IS_PGM;
		image = get_binary(file_in, h, w, 1);
		break;
	case 8:
		if (is_ppm != NULL) *is_ppm = QX_DEF_IS_PPM;
		image = get_binary(file_in, h, w, 3);
		break;
	default:
		break;
	}
	fclose(file_in);
	return image;
}

int loadimage(char* file_name, float *image, int &h, int &w, int *nr_channel)
{
	FILE * file_in; int is_ppm; int nrc;
	char line[LEN_MAX];
	int	i; float imax;
	is_ppm = QX_DEF_IS_PGM;
	file_in=fopen(file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			sscanf(line, "%d %d\n", &w, &h);
			break;
		}
	}
	switch (i)
	{
	case 7:
		fscanf(file_in, "%f\n", &imax);
		if (nr_channel != NULL) (*nr_channel) = 1;
		memset(image, 0, sizeof(float)*h*w);
		fread(image, sizeof(float), h*w, file_in);
		is_ppm = QX_DEF_IS_PGM;
		break;
	case 8:
		fscanf(file_in, "%f\n", &imax);
		if (nr_channel != NULL) (*nr_channel) = 3;
		memset(image, 0, sizeof(float)*h*w * 3);
		fread(image, sizeof(float), h*w * 3, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	case 9:
		fscanf(file_in, "%d\n", &nrc);
		if (nr_channel != NULL) (*nr_channel) = nrc;
		fscanf(file_in, "%f\n", &imax);
		fread(image, sizeof(float), h*w*nrc, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	default:
		break;
	}
	fclose(file_in);
	return (is_ppm);
}
/*P2(PGM ASCII format), P3(is_ppm ASCII format), P5 (PGM Binary format), P6 (is_ppm Binary format)*/
byte* loadimage(char* file_name, int &h, int &w, bool &is_ppm)
{
	FILE * file_in;
	char line[LEN_MAX];
	int	i, imax;
	byte *image = NULL;
	file_in=fopen(file_name,"rb");
	//fopen(&file_in, file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		getchar();
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		//fscanf(file_in,"%d\n",&i);
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		getchar();
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf(line, "%d %d\n",	&w,	&h);
			sscanf(line, "%d %d\n", &w, &h);
			break;
		}
	}
	//fscanf(file_in, "%d\n", &imax);
	fscanf(file_in, "%d\n", &imax);
	switch (i)
	{
	case 2:
		image = get_ascii_pgm(file_in, w, h);
		is_ppm = QX_DEF_IS_PGM;
		break;
	case 3:
		image = get_ascii_ppm(file_in, w, h);
		is_ppm = QX_DEF_IS_PPM;
		break;
	case 5:
		image = get_binary_pgm(file_in, w, h);
		is_ppm = QX_DEF_IS_PGM;
		break;
	case 6:
		image = get_binary_ppm(file_in, w, h);
		is_ppm = QX_DEF_IS_PPM;
		break;
	default:
		break;
	}
	fclose(file_in);
	return image;
}
int loadimage(byte* out, char* file_name, int h_in, int w_in)
{
	FILE * file_in;
	char line[LEN_MAX];
	int	i, imax, h, w;
	file_in=fopen(file_name,"rb");
	//fopen(&file_in, file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		getchar();
		return(-1);
	}
	if (fgetc(file_in) == 'P')
		//fscanf(file_in,"%d\n",&i);
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		getchar();
		return(-1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf(line, "%d %d\n",	&w,	&h);
			sscanf(line, "%d %d\n", &w, &h);
			break;
		}
	}
	if (h != h_in || w != w_in)
	{
		printf("The allocated memory for the image is not correct!! It should be [%dx%d].\n", h, w);
		getchar();
		return(-1);
	}
	//fscanf(file_in, "%d\n", &imax);
	fscanf(file_in, "%d\n", &imax);
	switch (i)
	{
	case 2:
		get_ascii_pgm(out, file_in, w, h);
		break;
	case 3:
		get_ascii_ppm(out, file_in, w, h);
		break;
	case 5:
		fread(out, sizeof(byte), w*h, file_in);
		break;
	case 6:
		fread(out, sizeof(byte), h*w * 3, file_in);
		break;
	default:
		break;
	}
	fclose(file_in);
	return(0);
}
int loadimage(float* out, byte *out_u, char* file_name, int h_in, int w_in)
{
	FILE * file_in;
	char line[LEN_MAX];
	int	i, imax, h, w;
	file_in=fopen(file_name,"rb");
	//fopen(&file_in, file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		getchar();
		return(-1);
	}
	if (fgetc(file_in) == 'P')
		//fscanf(file_in,"%d\n",&i);
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		getchar();
		return(-1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf(line, "%d %d\n",	&w,	&h);
			sscanf(line, "%d %d\n", &w, &h);
			break;
		}
	}
	if (h != h_in || w != w_in)
	{
		printf("The allocated memory for the image is not correct!! It should be [%dx%d].\n", h, w);
		getchar();
		return(-1);
	}
	//fscanf(file_in, "%d\n", &imax);
	fscanf(file_in, "%d\n", &imax);
	switch (i)
	{
	case 2:
		get_ascii_pgm(out_u, file_in, w, h);
		for (i = 0; i<h_in*w_in; i++) out[i] = (float)out_u[i];
		break;
	case 3:
		get_ascii_ppm(out_u, file_in, w, h);
		for (i = 0; i<h_in*w_in * 3; i++) out[i] = (float)out_u[i];
		break;
	case 5:
		fread(out_u, sizeof(byte), w*h, file_in);
		for (i = 0; i<h_in*w_in; i++) out[i] = (float)out_u[i];
		break;
	case 6:
		fread(out_u, sizeof(byte), h*w * 3, file_in);
		for (i = 0; i<h_in*w_in * 3; i++)
			out[i] = (float)out_u[i];
		break;
	default:
		break;
	}
	fclose(file_in);
	return(0);
}
byte* get_binary_ppm(FILE* in, int h, int w)
{
	byte *out = NULL;
	out = (byte*)malloc(h*w * 3);
	fread(out, sizeof(byte), h*w * 3, in);
	return out;
}
byte* get_binary_pgm(FILE* file_in, int h, int w)
{
	byte *image = NULL;
	image = (byte*)malloc(w*h);
	fread(image, sizeof(byte), w*h, file_in);
	return image;
}
byte* get_ascii_ppm(FILE* file_in, int h, int w)
{
	byte *image = NULL; int	y, x; int	color_r, color_g, color_b;
	image = (byte*)malloc(w*h * 3);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
	{
		//if(fscanf(file_in,"%d%d%d",&color_r,&color_g,&color_b)!=3)
		if (fscanf(file_in, "%d%d%d", &color_r, &color_g, &color_b) != 3)
		{
			printf("error in reading file.\n");
			getchar();
			exit(0);
		}
		image[(y*w + x) * 3] = (byte)color_r;
		image[(y*w + x) * 3 + 1] = (byte)color_g;
		image[(y*w + x) * 3 + 2] = (byte)color_b;
	}
	return image;
}
byte* get_ascii_pgm(FILE* file_in, int h, int w)
{
	byte *image = NULL;
	int	y, x, lum;
	image = (byte*)malloc(w*h);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
	{
		//if(fscanf(file_in,"%d",&lum)!=1)
		if (fscanf(file_in, "%d", &lum) != 1)
		{
			printf("error in reading file.\n");
			getchar();
			exit(0);
		}
		image[y*w + x] = (byte)lum;
	}
	return image;
}
void get_ascii_ppm(byte *image, FILE* file_in, int h, int w)
{
	int	y, x; int color_r, color_g, color_b;
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
	{
		//if(fscanf(file_in,"%d%d%d",&color_r,&color_g,&color_b)!=3)
		if (fscanf(file_in, "%d%d%d", &color_r, &color_g, &color_b) != 3)
		{
			printf("error in reading file.\n");
			getchar();
			exit(0);
		}
		image[(y*w + x) * 3] = (byte)color_r;
		image[(y*w + x) * 3 + 1] = (byte)color_g;
		image[(y*w + x) * 3 + 2] = (byte)color_b;
	}
}
void get_ascii_pgm(byte *image, FILE* file_in, int h, int w)
{
	int	y, x, lum;
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
	{
		//if(fscanf(file_in,"%d",&lum)!=1)
		if (fscanf(file_in, "%d", &lum) != 1)
		{
			printf("error in reading file.\n");
			getchar();
			exit(0);
		}
		image[y*w + x] = (byte)lum;
	}
}
void saveimage_ppm(char *file_name, byte *image, int h, int w, bool is_binary)
{
	if (is_binary == true) write_binary_ppm(file_name, image, h, w);
	else write_ascii_ppm(file_name, image, h, w);
}
void saveimage_pgm(char *file_name, byte *image, int h, int w, bool is_binary)
{
	if (is_binary == true) write_binary_pgm(file_name, image, h, w);
	else write_ascii_pgm(file_name, image, h, w);
}
void write_binary_ppm(char* file_name, byte *in, int h, int w)
{
	FILE * file_out;
	file_out=fopen(file_name,"wb");
	//fopen(&file_out, file_name, "wb");
	fprintf(file_out, "P6\n%d %d\n%d\n", w, h, 255);
	fwrite(in, sizeof(byte), h*w * 3, file_out);
	fclose(file_out);
}
void write_binary_pgm(char* file_name, byte *image, int h, int w)
{
	FILE* file_out;
	file_out=fopen(file_name,"wb");
	//fopen(&file_out, file_name, "wb");
	if (file_out == NULL)
	{
		printf("Can't open the file %s, exit!\n", file_name);
		getchar();
		exit(0);
	}
	fprintf(file_out, "P5\n%d %d\n%d\n", w, h, 255);
	fwrite(image, sizeof(byte), w*h, file_out);
	fclose(file_out);
}
void write_ascii_ppm(char* file_name, byte *image, int h, int w)
{
	FILE* file_out; int y, x;
	file_out=fopen(file_name,"w");
	//fopen(&file_out, file_name, "w");
	fprintf(file_out, "P3\n%d %d\n%d\n", w, h, 255);
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			fprintf(file_out, "%d", image[(y*w + x) * 3 + 0]);
			fprintf(file_out, "%c", ' ');
			fprintf(file_out, "%d", image[(y*w + x) * 3 + 1]);
			fprintf(file_out, "%c", ' ');
			fprintf(file_out, "%d", image[(y*w + x) * 3 + 2]);
			fprintf(file_out, "%c", ' ');
		}
	}
	fclose(file_out);
}
void write_ascii_pgm(char* file_name, byte *image, int h, int w)
{
	FILE* file_out; int y, x;
	file_out=fopen(file_name,"w");
	//fopen(&file_out, file_name, "w");
	fprintf(file_out, "P2\n%d %d\n%d\n", w, h, 255);
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			fprintf(file_out, "%d", image[y*w + x]);
			fprintf(file_out, "%c", ' ');
		}
	}
	fclose(file_out);
}
/*extended function*/
void qx_image_size(char* file_name, int &h, int &w, int *nr_channel)
{
	FILE * file_in; int nrc;
	char line[LEN_MAX];
	int	i;
	byte *image = NULL;
	file_in=fopen(file_name,"rb");
	//fopen(&file_in, file_name, "rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		getchar();
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		//fscanf(file_in,"%d\n",&i);
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		getchar();
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf(line, "%d %d\n",	&w,	&h);
			sscanf(line, "%d %d\n", &w, &h);
			if (i == 9 && nr_channel != NULL)
			{
				fscanf(file_in, "%d\n", &nrc);
				(*nr_channel) = nrc;
			}
			break;
		}
	}
}
byte* verticalflip(int	width, int imageHeight, byte* pixel)
{
	byte *textureFlip = NULL;
	int	y, x;
	textureFlip = (byte*)malloc(width*imageHeight * 3);
	for (y = 0; y<imageHeight; y++)
	for (x = 0; x<width; x++)
	{
		textureFlip[(y*width + x) * 3] = pixel[((imageHeight - y - 1)*width + x) * 3];
		textureFlip[(y*width + x) * 3 + 1] = pixel[((imageHeight - y - 1)*width + x) * 3 + 1];
		textureFlip[(y*width + x) * 3 + 2] = pixel[((imageHeight - y - 1)*width + x) * 3 + 2];
	}
	free(pixel);
	return textureFlip;
}
float ** loadimage_pgm(char *file_name, int &h, int &w)
{
	float ** image_2d; byte* image; bool is_ppm; int y, x, counter;
	is_ppm = QX_DEF_IS_PPM; counter = 0;
	image = loadimage(file_name, h, w, is_ppm);
	image_2d = qx_allocf(h, w);
	if (is_ppm == QX_DEF_IS_PPM)
	{
		for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		{
			image_2d[y][x] = (0.299F)*float(int(image[counter])) + (0.587F)*float(int(image[counter + 1])) + (0.114F)*float(int(image[counter + 2]));
			counter = counter + 3;
		}
	}
	else
	{
		for (y = 0; y<h; y++) for (x = 0; x<w; x++)
			image_2d[y][x] = float(int(image[counter++]));
	}
	free(image);
	return image_2d;
}
float *** loadimage_ppm(char *file_name, int &h, int &w)
{
	float *** image_3d; byte* image; bool is_ppm; int y, x, k, counter;
	is_ppm = QX_DEF_IS_PPM; counter = 0;
	image = loadimage(file_name, h, w, is_ppm);
	image_3d = qx_allocf_3(h, w, 3);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
	for (k = 0; k<3; k++) image_3d[y][x][k] = float(int(image[counter++]));
	free(image);
	return image_3d;
}
byte *** loadimage_ppm_u(char *file_name, int &h, int &w)
{
	byte *** image_3d; byte* image; bool is_ppm; int y, x, k, counter;
	is_ppm = QX_DEF_IS_PPM; counter = 0;
	image = loadimage(file_name, h, w, is_ppm);
	image_3d = qx_allocu_3(h, w, 3);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++) for (k = 0; k<3; k++)
		image_3d[y][x][k] = image[counter++];
	free(image);
	return image_3d;
}
byte **loadimage_pgm_u(char *file_name, int &h, int &w)
{
	byte ** image_2d; byte* image; bool is_ppm; int y, x, counter; double lum;
	is_ppm = QX_DEF_IS_PPM; counter = 0;
	image = loadimage(file_name, h, w, is_ppm);
	image_2d = qx_allocu(h, w);
	if (is_ppm == QX_DEF_IS_PPM)
	{
		for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		{
			lum = 0.299*double(image[counter]) + 0.587*double(image[counter + 1]) + 0.114*double(image[counter + 2]);
			image_2d[y][x] = byte(lum + 0.5);
			counter = counter + 3;
		}
	}
	else
	{
		for (y = 0; y<h; y++) for (x = 0; x<w; x++)
			image_2d[y][x] = image[counter++];
	}
	free(image);
	return image_2d;
}

int **loadimage_pgm_i(char *file_name, int &h, int &w)
{
	int **image; byte **image_2d; int y, x;
	image_2d = loadimage_pgm_u(file_name, h, w);
	image = qx_alloci(h, w);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image[y][x] = (int)image_2d[y][x];
	qx_freeu(image_2d);
	return image;
}
void saveimage_pgm(char *file_name, float **image, int h, int w, int scale)
{
	byte *image_u; int counter, y, x;
	counter = 0;
	image_u = new byte[h*w];
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image_u[counter++] = byte(image[y][x] * scale + 0.5f);
	write_binary_pgm(file_name, image_u, h, w);
	delete image_u;
}
void saveimage_pgm_ascii(char *file_name, float **image, int h, int w, int scale)
{
	byte *image_u; int counter, y, x;
	counter = 0;
	image_u = new byte[h*w];
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image_u[counter++] = byte(image[y][x] * scale + 0.5f);
	write_ascii_pgm(file_name, image_u, h, w);
	delete image_u;
}
void saveimage_pgm_ascii(char *file_name, int **image, int h, int w, int scale)
{
	byte *image_u; int counter, y, x;
	counter = 0;
	image_u = new byte[h*w];
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image_u[counter++] = byte(image[y][x] * scale);
	write_ascii_pgm(file_name, image_u, h, w);
	delete image_u;
}
void saveimage_pgm(char *file_name, byte **image, int h, int w, int scale)
{
	byte *image_u; int counter, y, x;
	counter = 0;
	image_u = new byte[h*w];
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image_u[counter++] = image[y][x] * scale;
	write_binary_pgm(file_name, image_u, h, w);
	delete image_u;
}
void saveimage_pgm(char *file_name, int **image, int h, int w, int scale)
{
	byte *image_u; int counter, y, x;
	counter = 0;
	image_u = new byte[h*w];
	for (y = 0; y<h; y++) for (x = 0; x<w; x++)
		image_u[counter++] = byte(image[y][x] * scale + 0.5f);
	write_binary_pgm(file_name, image_u, h, w);
	delete image_u;
}
void saveimage_ppm(char *file_name, float ***image, int h, int w, int scale)
{
	int len = h*w * 3; float *image_ = image[0][0];
	byte *in = new byte[len];
	for (int i = 0; i<len; i++)
		in[i] = (byte)(image_[i] * scale);
	write_binary_ppm(file_name, in, h, w);
	delete[] in;
}
void saveimage_ppm(char *file_name, double ***image, int h, int w, int scale)
{
	int len = h*w * 3; double *image_ = image[0][0];
	byte *in = new byte[len];
	for (int i = 0; i<len; i++)
		in[i] = (byte)(image_[i] * scale);
	write_binary_ppm(file_name, in, h, w);
	delete[] in;
}
void saveimage_ppm(char *file_name, byte ***image, int h, int w, int scale)
{
	write_binary_ppm(file_name, image[0][0], h, w);
}
float ***loadimage_ftif(char *file_name, int &h, int &w, int &nr_channel)
{
	FILE* file_in;
	if ((file_in = fopen( file_name, "r") )!= 0)
	{
		printf("Please check input file_name: %s\n", file_name);
		exit(0);
	}
	float ***image = NULL; int y, x, c; float color;
	if (fscanf(file_in, "%d%d%d", &h, &w, &nr_channel) != 3)
	{
		printf("error in reading file.\n");
		getchar();
		exit(0);
	}
	image = qx_allocf_3(h, w, nr_channel);
	for (y = 0; y<h; y++) for (x = 0; x<w; x++) for (c = 0; c<nr_channel; c++)
	{
		fscanf(file_in, "%f", &color);
		image[0][0][(y*w + x)*nr_channel + c] = color;
	}
	fclose(file_in);
	return(image);
}

void qx_saveimage(char* filename, byte *image, int h, int w, int channel)
{
	FILE* file_out; byte maxx = 255;
	file_out=fopen( filename, "wb");
	if (channel == 1) fprintf(file_out, "P5\n%d %d\n%d\n", w, h, maxx);
	//else if(channel==3) fprintf(file_out,"P6\n%d %d\n%d\n",w,h,maxx);
	else if (channel == 3) fprintf(file_out, "P6\n%d %d\n%d\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%d\n", w, h, channel, maxx);
	fwrite(image, sizeof(byte), h*w*channel, file_out);

	fclose(file_out);
}
void qx_saveimage(char* filename, float *image, int h, int w, int channel)
{
	FILE* file_out; int i; float maxx;
	file_out=fopen(filename,"wb");
	//fopen(&file_out, filename, "wb");
	maxx = image[0];
	for (i = 0; i<h*w*channel; i++) maxx = max(maxx, image[i]);
	if (channel == 1) fprintf(file_out, "P7\n%d %d\n%f\n", w, h, maxx);
	else if (channel == 3) fprintf(file_out, "P8\n%d %d\n%f\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%f\n", w, h, channel, maxx);
	fwrite(image, sizeof(float), w*h*channel, file_out);
	fclose(file_out);
}
void qx_saveimage(char*filename, short*image, int h, int w, int channel)
{
	FILE* file_out; int i; short maxx;
	file_out=fopen(filename,"wb");
	//fopen(&file_out, filename, "wb");
	maxx = image[0];
	for (i = 0; i<h*w*channel; i++) maxx = max(maxx, image[i]);
	if (channel == 1) fprintf(file_out, "P7\n%d %d\n%f\n", w, h, maxx);
	else if (channel == 3) fprintf(file_out, "P8\n%d %d\n%f\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%f\n", w, h, channel, maxx);
	fwrite(image, sizeof(short), w*h*channel, file_out);
	fclose(file_out);
}
void qx_saveimage(char* filename, double *image, int h, int w, int channel)
{
	FILE* file_out; int i; float maxx;
	float *image_f = new float[h*w*channel];
	for (i = 0; i<h*w*channel; i++) image_f[i] = (float)image[i];
	file_out=fopen(filename,"wb");
	//fopen(&file_out, filename, "wb");
	maxx = image_f[0];
	for (i = 0; i<h*w*channel; i++) maxx = max(maxx, image_f[i]);
	if (channel == 1) fprintf(file_out, "P7\n%d %d\n%f\n", w, h, maxx);
	else if (channel == 3) fprintf(file_out, "P8\n%d %d\n%f\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%f\n", w, h, channel, maxx);
	fwrite(image_f, sizeof(float), w*h*channel, file_out);
	fclose(file_out);
	delete[] image_f; image_f = NULL;
}
int qx_loadimage(char* filename, byte *image, int h, int w, int *nr_channel)
{
	FILE * file_in; int nrc;
	char line[LEN_MAX];
	int	i; int imax, hc, wc;
	byte *image_ = image;
	file_in=fopen( filename, "rb");
	if (!file_in)
	{
		printf("Please check input filename: %s\n", filename);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			sscanf(line, "%d %d\n", &wc, &hc);
			break;
		}
	}
	char str_tmp[100];
	switch (i)
	{
	case 5:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 1;
		memset(image, 0, sizeof(byte)*h*w);
		fread(image, sizeof(byte), h*w, file_in);
		break;
	case 6:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 3;
		memset(image, 0, sizeof(byte)*h*w * 3);
		fread(image, sizeof(byte), h*w * 3, file_in);
		break;
	case 2:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++)
		{
			//if(fscanf(file_in,"%d",&imax)!=1){printf("error in reading file.\n");getchar();exit(0);}
			fscanf(file_in, "%d", &imax);
			*image_++ = imax;
		}
		break;
	case 3:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		int cr, cg, cb;
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++)
		{
			//if(fscanf(file_in,"%d%d%d",&cr,&cg,&cb)!=3){printf("error in reading file.\n");getchar();exit(0);}
			fscanf(file_in, "%d%d%d", &cr, &cg, &cb);
			*image_++ = cr; *image_++ = cg; *image_++ = cb;
		}
		break;
	case 9:
		fgets(str_tmp, 100, file_in);
		nrc = atoi(str_tmp);
		fscanf(file_in, "%d\n", &nrc);
		if (nr_channel != NULL) (*nr_channel) = nrc;
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		fread(image, sizeof(byte), h*w*nrc, file_in);
		break;
	default:
		printf("Can not open image [%s]!!\n", filename);
		break;
	}
	fclose(file_in);
	return (0);
}
int qx_loadimage(char* filename, float*image, int h, int w, int *nr_channel)
{
	FILE * file_in; int is_ppm; int nrc;
	char line[LEN_MAX];
	int	i; double imax; int hi, wi;
	is_ppm = QX_DEF_IS_PGM;
	file_in=fopen(filename, "rb");
	if (!file_in)
	{
		printf("Please check input filename: %s\n", filename);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			sscanf(line, "%d %d\n", &wi, &hi);
			break;
		}
	}
	char str_tmp[100];
	switch (i)
	{
	case 7:
		fgets(str_tmp, 100, file_in);
		imax = atof(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 1;
		memset(image, 0, sizeof(float)*h*w);
		fread(image, sizeof(float), h*w, file_in);
		is_ppm = QX_DEF_IS_PGM;
		break;
	case 8:
		fgets(str_tmp, 100, file_in);
		imax = atof(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 3;
		memset(image, 0, sizeof(float)*h*w * 3);
		fread(image, sizeof(float), h*w * 3, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	case 9:
		fgets(str_tmp, 100, file_in);
		nrc = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = nrc;
		fgets(str_tmp, 100, file_in);
		imax = atof(str_tmp);
		fread(image, sizeof(float), h*w*nrc, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	default:
		break;
	}
	fclose(file_in);
	return (is_ppm);
}
int qx_loadimage(char* filename, short*image, int h, int w, int *nr_channel)
{
	FILE * file_in; int is_ppm; int nrc;
	char line[LEN_MAX];
	int	i; int imax; int hi, wi;
	is_ppm = QX_DEF_IS_PGM;
	file_in=fopen( filename, "rb");
	if (!file_in)
	{
		printf("Please check input filename: %s\n", filename);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		fscanf(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			sscanf(line, "%d %d\n", &wi, &hi);
			break;
		}
	}
	char str_tmp[100];
	switch (i)
	{
	case 7:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 1;
		memset(image, 0, sizeof(short)*h*w);
		fread(image, sizeof(short), h*w, file_in);
		is_ppm = QX_DEF_IS_PGM;
		break;
	case 8:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 3;
		memset(image, 0, sizeof(short)*h*w * 3);
		fread(image, sizeof(short), h*w * 3, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	case 9:
		fgets(str_tmp, 100, file_in);
		nrc = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = nrc;
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		fread(image, sizeof(short), h*w*nrc, file_in);
		is_ppm = QX_DEF_IS_PPM;
		break;
	default:
		break;
	}
	fclose(file_in);
	return (is_ppm);
}


qx_tree_filter::qx_tree_filter()
{
}
qx_tree_filter::~qx_tree_filter()
{
	clean();
}
void qx_tree_filter::clean()
{
}
int qx_tree_filter::init(int h, int w, int nr_channel, double sigma_range, int nr_neighbor)
{
	m_h = h; m_w = w; m_nr_pixel = m_h*m_w;
	m_mst.init(m_h, m_w, nr_channel, nr_neighbor);
	update_table(sigma_range);
	return(0);
}
void qx_tree_filter::update_table(double sigma_range)
{
	sigma_range = max(0.01, sigma_range);
	for (int i = 0; i <= QX_DEF_CHAR_MAX; i++) m_table[i] = exp(-double(i) / (QX_DEF_CHAR_MAX*sigma_range));//weight table
}
int qx_tree_filter::build_tree(byte*texture)
{
	m_mst.mst(texture);
	m_mst_parent = m_mst.get_parent();
	m_mst_nr_child = m_mst.get_nr_child();
	m_mst_children = m_mst.get_children();
	m_mst_rank = m_mst.get_rank();
	m_mst_weight = m_mst.get_weight();
	m_node_id = m_mst.get_node_id();

	return(0);
}
template<typename T>void qx_tree_filter::combine_tree(T*image_filtered)
{
	//double*value = m_mst_value_sum_aggregated_from_parent_to_child;
	//for (int i = 0; i<m_nr_pixel; i++)
	//{
	//	*image_filtered++ = *value++;//every slices will have the same weight, thus don't need normalization.
	//}
}
template<typename T>void qx_tree_filter::init_tree_value(T*image, bool compute_weight)
{
	//memset(m_mst_value_sum_aggregated_from_child_to_parent, 0, sizeof(double)*m_nr_pixel);
	//memset(m_mst_value_sum_aggregated_from_parent_to_child, 0, sizeof(double)*m_nr_pixel);
	//if (compute_weight)
	//{
	//	memset(m_mst_weight_sum_aggregated_from_child_to_parent, 0, sizeof(double)*m_nr_pixel);
	//	memset(m_mst_weight_sum_aggregated_from_parent_to_child, 0, sizeof(double)*m_nr_pixel);/*nodes*/
	//}
	//double*value_to_be_filtered = m_mst_value_to_be_filtered;
	//for (int i = 0; i<m_nr_pixel; i++)
	//{
	//	*value_to_be_filtered++ = double(*image++);
	//}
}
int qx_tree_filter::filter(double*cost, double*cost_backup, int nr_plane)
{
	memcpy(cost_backup, cost, sizeof(double)*m_h*m_w*nr_plane);
	int*node_id = m_node_id;
	int*node_idt = &(node_id[m_nr_pixel - 1]);
	for (int i = 0; i<m_nr_pixel; i++)
	{
		int id = *node_idt--;
		int id_ = id*nr_plane;
		int nr_child = m_mst_nr_child[id];
		if (nr_child>0)
		{
			double*value_sum = &(cost_backup[id_]);
			for (int j = 0; j<nr_child; j++)
			{
				int id_child = m_mst_children[id][j];
				int id_child_ = id_child*nr_plane;
				double weight = m_table[m_mst_weight[id_child]];
				//value_sum+=m_mst_value_sum_aggregated_from_child_to_parent[id_child]*weight;
				double*value_child = &(cost_backup[id_child_]);
				for (int k = 0; k<nr_plane; k++)
				{
					value_sum[k] += (*value_child++)*weight;
				}
			}
		}
		//else
		//{
		//	memcpy(&(cost_backup[id_]),&(cost[id_]),sizeof(double)*nr_plane);
		//}
		//printf("[id-value-weight]: [%d - %3.3f - %3.3f]\n",id,m_mst_[id].value_sum_aggregated_from_child_to_parent,m_mst_[id].weight_sum_aggregated_from_child_to_parent);
	}
	int*node_id0 = node_id;
	int tree_parent = *node_id0++;
	int tree_parent_ = tree_parent*nr_plane;
	memcpy(&(cost[tree_parent_]), &(cost_backup[tree_parent_]), sizeof(double)*nr_plane);
	for (int i = 1; i<m_nr_pixel; i++)//K_00=f(0,00)[K_0-f(0,00)J_00]+J_00, K_00: new value, J_00: old value, K_0: new value of K_00's parent
	{
		int id = *node_id0++;
		int id_ = id*nr_plane;
		int parent = m_mst_parent[id];
		int parent_ = parent*nr_plane;

		double*value_parent = &(cost[parent_]);//K_0
		double*value_current = &(cost_backup[id_]);//J_00
		double*value_final = &(cost[id_]);
		double weight = m_table[m_mst_weight[id]];//f(0,00)

		for (int k = 0; k<nr_plane; k++)
		{
			double vc = *value_current++;
			*value_final++ = weight*((*value_parent++) - weight*vc) + vc;
		}
		//printf("Final: [id-value-weight]: [%d - %3.3f - %3.3f]\n",id,m_mst_[id].value_sum_aggregated_from_parent_to_child,m_mst_[id].weight_sum_aggregated_from_parent_to_child);
	}
	return(0);
}


