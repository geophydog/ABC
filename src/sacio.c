/*******************************************************************************
 *                                  sacio.c                                    *
 *  SAC I/O functions:                                                         *
 *      read_sac_head    read SAC header                                       *
 *      read_sac         read SAC binary data                                  *
 *      read_sac_xy      read SAC binary XY data                               *
 *      read_sac_pdw     read SAC data in a partial data window (cut option)   *
 *      write_sac        Write SAC binary data                                 *
 *      write_sac_xy     Write SAC binary XY data                              *
 *      new_sac_head     Create a new minimal SAC header                       *
 *      sac_head_index   Find the offset of specified SAC head fields          *
 *      issac            Check if a file in in SAC format                      *
 *                                                                             *
 *  Author: Dongdong Tian @ USTC                                               *
 *                                                                             *
 *  Revisions:                                                                 *
 *      2014-03-19  Dongdong Tian   Modified from Prof. Lupei Zhu's code       *
 *      2014-08-02  Dongdong Tian   Better function naming and coding style    *
 *      2014-08-13  Dongdong Tian   Add new funtions:                          *
 *                                  - read_sac_xy                              *
 *                                  - write_sac_xy                             *
 *                                  - sac_head_index                           *
 *      2016-03-01  Dongdong Tian   Add new function: issac                    *
 *                                                                             *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fftw3.h>
//#include </home/feng_xuping/MY_LIB/FFTW3/include/fftw3.h>
#include "sacio.h"

/* function prototype for local use */
static void    byte_swap       (char *pt, size_t n);
static int     check_sac_nvhdr (const int nvhdr);
static void    map_chdr_in     (char *memar, char *buff);
static int     read_head_in    (const char *name, SACHEAD *hd, FILE *strm);
static void    map_chdr_out    (char *memar, char *buff);
static int     write_head_out  (const char *name, SACHEAD hd, FILE *strm);

/* a SAC structure containing all null values */
static SACHEAD sac_null = {
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }
};

/*
 *  read_sac_head
 *
 *  Description: Read binary SAC header from file.
 *
 *  IN:
 *      const char *name : File name
 *  OUT:
 *      SACHEAD    *hd   : SAC header
 *
 *  Return: 0 if success, -1 if failed
 *
 */
int read_sac_head(const char *name, SACHEAD *hd)
{
    FILE    *strm;
    int     lswap;

    if ((strm = fopen(name, "rb")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", name);
        return -1;
    }

    lswap = read_head_in(name, hd, strm);

    fclose(strm);

    return ((lswap == -1) ? -1 : 0);
}

/*
 *  read_sac
 *
 *  Description: Read binary SAC data from file.
 *
 *  IN:
 *      const char *name : file name
 *  OUT:
 *      SACHEAD    *hd   : SAC header to be filled
 *  Return: float pointer to the data array, NULL if failed.
 *
 */
float *read_sac(const char *name, SACHEAD *hd)
{
    FILE    *strm;
    float   *ar;
    int     lswap;
    size_t  sz;

    if ((strm = fopen(name, "rb")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", name);
        return NULL;
    }

    lswap = read_head_in(name, hd, strm);

    if (lswap == -1) {
        fclose(strm);
        return NULL;
    }

    sz = (size_t) hd->npts * SAC_DATA_SIZEOF;
    if (hd->iftype == IXY) sz *= 2;

    if ((ar = (float *)malloc(sz)) == NULL) {
        fprintf(stderr, "Error in allocating memory for reading %s\n", name);
        fclose(strm);
        return NULL;
    }

    if (fread((char*)ar, sz, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC data %s\n", name);
        free(ar);
        fclose(strm);
        return NULL;
    }
    fclose(strm);

    if (lswap == TRUE) byte_swap((char*)ar, sz);

    return ar;
}

/*
 *  read_sac_xy
 *
 *  Description:    read SAC XY binary file
 *
 *  IN:
 *      const char *name    :   file name
 *  OUT:
 *      SACHEAD *hd         :   SAC head to be filled
 *      float *xdata        :   pointer for X
 *      float *ydata        :   pointer for Y
 *
 *  Return: 0 for success, -1 for fail
 *
 */
int read_sac_xy(const char *name, SACHEAD *hd, float *xdata, float *ydata)
{
    float *data;
    size_t npts;

    if ((data = read_sac(name, hd)) == NULL)  return -1;

    npts = (size_t)hd->npts;
    if ((xdata = (float *)malloc(npts*SAC_DATA_SIZEOF)) == NULL) {
        fprintf(stderr, "Error in allocating memory for %s\n", name);
        return -1;
    }
    if ((ydata = (float *)malloc(npts*SAC_DATA_SIZEOF)) == NULL) {
        fprintf(stderr, "Error in allocating memory for %s\n", name);
        return -1;
    }

    memcpy(xdata, data     , npts*SAC_DATA_SIZEOF);
    memcpy(ydata, data+npts, npts*SAC_DATA_SIZEOF);

    free(data);
    return 0;
}

/*
 *  write_sac
 *
 *  Description:    write binary SAC data
 *
 *  IN:
 *      const char *name    :   file name
 *      SACHEAD     hd      :   header
 *      const float *ar     :   float data array
 *
 *  Return:
 *      -1  :   fail
 *      0   :   succeed
 *
 */
int write_sac(const char *name, SACHEAD hd, const float *ar)
{
    FILE    *strm;
    size_t  sz;

    if ((strm = fopen(name, "wb")) == NULL) {
        fprintf(stderr, "Error in opening file for writing %s\n", name);
        return -1;
    }

    if (write_head_out(name, hd, strm) == -1) {
        fclose(strm);
        return -1;
    }

    sz = (size_t)hd.npts * SAC_DATA_SIZEOF;
    if (hd.iftype == IXY) sz *= 2;

    if (fwrite(ar, sz, 1, strm) != 1) {
        fprintf(stderr, "Error in writing SAC data for writing %s\n", name);
        fclose(strm);
        return -1;
    }
    fclose(strm);
    return 0;
}

/*
 *  write_sac_xy
 *
 *  Description:    write binary SAC XY data
 *
 *  IN:
 *      const char *name    :   file name
 *      SACHEAD     hd      :   header
 *      const float *xdata  :   float data array for X
 *      const float *ydata  :   float data array for Y
 *
 *  Return:
 *      -1  :   fail
 *      0   :   succeed
 *
 */
int write_sac_xy(const char *name, SACHEAD hd, const float *xdata, const float *ydata)
{
    float *ar;
    int npts;
    int error;
    size_t sz;

    npts = hd.npts;
    sz = (size_t)npts * SAC_DATA_SIZEOF;

    if ((ar = (float *)malloc(sz*2)) == NULL) {
        fprintf(stderr, "Error in allocating memory for file %s\n", name);
        return -1;
    }
    memcpy(ar,      xdata, sz);
    memcpy(ar+npts, ydata, sz);

    /* needed for XY data */
    hd.iftype = IXY;
    hd.leven = FALSE;

    error = write_sac(name, hd, ar);

    free(ar);
    return error;
}

/*
 *  read_sac_pdw
 *
 *  Description:
 *      Read portion of data from file.
 *
 *  Arguments:
 *      const char  *name   :   file name
 *      SACHEAD     *hd     :   SAC header to be filled
 *      int         tmark   :   time mark in SAC header
 *                                  -5  ->  b;
 *                                  -4  ->  e;
 *                                  -3  ->  o;
 *                                  -2  ->  a;
 *                                  0-9 ->  Tn;
 *                                  others -> t=0;
 *      float       t1      :   begin time is tmark + t1
 *      float       t2      :   end time is tmark + t2
 *
 *  Return:
 *      float pointer to the data array, NULL if failed.
 *
 */
float *read_sac_pdw(const char *name, SACHEAD *hd, int tmark, float t1, float t2)
{
    FILE    *strm;
    int     lswap;
    float   tref;
    int     nt1, nt2, npts, nn;
    float   *ar, *fpt;

    if ((strm = fopen(name, "rb")) == NULL) {
        fprintf(stderr, "Error in opening %s\n", name);
        return NULL;
    }

    lswap = read_head_in(name, hd, strm);

    if (lswap == -1) {
        fclose(strm);
        return NULL;
    }

    nn = (int)((t2-t1)/hd->delta);
    if (nn<=0 || (ar = (float *)calloc((size_t)nn, SAC_DATA_SIZEOF)) == NULL) {
        fprintf(stderr, "Errorin allocating memory for reading %s n=%d\n", name, nn);
        fclose(strm);
        return NULL;
    }

    tref = 0.;
    if ((tmark>=-5&&tmark<=-2) || (tmark>=0 && tmark<=9)) {
        tref = *((float *) hd + TMARK + tmark);
        if (fabs(tref+12345.)<0.1) {
            fprintf(stderr, "Time mark undefined in %s\n", name);
            free(ar);
            fclose(strm);
            return NULL;
        }
    }
    t1 += tref;
    nt1 = (int)((t1 - hd->b) / hd->delta);
    nt2 = nt1 + nn;
    npts = hd->npts;
    hd->npts = nn;
    hd->b   = t1;
    hd->e   = t1 + nn * hd->delta;

    if (nt1>npts || nt2 <0) return ar;    /* return zero filled array */
    /* maybe warnings are needed! */

    if (nt1<0) {
        fpt = ar - nt1;
        nt1 = 0;
    } else {
        if (fseek(strm, nt1*SAC_DATA_SIZEOF, SEEK_CUR) < 0) {
            fprintf(stderr, "Error in seek %s\n", name);
            free(ar);
            fclose(strm);
            return NULL;
        }
        fpt = ar;
    }
    if (nt2>npts) nt2 = npts;
    nn = nt2 - nt1;

    if (fread((char *)fpt, (size_t)nn * SAC_DATA_SIZEOF, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC data %s\n", name);
        free(ar);
        fclose(strm);
        return NULL;
    }
    fclose(strm);

    if (lswap == TRUE) byte_swap((char*)ar, (size_t)nn*SAC_DATA_SIZEOF);

    return ar;
}

/*
 *  new_sac_head
 *
 *  Description: create a new SAC header with required fields
 *
 *  IN:
 *      float   dt  :   sample interval
 *      int     ns  :   number of points
 *      float   b0  :   starting time
 */
SACHEAD new_sac_head(float dt, int ns, float b0)
{
    SACHEAD hd = sac_null;
    hd.delta    =   dt;
    hd.npts     =   ns;
    hd.b        =   b0;
    hd.o        =   0.;
    hd.e        =   b0+(ns-1)*dt;
    hd.iztype   =   IO;
    hd.iftype   =   ITIME;
    hd.leven    =   TRUE;
    hd.nvhdr    =   SAC_HEADER_MAJOR_VERSION;
    return hd;
}

/*
 *  sac_head_index
 *
 *  Description: return the index of a specified sac head field
 *
 *  In:
 *      const char *name    :   name of sac head field
 *  Return:
 *      index of a specified field in sac head
 *
 */
int sac_head_index(const char *name)
{
    const char fields[SAC_HEADER_NUMBERS+SAC_HEADER_STRINGS][10] = {
        "delta",    "depmin",   "depmax",   "scale",    "odelta",
        "b",        "e",        "o",        "a",        "internal1",
        "t0",       "t1",       "t2",       "t3",       "t4",
        "t5",       "t6",       "t7",       "t8",       "t9",
        "f",        "resp0",    "resp1",    "resp2",    "resp3",
        "resp4",    "resp5",    "resp6",    "resp7",    "resp8",
        "resp9",    "stla",     "stlo",     "stel",     "stdp",
        "evla",     "evlo",     "evel",     "evdp",     "mag",
        "user0",    "user1",    "user2",    "user3",    "user4",
        "user5",    "user6",    "user7",    "user8",    "user9",
        "dist",     "az",       "baz",      "gcarc",    "internal2",
        "internal3","depmen",   "cmpaz",    "cmpinc",   "xminimum",
        "xmaximum", "yminimum", "ymaximum", "unused1",  "unused2",
        "unused3",  "unused4",  "unused5",  "unused6",  "unused7",
        "nzyear",   "nzjday",   "nzhour",   "nzmin",    "nzsec",
        "nzmsec",   "nvhdr",    "norid",    "nevid",    "npts",
        "internal4","nwfid",    "nxsize",   "nysize",   "unused8",
        "iftype",   "idep",     "iztype",   "unused9",  "iinst",
        "istreg",   "ievreg",   "ievtyp",   "iqual",    "isynth",
        "imagtyp",  "imagsrc",  "unused10", "unused11", "unused12",
        "unused13", "unused14", "unused15", "unused16", "unused17",
        "leven",    "lpspol",   "lovrok",   "lcalda",   "unused18",
        "kstnm",    "kevnm",    "kevnmmore",
        "khole",    "ko",       "ka",
        "kt0",      "kt1",      "kt2",
        "kt3",      "kt4",      "kt5",
        "kt6",      "kt7",      "kt8",
        "kt9",      "kf",       "kuser0",
        "kuser1",   "kuser2",   "kcmpnm",
        "knetwk",   "kdatrd",   "kinst",
    };
    int i;

    for (i=0; i<SAC_HEADER_NUMBERS+SAC_HEADER_STRINGS; i++)
        if ((strcasecmp(name, fields[i]) == 0))  return i;

    return -1;
}

/*
 *  issac
 *
 *  Description: check if a file is in SAC format
 *
 *  In:
 *      const char *name    :   sac filename
 *  Return:
 *      -1 : fail
 *      TRUE  : is a SAC file
 *      FALSE  : not a SAC file
 *
 */
int issac(const char *name)
{
    FILE *strm;
    int nvhdr;

    if ((strm = fopen(name, "rb")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", name);
        return -1;
    }

    if (fseek(strm, SAC_VERSION_LOCATION * SAC_DATA_SIZEOF, SEEK_SET)) return FALSE;
    if (fread(&nvhdr, sizeof(int), 1, strm) != 1) return FALSE;
    if (check_sac_nvhdr(nvhdr) == -1) return FALSE;
    else return TRUE;
}

/******************************************************************************
 *                                                                            *
 *              Functions below are only for local use!                       *
 *                                                                            *
 ******************************************************************************/

/*
 *  byte_swap : reverse the byte order of 4 bytes int/float.
 *
 *  IN:
 *      char    *pt : pointer to byte array
 *      size_t   n  : number of bytes
 *  Return: none
 *
 *  Notes:
 *      For 4 bytes,
 *      byte swapping means taking [0][1][2][3],
 *      and turning it into [3][2][1][0]
 */
static void byte_swap(char *pt, size_t n)
{
    size_t  i   ;
    char    tmp ;
    for (i=0; i<n; i+=4) {
        tmp     =   pt[i+3];
        pt[i+3] =   pt[i];
        pt[i]   =   tmp;

        tmp     =   pt[i+2];
        pt[i+2] =   pt[i+1];
        pt[i+1] =   tmp;
    }
}

/*
 *  check_sac_nvhdr
 *
 *  Description: Determine the byte order of the SAC file
 *
 *  IN:
 *      const int nvhdr : nvhdr from header
 *
 *  Return:
 *      FALSE   no byte order swap is needed
 *      TRUE    byte order swap is needed
 *      -1      not in sac format ( nvhdr != SAC_HEADER_MAJOR_VERSION )
 *
 */
static int check_sac_nvhdr(const int nvhdr)
{
    int lswap = FALSE;

    if (nvhdr != SAC_HEADER_MAJOR_VERSION) {
        byte_swap((char*) &nvhdr, SAC_DATA_SIZEOF);
        if (nvhdr == SAC_HEADER_MAJOR_VERSION)
            lswap = TRUE;
        else
            lswap = -1;
    }
    return lswap;
}

/*
 *  map_chdr_in:
 *       map strings from buffer to memory
 */
static void map_chdr_in(char *memar, char *buff)
{
    char    *ptr1;
    char    *ptr2;
    int     i;

    ptr1 = memar;
    ptr2 = buff;

    memcpy(ptr1, ptr2, 8);
    *(ptr1+8) = '\0';
    ptr1 += 9;
    ptr2 += 8;

    memcpy(ptr1, ptr2, 16);
    *(ptr1+16) = '\0';
    ptr1 += 18;
    ptr2 += 16;

    for (i=0; i<21; i++) {
        memcpy(ptr1, ptr2, 8);
        *(ptr1+8) = '\0';
        ptr1 += 9;
        ptr2 += 8;
    }
}

/*
 *  read_head_in:
 *      read sac header in and deal with possible byte swap.
 *
 *  IN:
 *      const char *name : file name, only for debug
 *      SACHEAD    *hd   : header to be filled
 *      FILE       *strm : file handler
 *
 *  Return:
 *      0   :   Succeed and no byte swap
 *      1   :   Succeed and byte swap
 *     -1   :   fail.
 */
static int read_head_in(const char *name, SACHEAD *hd, FILE *strm)
{
    char   *buffer;
    int     lswap;

    if (sizeof(float) != SAC_DATA_SIZEOF || sizeof(int) != SAC_DATA_SIZEOF) {
        fprintf(stderr, "Mismatch in size of basic data type!\n");
        return -1;
    }

    /* read numeric parts of the SAC header */
    if (fread(hd, SAC_HEADER_NUMBERS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC header %s\n", name);
        return -1;
    }

    /* Check Header Version and Endian  */
    lswap = check_sac_nvhdr(hd->nvhdr);
    if (lswap == -1) {
        fprintf(stderr, "Warning: %s not in sac format.\n", name);
        return -1;
    } else if (lswap == TRUE) {
        byte_swap((char *)hd, SAC_HEADER_NUMBERS_SIZE);
    }

    /* read string parts of the SAC header */
    if ((buffer = (char *)malloc(SAC_HEADER_STRINGS_SIZE)) == NULL) {
        fprintf(stderr, "Error in allocating memory %s\n", name);
        return -1;
    }
    if (fread(buffer, SAC_HEADER_STRINGS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC header %s\n", name);
        free(buffer);
        return -1;
    }
    map_chdr_in((char *)(hd)+SAC_HEADER_NUMBERS_SIZE, buffer);
    free(buffer);

    return lswap;
}

/*
 *   map_chdr_out:
 *      map strings from memory to buffer
 */
static void map_chdr_out(char *memar, char *buff)
{
    char    *ptr1;
    char    *ptr2;
    int     i;

    ptr1 = memar;
    ptr2 = buff;

    memcpy(ptr2, ptr1, 8);
    ptr1 += 9;
    ptr2 += 8;

    memcpy(ptr2, ptr1, 16);
    ptr1 += 18;
    ptr2 += 16;

    for (i=0; i<21; i++) {
        memcpy(ptr2, ptr1, 8);
        ptr1 += 9;
        ptr2 += 8;
    }
}
/*
 *
 *  write_head_out
 *
 *  IN:
 *      const char *name : file name, only for debug
 *      SACHEAD     hd   : header to be written
 *      FILE       *strm : file handler
 *
 *  Return:
 *      -1  :   failed.
 *       0  :   success.
 */
static int write_head_out(const char *name, SACHEAD hd, FILE *strm)
{
    char *buffer;

    if (sizeof(float) != SAC_DATA_SIZEOF || sizeof(int) != SAC_DATA_SIZEOF) {
        fprintf(stderr, "Mismatch in size of basic data type!\n");
        return -1;
    }

    if (fwrite(&hd, SAC_HEADER_NUMBERS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in writing SAC data for writing %s\n", name);
        return -1;
    }

    if ((buffer = (char *)malloc(SAC_HEADER_STRINGS_SIZE)) == NULL){
        fprintf(stderr, "Error in allocating memory %s\n", name);
        return -1;
    }
    map_chdr_out((char *)(&hd)+SAC_HEADER_NUMBERS_SIZE, buffer);

    if (fwrite(buffer, SAC_HEADER_STRINGS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in writing SAC data for writing %s\n", name);
        return -1;
    }
    free(buffer);

    return 0;
}

/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-Xuping Feng+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/********************************************************************************************************************************/
/*+++++++++++++++++++++++++++++++++calculate cross-correlation of two inputting SAC files+++++++++++++++++++++++++*/
void cor (char *sac1, char *sac2, float lag_time, char *sac_cor) {
    int i, npts, tap_npts, lag_npts;
    float *data1, *data2, *data_cor, scale;
    fftw_complex *in1, *in2, *out1, *out2, *cor_in, *cor_out;
    fftw_plan p1, p2, p3;
    SACHEAD hd1, hd2;

    data1 = read_sac(sac1, &hd1);
    data2 = read_sac(sac2, &hd2);
    npts = 2*hd1.npts-1;
    lag_npts = (int)(lag_time/hd1.delta);
    tap_npts = (int)(lag_npts/50.);
    scale = (float)npts;

    /*-------------------------allocating dynamic memories to execute FFT transforam--------------------*/
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    cor_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    cor_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    data_cor = (float*) malloc(sizeof(float) * (2*lag_npts-1));

    for (i = 0; i < npts; i ++) {
        if ( i < hd1.npts ) {
            in1[i][0] = data1[hd1.npts-i-1];
            in2[i][0] = data2[hd2.npts-i-1];
        }
        else {
            in1[i][0] = data1[npts-i-1];
            in2[i][0] = data2[npts-i-1];
        }
        in1[i][1] = 0.;
        in2[i][1] = 0.;
    }

    p1 = fftw_plan_dft_1d(npts, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_1d(npts, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p1); fftw_execute(p2);
    fftw_destroy_plan(p1); fftw_destroy_plan(p2);

    /*-------------------------executing cross-correlation in frequency domain---------------------------*/
    for ( i = 0; i < npts; i ++ ) {
        cor_in[i][0] = out1[i][0]*out2[i][0] + out1[i][1]*out2[i][1];
        cor_in[i][1] = out1[i][1]*out2[i][0] - out1[i][0]*out2[i][1];
    }

    fftw_free(in1); fftw_free(in2); fftw_free(out1); fftw_free(out2);

    p3 = fftw_plan_dft_1d(npts, cor_in, cor_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p3);
    fftw_destroy_plan(p3);

    /*----------------------------to get cross-correlation according to lag time length------------------*/
    for ( i = 0; i < (2*lag_npts-1); i ++ ) data_cor[i] = cor_out[i+npts/2-lag_npts][0]/scale;
    hd1.npts = 2*lag_npts - 1; hd1.b = 0. - lag_time; hd1.e = lag_time;
    hd1.evlo = hd1.stlo; hd1.evla = hd1.stla; hd1.stlo = hd2.stlo; hd1.stla = hd2.stla;

    /*------------------------------------write cross-correlation file-----------------------------------*/
    write_sac(sac_cor, hd1, data_cor);
    fftw_free(cor_in); fftw_free(cor_out); free(data1); free(data2); free(data_cor);
}


/*++++++++++++++++++++++++++++++++++++++++++++++normalization in time domain++++++++++++++++++++++++++++++++++++++++*/
void norm( char *sacin, char *sacout, int npts ) {
    FILE *ff = fopen( "norm.sh", "w" );
    fprintf(ff, "${SAC_DISPLAY_COPYRIGHT}=0\n");
    fprintf(ff, "sac<<END\n");
    fprintf(ff, "r %s\n", sacin);
    fprintf(ff, "abs\n");
    fprintf(ff, "smooth mean h %d\n", npts);
    fprintf(ff, "w a.avg\n");
    fprintf(ff, "r %s\n", sacin);
    fprintf(ff, "divf a.avg\n");
    fprintf(ff, "w %s\n", sacout);
    fprintf(ff, "q\n"); fprintf(ff, "END");
    fclose(ff); system("sh norm.sh"); system("rm norm.sh");
}

/*++++++++++++++++++++++++++++++++++++++++++++++normalization in time domain++++++++++++++++++++++++++++++++++++++++*/
void normal( char *sacin, char *sacout, int npts ) {
    int i;
    float tmp = 0, *data, *mean;
    SACHEAD hd;
    data = read_sac(sacin, &hd);
    mean = (float *) malloc( sizeof(float) * hd.npts );

    for ( i = 0; i < (2*npts+1); i ++ ) tmp += fabs(data[i]) / (2*npts+1);

/*-----------------------------------run absolute mean smooth in middle part-----------------------------------------*/

    for ( i = npts; i < (hd.npts-npts-1); i ++ ) {
        mean[i] = data[i] / tmp;
        tmp = tmp - fabs(data[i-npts])/(2*npts+1) + fabs(data[i+npts+1])/(2*npts+1);
    }
/*------------------------------------run absolute mean smooth on each side------------------------------------------*/

    for ( i = 0; i < npts; i ++ ) {
        mean[i] = data[i] / tmp;
        mean[hd.npts-npts+i] = data[hd.npts-npts+i] / tmp;
    }
    write_sac(sacout, hd, mean);
    free(data); free(mean);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++get 2's integral power+++++++++++++++++++++++++++++++++++++++++*/
int pow_next2( int n ) {
    int m;
    float f;
    f = log((float)n) / log(2.);
    if ( f == (int)f ) return n;
    else {
        m = (int)(pow(2, (int)f + 1));
        return m;
    }
}
void bp ( char *sacin, char *sacout, float f1, float f2, float f3, float f4, int npow ) {
    int i, j;
    float *datain, *dataout, sp, f, pi = 3.1415926535;
    fftw_complex *in, *out;
    fftw_plan p1, p2;
    SACHEAD hd;

    datain = read_sac( sacin, &hd );
    sp = 1./hd.delta;

    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * hd.npts);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * hd.npts);
    dataout = (float *) malloc(sizeof(float) * hd.npts);

    for ( i = 0; i < hd.npts; i ++ ) dataout[i] = 0.;
    for ( i = 0; i < hd.npts; i ++ ) {
        f = i*sp/hd.npts;
        if ( f < f1 ) continue;
        else if ( f >= f1 && f < f2 ) {
            dataout[i] = 1.;
            for ( j = 0; j < npow; j ++ )
                dataout[i] = dataout[i]*(1.+sin(pi/2.*(f-f1)/(f2-f1))) / 2.;
        }
        else if ( f >= f2 && f < f3 ) dataout[i] = 1.;
        else if ( f >= f3 && f <= f4 ) {
            dataout[i] = 1.;
            for ( j = 0; j < npow; j ++ )
                dataout[i] = dataout[i]*(1.+cos(pi/2.*(f-f3)/(f4-f3))) / 2.;
        }
        else continue;
    }

    for ( i = 0; i < hd.npts; i ++ ) {
        in[i][0] = datain[i];
        in[i][1] = 0.;
    }
    p1 = fftw_plan_dft_1d( hd.npts, in , out, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    for ( i = 0; i < hd.npts; i ++ ) {
        out[i][0] = out[i][0] * dataout[i];
        out[i][1] = out[i][1] * dataout[i];
    }

    p2 = fftw_plan_dft_1d( hd.npts, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p2);
    fftw_destroy_plan(p2);
    for ( i = 0; i < hd.npts; i ++ ) dataout[i] = in[i][0]/hd.npts;

    write_sac(sacout, hd, dataout);
    fftw_free(in); fftw_free(out); free(datain); free(dataout);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++spectral whitening+++++++++++++++++++++++++++++++++++++++++++++++*/
void whiten_f ( char *sacin, char *sacout, int npts, float f1, float f2, float f3, float f4 ) {
    float *data, *sqr, *sout, sum = 0, sp, f;
    int i, j, k, index1, index2;
    fftw_complex *in, *out;
    fftw_plan p1, p2;
    SACHEAD hd;

    data = read_sac( sacin, &hd );
    sp = 1./hd.delta;

    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * hd.npts);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * hd.npts);
    sqr = (float *) malloc(sizeof(float) * hd.npts );
    sout = (float *) malloc(sizeof(float) * hd.npts );

    for ( i = 0; i < hd.npts; i ++ ) {
        in[i][0] = data[i]; in[i][1] = 0.; sout[i] = 0.;
    }
    p1 = fftw_plan_dft_1d( hd.npts, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(p1);
    fftw_destroy_plan(p1);

    for ( i = 0; i < hd.npts; i ++ ) sqr[i] = sqrt( pow(out[i][0],2.) + pow(out[i][1],2.) );

    for ( i = 0; i < hd.npts; i ++ ) {
        f = i * sp/hd.npts;
        if ( f >= f1 && f <= f4 ) {
            sum = 0.;
            for ( j = -npts; j <= npts; j ++ ) {
                k = i + j;
                sum += sqr[k];
            }
            sout[i] = sum/(2*npts+1);
        }
        else sout[i] = sqr[i];
    }

    for ( i = 0; i < hd.npts; i ++ ) {
        f = i * sp / hd.npts;
        if ( f >= f1 && f <= f4 ) {
            out[i][0] /= sout[i]; out[i][1] /= sout[i];
        }
        else {
            out[i][0] = 0.; out[i][1] = 0.;
        }
    }

    p2 = fftw_plan_dft_1d( hd.npts, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p2);
    fftw_destroy_plan(p2);
    for ( i = 0; i < hd.npts; i ++ ) data[i] = in[i][0]/hd.npts;
    fftw_free(in); fftw_free(out);

    write_sac( sacout, hd, data );
    free(data);
}

/*++++++++++++++++++++++++++++++++++++++++compute julian day from year-month-day+++++++++++++++++++++++++++++++++*/
int julian( int year, int mon, int day ) {
    int jday;

    if ( (year % 100 == 0 && year % 400 == 0) || (year%100 != 0 && year % 4 == 0)  ) switch( mon ) {
        case 1: jday = day; break;
        case 2: jday = 31 + day; break;
        case 3: jday = 31 + 29 + day; break;
        case 4: jday = 31 + 29 + 31 + day; break;
        case 5: jday = 31 + 29 + 31 + 30 + day; break;
        case 6: jday = 31 + 29 + 31 + 30 + 31 + day; break;
        case 7: jday = 31 + 29 + 31 + 30 + 31 + 30 + day; break;
        case 8: jday = 31 + 29 + 31 + 30 + 31 + 30 + 31 + day; break;
        case 9: jday = 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + day; break;
        case 10: jday = 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + day; break;
        case 11: jday = 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + day; break;
        case 12: jday = 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + day; break;
        default: break;
    }
    else switch( mon ) {
        case 1: jday = day; break;
        case 2: jday = 31 + day; break;
        case 3: jday = 31 + 28 + day; break;
        case 4: jday = 31 + 28 + 31 + day; break;
        case 5: jday = 31 + 28 + 31 + 30 + day; break;
        case 6: jday = 31 + 28 + 31 + 30 + 31 + day; break;
        case 7: jday = 31 + 28 + 31 + 30 + 31 + 30 + day; break;
        case 8: jday = 31 + 28 + 31 + 30 + 31 + 30 + 31 + day; break;
        case 9: jday = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + day; break;
        case 10: jday = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + day; break;
        case 11: jday = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + day; break;
        case 12: jday = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + day; break;
        default: break;
    }

    return jday;
}

/*++++++++++++++++++++++++++++++calculate abslute time relative to 1970-01-01T00:00:00++++++++++++++++++++++++++++++++*/
float abs_time(int year, int jday, int hour, int min, int sec, float msec) {
    int yr, mn, nyday = 0;
    float abssec;

    for ( yr = 1971; yr < year; yr ++ ) {
        if ( 4*(yr/4) == yr ) nyday += 366;
        else nyday += 365;
    }
    abssec = 24. * 3600. * (nyday + jday - 1 ) + 3600. * hour + 60. * min + 1.* sec + 0.001 * msec;
    return abssec;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++cut SAC foramt file++++++++++++++++++++++++++++++++++++++++++++*/
void cut_sac(char *sacin, char *sacout, float evt0, float startt0, int npts) {
    int start_index, i;
    float *cut_data, *data, sact0;
    SACHEAD hd;
    cut_data = (float *) malloc( sizeof(float) * npts );
    data = read_sac(sacin, &hd);
    sact0 = abs_time(hd.nzyear, hd.nzjday, hd.nzhour, hd.nzmin, hd.nzsec, hd.nzmsec);
    start_index = (int) ( (evt0 - sact0 + startt0 )/hd.delta);
    for (i = 0; i < npts; i ++) cut_data[i] = data[i+start_index];
    hd.b = 0.; hd.e = (npts-1) * hd.delta; hd.npts = npts;
    write_sac(sacout, hd, cut_data);
    free(data); free(cut_data);
}

/*+++++++++++++++++++++++++Spectral whitening: number of FFT points is 2^n(n is an integer)+++++++++++++++++++++++++*/
void spe_whi ( char *sacin, char *sacout, int npts, float f1, float f2, float f3, float f4 ) {
    float *data, *sqr, *sout, sum = 0, sp, f;
    int i, j, k, f1_index, f4_index, fftn;
    fftw_complex *in, *out;
    fftw_plan p;
    SACHEAD hd;

    data = read_sac( sacin, &hd );
    sp = 1./hd.delta;
    fftn = pow_next2(hd.npts);
    f1_index = (int)(f1*fftn*hd.delta); f4_index = (int)(f4*fftn*hd.delta);

    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftn);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftn);
    sqr = (float *) malloc(sizeof(float) * hd.npts );
    sout = (float *) malloc(sizeof(float) * hd.npts );

    for ( i = 0; i < fftn; i ++ ) {
        if ( i < hd.npts ) in[i][0] = data[i];
        else in[i][0] = 0.; in[i][1] = 0.;
    }
    p = fftw_plan_dft_1d( fftn, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(p);
    fftw_destroy_plan(p);

    for ( i = 0; i < hd.npts; i ++ ) {
        sqr[i] = sqrt( pow(out[i][0],2.) + pow(out[i][1],2.) );
        sout[i] = 0.;
    }
    for ( i = (f1_index - npts); i < (f1_index + npts); i ++ ) sum += sqr[i];
    for ( i = f1_index; i <= f4_index; i ++ ) {
        sout[i] = sum/(2*npts+1);
        sum = sum + sqr[i+npts] - sqr[i-npts];
    }

    for ( i = 0; i < fftn; i ++ ) {
        if ( i >= f1_index && i <= f4_index ) {
            out[i][0] /= sout[i]; out[i][1] /= sout[i];
        }
        else {
            out[i][0] = 0.; out[i][1] = 0.;
        }
    }

    p = fftw_plan_dft_1d( fftn, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p);
    fftw_destroy_plan(p);
    for ( i = 0; i < hd.npts; i ++ ) data[i] = in[i][0]/fftn;
    fftw_free(in); fftw_free(out);

    write_sac( sacout, hd, data );
    free(data);
}
