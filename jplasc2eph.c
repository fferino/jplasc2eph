/*
 jplasc2eph.c


 The MIT License (MIT)

 Copyright (c) 2015 Fabrice Ferino

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>

#define OFFSET_IN_FILE(f)  k##f##Offset
#define SIZE_IN_FILE(f) k##f##Size

#define DECLARE_SIZE_OFFSET( f, s , p)  \
 const long OFFSET_IN_FILE(f) = OFFSET_IN_FILE(p) + SIZE_IN_FILE(p); \
 const long SIZE_IN_FILE(f) = s

const long OFFSET_IN_FILE(TTL) = 0;
const long SIZE_IN_FILE(TTL) = 3 * 84;
DECLARE_SIZE_OFFSET(CName1, 400*6, TTL);
DECLARE_SIZE_OFFSET(SS, 3 * sizeof(double), CName1);
DECLARE_SIZE_OFFSET(NCon, sizeof(int), SS);
DECLARE_SIZE_OFFSET(AU, sizeof(double), NCon);
DECLARE_SIZE_OFFSET(EmRatio, sizeof(double), AU);
DECLARE_SIZE_OFFSET(IPT, 3 * 12 * sizeof(int), EmRatio);
DECLARE_SIZE_OFFSET(Denum, sizeof(int), IPT);
DECLARE_SIZE_OFFSET(LPT, 3 * sizeof(int), Denum);
DECLARE_SIZE_OFFSET(CName2, 600 * 6, LPT);

enum
{
    ERR_EOF = -1,
    ERR_FILE_CREATE = -2,
    ERR_FILE_OPEN = -3,
    ERR_FORMAT = -4,
    ERR_FILE_SEEK = -5,
    ERR_FILE_READ = -6,
    ERR_DATA_GAP = -7,
    ERR_FILE_WRITE = -8,
    ERR_INVALID_DATA_RANGE = -9,
    ERR_USAGE = -100
};

enum
{
    RECORD_CONTINUATION = 0, /* record is continuation of prev one => write */
    RECORD_DUPLICATION = 1,    /* record is duplicate of prev one => skip (files overlap) */
    RECORD_NOT_CHECKED = 2,
    RECORD_OUT_OF_SEQ = 3,   /* record is out of sequence with prev one => abort*/
};

/* type definitions */
typedef struct EphCtx
{
    int ksize;
    int nCoeff;
    double ss[3];
    int nCon;
    int denumIndex, auIndex, emRatioIndex;
    double firstJD; /* start JD */
    double lastJD; /* last JD */
} EphCtx;

/* globals */
char gLine[128]; /* current line being read from file */


/*---------------------------------------------------------------------------*/

void usage()
{
    printf("usage: jplasc2eph header_file [ascii_file ...]\n"
           "\tconvert JPL ASCII ephemerides files to binary\n");
}

/*---------------------------------------------------------------------------*/

int cmp( const void* a, const void* b)
{
    const char** sa = (const char**) a;
    const char** sb = (const char**) b;

    return strcmp(*sa,*sb);
}

/*---------------------------------------------------------------------------*/

long RecordSize( const EphCtx* ctx)
{
    return ctx->nCoeff * sizeof(double);
}

/*---------------------------------------------------------------------------*/

int GetNextGroup( FILE* f, int* groupNo)
{
    for (;;)
    {
        if (!fgets( gLine, sizeof(gLine), f))
        {
            return ERR_EOF;
        }

        if (1 == sscanf(gLine, "GROUP %d", groupNo))
        {
            return fgets( gLine, sizeof(gLine), f) ? 0 : ERR_FORMAT;
        }
    }
}

/*---------------------------------------------------------------------------*/

int GetSpecificGroup( FILE* f, int groupNo)
{
    int nextGroup;

    int retVal = GetNextGroup(f, &nextGroup);

    if (0 == retVal)
    {
        retVal = nextGroup == groupNo ? 0 : ERR_FORMAT;
    }
    return retVal;
}

/*---------------------------------------------------------------------------*/

int ScanNextLine(FILE* f, int numArgs, const char* format, ...)
{
    va_list args;
    int numRead;
    if ( !fgets(gLine, sizeof(gLine), f))
    {
        return ERR_EOF;
    }

    va_start(args, format);

    numRead = vsscanf(gLine, format, args);

    va_end(args);

    return (numArgs == numRead) ? 0 : ERR_FORMAT;
}

/*---------------------------------------------------------------------------*/

int ScanNextValuesLine(FILE* f, double values[/*3*/])
{
    int numRead;
    char* s;
    if ( !fgets(gLine, sizeof(gLine), f))
    {
        return ERR_EOF;
    }

    s = gLine;

    while (*s)
    {
        if (*s == 'D') { *s = 'E'; }
        ++s;
    }

    numRead = sscanf(gLine, " %lf %lf %lf",
                     values, values+1, values+2);

    return (3 == numRead) ? 0 : ERR_FORMAT;
}

/*---------------------------------------------------------------------------*/

int ProcessNamesAux( FILE* out, FILE* f, EphCtx* ctx, int pass,
                    int numNames)
{
    size_t len;
    int n;
    int currNameIndex = (0 == pass) ? 0 : 400;

    while (numNames > 0)
    {
        if ( !fgets(gLine, sizeof(gLine), f))
        {
            return ERR_EOF;
        }

        len = strlen(gLine);
        if (len != 81)
        {
            return ERR_FORMAT;
        }

        for (n = 0; n < 10; --numNames, ++n, ++currNameIndex)
        {
            const char* name = gLine + n * 8 + 2;

            if (0 == memcmp("AU    ", name, 6))
            {
                ctx->auIndex = currNameIndex;
            }
            else if (0 == memcmp("EMRAT ", name, 6))
            {
                ctx->emRatioIndex = currNameIndex;
            }
            else if (0 == memcmp("DENUM ", name, 6))
            {
                ctx->denumIndex = currNameIndex;
            }
            if (1 != fwrite( gLine + n * 8 + 2, 6, 1, out))
            {
                return ERR_FILE_WRITE;
            }
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------*/

/* read ctx->ncon names, names 6 characters long (padded with space)
 and preceded by 2 spaces  10 names per line of 80 chars*/
int ProcessNames(FILE* out, FILE* f, EphCtx* ctx)
{
    int retVal;

    /* sentinels to make sure we found these */
    ctx->denumIndex = ctx->emRatioIndex = ctx->auIndex = -1;

    int numNames = (400 > ctx->nCon) ? ctx->nCon : 400; /* maximum 400 names (2400 chars)
                                                         in first section */
    if ((retVal = ProcessNamesAux(out, f, ctx, 0, numNames)))
    {
        return retVal;
    }

    /* still some more names ?*/
    numNames = ctx->nCon - 400;

    /* if more than 400 ncon, write them after the rest */
    if (numNames > 0)
    {
        if (0 != fseek( out, OFFSET_IN_FILE(CName2), SEEK_SET))
        {
            return ERR_FILE_SEEK;
        }

        if (( retVal = ProcessNamesAux(out, f, ctx, 1, numNames)))
        {
            return retVal;
        }
    }

    if (-1 == ctx->auIndex  || -1 == ctx->denumIndex || -1 == ctx->emRatioIndex)
    {
        return ERR_FORMAT;
    }

    return 0;
}


/*---------------------------------------------------------------------------*/

int CheckWithPrevious( double prev1stVal[/*2*/], double new1stVal[/*2*/])
{
    if (prev1stVal[0] == prev1stVal[1])
    {
        /* nothing to check for */
        return RECORD_CONTINUATION;
    }

    if( prev1stVal[0] == new1stVal[0] && prev1stVal[1] == new1stVal[1])
    {
        return RECORD_DUPLICATION;
    }
    if (prev1stVal[1] == new1stVal[0])
    {
        return RECORD_CONTINUATION;
    }
    return RECORD_OUT_OF_SEQ;
}

/*---------------------------------------------------------------------------*/

int WriteRecord(FILE* out, FILE* f, int numValues, double prev1stValues[/*2*/])
{
    int i, retVal;
    int numLines = (numValues + 3) / 3;
    int numValuesLastLine = numValues % 3;

    /* 6: 9 / 3 = 3 lines; 6 % 3 = 0 values on the last line */
    /* 7: 10 / 3 = 3 lines; 7 % 3 = 1 value on the last line */
    /* 8: 11 / 3 = 3 lines; 8 % 3 = 2 value on the last line */
    /* 9: 12 / 3 = 4 lines; 9 % 3 = 0 values on the last line */

    double values[3];

    if (0 == numLines)
    {
        return 0;
    }

    /* read the first line separately -- 
     need to check and return the first two vals */
    if ((retVal = ScanNextValuesLine(f, values)))
    {
        return retVal;
    }

    if (prev1stValues)
    {
        retVal = CheckWithPrevious(prev1stValues, values);
    }
    else
    {
        retVal = RECORD_NOT_CHECKED;
    }

    switch (retVal)
    {
        case RECORD_CONTINUATION:
            /* update the prev values */
            prev1stValues[0] = values[0];
            prev1stValues[1] = values[1];
            /* fall through: */

        case RECORD_NOT_CHECKED:
            if (3 != fwrite(values, sizeof(double), 3, out))
            {
                return ERR_FILE_WRITE;
            }
            /* write the rest of the record */
            for (i = 1; i < numLines-1; ++i)
            {
                if ((retVal = ScanNextValuesLine(f, values)))
                {
                    return retVal;
                }
                if (3 != fwrite(values, sizeof(double), 3, out))
                {
                    return ERR_FILE_WRITE;
                }
            }
            /* last line:  write 0 to 2 values */
            if (numValuesLastLine > 0)
            {
                if ((retVal = ScanNextValuesLine(f, values)))
                {
                    return retVal;
                }
                if (numValuesLastLine != fwrite(values, sizeof(double),
                                                numValuesLastLine, out))
                {
                    return ERR_FILE_WRITE;
                }
            }
            break;

        case RECORD_DUPLICATION:
            /* no need to write anything */
            /* skip over the rest of the record */
            for (i = 1; i < numLines; ++i)
            {
                if (!fgets(gLine, sizeof(gLine), f))
                {
                    return ERR_EOF;
                }
            }
            break;

        default:
            return ERR_DATA_GAP;
    }
    return 0;
}

/*---------------------------------------------------------------------------*/

int ProcessGroup1050(FILE* out, FILE* f)
{
    int i, j;
    int ipt[12][3]; /* match fortran order, column first */
    int lpt[3];
    char* s;

    for (i = 0; i < 3; ++i)
    {
        if (!fgets(gLine, sizeof(gLine), f))
        {
            return ERR_EOF;
        }

        s = gLine;
        for (j = 0; *s && j < 12; ++j)
        {
            ipt[j][i] = (int) strtol(s, &s, 10);
        }
        if ( !*s)
        {
            return ERR_FORMAT;
        }
        else
        {
            lpt[i] = (int) strtol(s, NULL, 10);
        }
    }

    if (0 != fseek(out, OFFSET_IN_FILE(IPT), SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }
    if (3 * 12 != fwrite(ipt, sizeof(int), 3 * 12, out))
    {
        return ERR_FILE_WRITE;
    }

    if ( 0 != fseek(out, OFFSET_IN_FILE(LPT), SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }
    if (3 != fwrite(lpt, sizeof(int), 3, out))
    {
        return ERR_FILE_WRITE;
    }

    return 0;
}

/*---------------------------------------------------------------------------*/

int ProcessHeader(FILE* out, EphCtx* ctx, const char* headerFileName)
{
    int i, tmp, retVal;

    FILE* f = fopen(headerFileName, "r");
    if (!f)
    {
        fprintf(stderr, "Error:unable to open header file: %s\n",
                headerFileName);
        return ERR_FILE_OPEN;

    }

    if ((retVal = ScanNextLine(f, 2, "KSIZE= %d NCOEFF= %d",
                              &ctx->ksize, &ctx->nCoeff)))
    {
        goto exit;
    }

    /*************** GROUP 1010 *****************/
    if ((retVal = GetSpecificGroup(f, 1010)))
    {
        goto exit;
    }

    /* goofy fortran i/o: fill 3 lines of 84 chars from the file */
    for (i = 0; i < 3; ++i)
    {
        long numChars;

        if ( !fgets(gLine, sizeof(gLine), f))
        {
            retVal = ERR_EOF;
            goto exit;
        }

        numChars = strlen( gLine);

        while (iscntrl(gLine[numChars-1]))
        {
            --numChars;
        }

        memset( gLine + numChars, 0x20, 84 - numChars);
        if (84 != fwrite( gLine, 1, 84, out))
        {
            retVal = ERR_FILE_WRITE;
            goto exit;
        }
    }

    /*************** GROUP 1030 *****************/
    if ((retVal = GetSpecificGroup(f, 1030)))
    {
        goto exit;
    }

    if ((retVal = ScanNextLine(f, 3, "%lf %lf %lf",
                               ctx->ss, ctx->ss+1, ctx->ss+2)))
    {
        goto exit;
    }

    /*************** GROUP 1040 *****************/
    if ((retVal = GetSpecificGroup(f, 1040)))
    {
        goto exit;
    }

    if ((retVal = ScanNextLine(f, 1, "%d", &ctx->nCon)))
    {
        goto exit;
    }

    /* read ctx->ncon names, names 6 characters long (padded with space)
     and preceded by 2 spaces  10 names per line of 80 chars*/
    if ((retVal = ProcessNames(out, f, ctx)))
    {
        goto exit;
    }

    /*************** GROUP 1041 *****************/
    if ((retVal = GetSpecificGroup(f, 1041)))
    {
        goto exit;
    }

    if ((retVal = ScanNextLine(f, 1, "%d", &tmp)))
    {
        goto exit;
    }

    if (tmp != ctx->nCon)
    {
        retVal = ERR_FORMAT;
        goto exit;
    }
    /* go to data record 1 and write the values */
    if (0 != fseek(out, RecordSize(ctx), SEEK_SET))
    {
        retVal = ERR_FILE_SEEK;
        goto exit;
    }

    if ((retVal = WriteRecord(out, f, ctx->nCon, NULL)))
    {
        goto exit;
    }

    /*************** GROUP 1050 *****************/
    if ((retVal = GetSpecificGroup(f, 1050)))
    {
        goto exit;
    }

    if ((retVal = ProcessGroup1050(out, f)))
    {
        goto exit;
    }


    retVal = 0;

exit:

    fclose(f);
    return retVal;

}

/*---------------------------------------------------------------------------*/

int AppendAscii( FILE* out, EphCtx* ctx, int seq,
                const char* asciiFileName,
                double lastJD[/*2*/])
{
    int retVal = 0, recNo, numCoeff;

    FILE* f = fopen(asciiFileName, "r");
    if (!f)
    {
        fprintf(stderr, "Error:unable to open data file: %s\n",
                asciiFileName);
        return ERR_FILE_OPEN;
    }

    for (;;)
    {
        if ((retVal = ScanNextLine(f, 2, "%d %d", &recNo, &numCoeff)))
        {
            if (ERR_EOF == retVal)
            {
                retVal = 0; /* normal exit: read everything that was in the file*/
            }
            goto exit;
        }
        if (numCoeff != ctx->nCoeff)
        {
            retVal = ERR_FORMAT;
            goto exit;
        }
        if ((retVal = WriteRecord(out, f, numCoeff, lastJD)))
        {
            goto exit;
        }
    }

exit:
    fclose(f);
    return retVal;
}

/*---------------------------------------------------------------------------*/

int ReadDouble(FILE* f, long offset, double* val)
{
    if (0 != fseek(f, offset, SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }

    if (1 != fread(val, sizeof(double), 1, f))
    {
        return ERR_FILE_READ;
    }
    return 0;
}

/*---------------------------------------------------------------------------*/

int CopyDouble( FILE* f, long srcOffset, long destOffset)
{
    int retVal;
    double value;

    if (( retVal = ReadDouble(f, srcOffset, &value)))
    {
        return retVal;
    }

    if (0!=fseek(f, destOffset, SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }

    if (1 != fwrite(&value, sizeof(double), 1, f))
    {
        return ERR_FILE_WRITE;
    }
    return 0;
}


/*---------------------------------------------------------------------------*/

int CopyDouble2Int( FILE* f, long srcOffset, long destOffset)
{
    int retVal;
    double value;
    int value2;

    if (( retVal = ReadDouble(f, srcOffset, &value)))
    {
        return retVal;
    }

    value2 = (int) floor(value);

    if (0!=fseek(f, destOffset, SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }

    if (1 != fwrite(&value2, sizeof(int), 1, f))
    {
        return ERR_FILE_WRITE;
    }
    return 0;
}


/*---------------------------------------------------------------------------*/

int WriteRecord0( FILE* f, EphCtx* ctx, long recordSize)
{
    int retVal;

    /*** SS ***/
    if (0 != fseek(f, OFFSET_IN_FILE(SS), SEEK_SET))
    {
        return ERR_FILE_SEEK;
    }

    if (3 != fwrite(ctx->ss, sizeof(double), 3, f))
    {
        return ERR_FILE_WRITE;
    }
    /*** ncon ****/
    if ( 1 != fwrite(&ctx->nCon, sizeof(int), 1, f))
    {
    }
    
    /*** AU ***/
    if ((retVal = CopyDouble( f, recordSize + ctx->auIndex * sizeof(double),
                                OFFSET_IN_FILE(AU))))
    {
        return retVal;
    }

    /*** EMRATIO ***/
    if ((retVal = CopyDouble( f, recordSize + ctx->emRatioIndex * sizeof(double),
                             OFFSET_IN_FILE(EmRatio))))
    {
        return retVal;
    }

    /*** DENUM ***/
    if ((retVal = CopyDouble2Int( f, recordSize + ctx->denumIndex * sizeof(double),
                                  OFFSET_IN_FILE(Denum))))
    {
        return retVal;
    }

    return 0;
}

/*---------------------------------------------------------------------------*/

int asc2eph( const char* headerFile, int numAsciiFiles, const char* asciiFiles[])
{
    FILE* f = fopen( "JPLEPH", "w+");
    int i, retVal = 0;
    EphCtx ctx;
    long recordSize;
    double previousValues[2] = {0.0, 0.0};

    if (!f)
    {
        fprintf(stderr, "Error:unable to create binary ephemeris file\n");
        return 2;
    }

    retVal = ProcessHeader( f, &ctx, headerFile);
    if (retVal)
    {
        goto exit;
    }

    /* make sure files are processed in order */
    if (numAsciiFiles > 1)
    {
        qsort( asciiFiles, numAsciiFiles, sizeof(asciiFiles[0]), cmp);
    }

    recordSize = RecordSize(&ctx);
    /* go to data record 2  -- data will go there */
    if (0 != fseek(f, 2 * recordSize, SEEK_SET))
    {
        retVal = ERR_FILE_SEEK;
        goto exit;
    }

    for (i = 0; i < numAsciiFiles; ++i)
    {
        retVal = AppendAscii( f, &ctx, i, asciiFiles[i], previousValues);
        if (retVal)
        {
            goto exit;
        }
    }

    /* WRAP EVERYTHING by writing missing part of record 0: SS,  NCon,
     AU, EM_RATIO, DENUM */
    /* first JD: read from first data record, i.e. record 2 */
    if (0 != fseek(f, 2 * recordSize, SEEK_SET))
    {
        retVal = ERR_FILE_SEEK;
        goto exit;
    }
    if (1 != fread(previousValues, sizeof(double), 1, f))
    {
        retVal = ERR_FILE_READ;
        goto exit;
    }
    /* previousValues[0] = firstJD, previousValues[1] = lastJD */

    /* make sure minJD <= firstJD and lastJD <= maxJD  
      firstJD <= lastJD is checked when processing the ascii files */
    if (previousValues[0] < ctx.ss[0] || previousValues[1] > ctx.ss[1])
    {
        retVal = ERR_INVALID_DATA_RANGE;
        goto exit;
    }
    ctx.ss[0] = previousValues[0];
    ctx.ss[1] = previousValues[1];


    retVal = WriteRecord0(f, &ctx, recordSize);

exit:

    fclose(f);
    return retVal;
}

/*---------------------------------------------------------------------------*/

int main(int argc, const char * argv[])
{
    if (argc < 3)
    {
        usage();
        return 1;
    }

    return asc2eph( argv[1], argc-2, argv+2);
}
