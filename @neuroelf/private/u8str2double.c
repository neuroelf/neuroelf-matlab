/* u8str2double  - convert tabular string (from file) to double matrix

FORMAT:       numtab = u8str2double(asctab)

Input fields:

      asctab      ASCII representation of table

Output fields:

      numtab      double matrix

% Version:  v0.9b
% Build:    11062416
% Date:     Jun-15 2010, 12:37 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, Jochen Weber
All rights reserved.

Actual conversion function (strtodmod, see below) originally
Copyright (c) 2002, Michael Ringgaard
Modified version
Copyright (c) 2002, Michael Ringgaard
Copyright (c) 2010, Jochen Weber
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Columbia University nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "mex.h"
#include "math.h"
#include "ctype.h"
#include <float.h>

#define SSTR (*sstr & 0xff)

/* array to detect any valid digit character:
   - '+' (0x2b) (+ := sign character)
   - '-' (0x2d) (- := sign character)
   - '.' (0x2e) (. := decimal point)
   - '0' (0x30) (0 := digit)
   - '1' (0x31) (1 := digit)
   - '2' (0x32) (2 := digit)
   - '3' (0x33) (3 := digit)
   - '4' (0x34) (4 := digit)
   - '5' (0x35) (5 := digit)
   - '6' (0x36) (6 := digit)
   - '7' (0x37) (7 := digit)
   - '8' (0x38) (8 := digit)
   - '9' (0x39) (9 := digit)
   - 'E' (0x45) (E := exponent marker)
   - 'e' (0x65) (e := exponent marker) */
const bool isanydigit[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x00 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x20 */
    0, 0, 0, 1, 0, 1, 1, 0,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, /* 0x40 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, /* 0x60 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x80 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xa0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xc0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xe0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0};

/* array to detect any valid space character:
   - \0  (0x00) (null-terminator := space)
   - \t  (0x09) (horiz. tab := space)
   - \n  (0x0a) (new-line := linebreak)
   - \b  (0x0b) (vert. tab := linebreak)
   - \f  (0x0c) (form-feed := linebreak)
   - \r  (0x0d) (carriage-return := linebreak)
   - ' ' (0x20) (space character := space)
   - ',' (0x2c) (comma := space)
   - ';' (0x3b) (semicolon := linebreak)
   - '[' (0x5b) (opening bracket := space)
   - ']' (0x5d) (closing bracket := space)
   - '|' (0x7c) (vertical line := space)
   - ' ' (0xa0) (odd implementation of special space := space) */
const bool isanyspace[256] = {
    1, 0, 0, 0, 0, 0, 0, 0, /* 0x00 */
    0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, /* 0x20 */
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x40 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x60 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x80 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, /* 0xa0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xc0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xe0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0};

/* array to detect any valid digit character:
   - '+' (0x2b) (0 := digit)
   - '-' (0x2d) (0 := digit)
   - '.' (0x2e) (0 := digit)
   - '0' (0x30) (0 := digit)
   - '1' (0x31) (0 := digit)
   - '2' (0x32) (0 := digit)
   - '3' (0x33) (0 := digit)
   - '4' (0x34) (0 := digit)
   - '5' (0x35) (0 := digit)
   - '6' (0x36) (0 := digit)
   - '7' (0x37) (0 := digit)
   - '8' (0x38) (0 := digit)
   - '9' (0x39) (0 := digit) */
const bool isfirstdigit[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x00 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x20 */
    0, 0, 0, 1, 0, 1, 1, 0,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x40 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x60 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x80 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xa0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xc0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xe0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0};

/* array to detect valid linebreaks:
   - \n  (0x0a) (new-line := linebreak)
   - \b  (0x0b) (vert. tab := linebreak)
   - \f  (0x0c) (form-feed := linebreak)
   - \r  (0x0d) (carriage-return := linebreak)
   - ';' (0x3b) (semicolon := linebreak) */
const bool islinebreak[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x00 */
    0, 0, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x20 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x40 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x60 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x80 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xa0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xc0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xe0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0};

/* array to detect valid space characters:
   - \0  (0x00) (null-terminator := space)
   - \t  (0x09) (horiz. tab := space)
   - ' ' (0x20) (space character := space)
   - ',' (0x2c) (comma := space)
   - '[' (0x5b) (opening bracket := space)
   - ']' (0x5d) (closing bracket := space)
   - '|' (0x7c) (vertical line := space)
   - ' ' (0xa0) (odd implementation of special space := space) */
const bool isspacev[256] = {
    1, 0, 0, 0, 0, 0, 0, 0, /* 0x00 */
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, /* 0x20 */
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x40 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x60 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0x80 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, /* 0xa0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xc0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, /* 0xe0 */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0};

/* helper function strtodmod

   taken from http://www.jbox.dk/sanos/source/lib/strtod.c.html

   modified to work without second parameter

   strtod.c

   Convert string to double

 Copyright (C) 2002 Michael Ringgaard. All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. Neither the name of the project nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 SUCH DAMAGE.
 */
int strtodmod(const unsigned char *str, int ne, double *storage)
{
    double number;
    int exponent;
    int negative;
    double p10;
    int n;
    int num_digits;
    int num_decimals;

    /* handle optional sign */
    negative = 0;
    switch (*str) {
        case '-': negative = 1; /* fall through to increment position */
        case '+': --ne; ++str;
    }

    number = 0.;
    num_digits = 0;

    /* process string of digits */
    while (ne > 0 && isdigit(*str)) {
        number = number * 10. + (*str++ - '0');
        --ne;
        ++num_digits;
    }

    exponent = 0;
    num_decimals = 0;

    /* already reached the end? */
    if (ne == 0) {
        if (negative)
            *storage = -number;
        else
            *storage = number;
        return ne;
    }

    /* process decimal part */
    if (*str == '.') {
        --ne;
        ++str;

        while (ne > 0 && isdigit(*str)) {
            --ne;
            number = number * 10. + (*str++ - '0');
            num_digits++;
            num_decimals++;
        }
        exponent -= num_decimals;
    }

    /* test for inf otherwise nan */
    if (num_digits == 0) {
        if ((ne > 2) &&
            ((*str == 'I') ||
             (*str == 'i')) &&
            ((str[1] == 'n') ||
             (str[1] == 'N')) &&
            ((str[2] == 'f') ||
             (str[2] == 'F'))) {
            if (negative == 1)
                *storage = -mxGetInf();
            else
                *storage = mxGetInf();
        } else
            *storage = mxGetNaN();
        while ((ne > 0) && !isanyspace[*str]) { ++str; --ne; }
        return ne;
    }

    /* correct for sign */
    if (negative) number = -number;

    /* process an exponent string */
    if (ne > 0 && (*str == 'e' || *str == 'E')) {
        /* handle optional sign */
        negative = 0;
        --ne;
        if (ne > 0)
            switch(*++str) {
                case '-': negative = 1; /* fall through to increment pos */
                case '+': --ne; ++str;
            }

        /* process string of digits */
        n = 0;
        while (ne > 0 && isdigit(*str)) {
            --ne;
            n = n * 10 + (*str++ - '0');
        }

        if (negative)
            exponent -= n;
        else
            exponent += n;
    }

    /* range check for exponent, FIX: 0e400 is 0! */
    if (exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP) {
        if (number > 0)
            *storage = mxGetInf();
        else if (number < 0)
            *storage = -mxGetInf();
        else
            *storage = 0.0;
        return ne;
    }

    /* scale the result */
    p10 = 10.;
    n = exponent;
    if (n < 0) n = -n;
    while (n) {
        if (n & 1) {
            if (exponent < 0)
                number /= p10;
            else
                number *= p10;
        }
        n >>= 1;
        p10 *= p10;
    }

    /* store result */
    *storage = number;

    /* cut the crap */
    while (ne > 0 && !isanyspace[*str]) { --ne; ++str; }

    /* but return number of remaining chars */
    return ne;
}
int strtodmodc(const mxChar *sstr, int ne, double *storage)
{
    double number;
    int exponent;
    int negative;
    double p10;
    int n;
    int num_digits;
    int num_decimals;

    /* handle optional sign */
    negative = 0;
    switch (SSTR) {
        case '-': negative = 1; /* fall through to increment position */
        case '+': --ne; ++sstr;
    }

    number = 0.;
    num_digits = 0;

    /* process string of digits */
    while (ne > 0 && isdigit((unsigned char) SSTR)) {
        number = number * 10. + (*sstr++ - '0');
        --ne;
        ++num_digits;
    }

    exponent = 0;
    num_decimals = 0;

    /* already reached the end? */
    if (ne == 0) {
        if (negative)
            *storage = -number;
        else
            *storage = number;
        return ne;
    }

    /* process decimal part */
    if (*sstr == '.') {
        --ne;
        ++sstr;

        while (ne > 0 && isdigit((unsigned char) SSTR)) {
            --ne;
            number = number * 10. + (*sstr++ - '0');
            num_digits++;
            num_decimals++;
        }
        exponent -= num_decimals;
    }

    /* test for inf otherwise nan */
    if (num_digits == 0) {
        if ((ne > 2) &&
            ((*sstr == 'I') ||
             (*sstr == 'i')) &&
            ((sstr[1] == 'n') ||
             (sstr[1] == 'N')) &&
            ((sstr[2] == 'f') ||
             (sstr[2] == 'F'))) {
            if (negative == 1)
                *storage = -mxGetInf();
            else
                *storage = mxGetInf();
        } else
            *storage = mxGetNaN();
        while ((ne > 0) && !isanyspace[SSTR]) { ++sstr; --ne; }
        return ne;
    }

    /* correct for sign */
    if (negative) number = -number;

    /* process an exponent string */
    if (ne > 0 && (*sstr == 'e' || *sstr == 'E')) {
        /* handle optional sign */
        negative = 0;
        --ne;
        if (ne > 0)
            switch(*++sstr) {
                case '-': negative = 1; /* fall through to increment pos */
                case '+': --ne; ++sstr;
            }

        /* process string of digits */
        n = 0;
        while (ne > 0 && isdigit((unsigned char) SSTR)) {
            --ne;
            n = n * 10 + (*sstr++ - '0');
        }

        if (negative)
            exponent -= n;
        else
            exponent += n;
    }

    /* range check for exponent, FIX: 0e400 is 0! */
    if (exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP) {
        if (number > 0)
            *storage = mxGetInf();
        else if (number < 0)
            *storage = -mxGetInf();
        else
            *storage = 0.0;
        return ne;
    }

    /* scale the result */
    p10 = 10.;
    n = exponent;
    if (n < 0) n = -n;
    while (n) {
        if (n & 1) {
            if (exponent < 0)
                number /= p10;
            else
                number *= p10;
        }
        n >>= 1;
        p10 *= p10;
    }

    /* store result */
    *storage = number;

    /* cut the crap */
    while (ne > 0 && !isanyspace[SSTR]) { --ne; ++sstr; }

    /* but return number of remaining chars */
    return ne;
}


/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    const int *idim;
    const unsigned char *str = NULL;
    const mxChar *sstr = NULL;
    int ne = 0, nen = 0, nc = 0, nr = 1, col = 0, row = 0, skip = 0;
    double *table;
    double **column;
    double pinfv, ninfv, nanv;

    /* test input arguments */
    if (nrhs < 1 || nlhs > 3)
        mexErrMsgTxt("u8str2double requires at least one input and up to three ouputs.");
    if ((mxGetClassID(*prhs) != mxUINT8_CLASS) &&
        (mxGetClassID(*prhs) != mxCHAR_CLASS))
        mexErrMsgTxt("the input to u8str2double must be of type uint8/char.");
    if (mxGetNumberOfDimensions(*prhs) != 2)
        mexErrMsgTxt("the input to u8str2double must be 2D.");
    idim = mxGetDimensions(*prhs);
    if (*idim++ != 1)
        mexErrMsgTxt("the input to u8str2double must be a row vector.");
    if ((nrhs > 2) &&
        (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[1]) == 1) &&
        (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {
        pinfv = *((const double*) mxGetData(prhs[1]));
        ninfv = *((const double*) mxGetData(prhs[2]));
        if ((!mxIsInf(pinfv)) &&
            (!mxIsNaN(pinfv)) &&
            (!mxIsInf(ninfv)) &&
            (!mxIsNaN(ninfv)) &&
            (pinfv >= 1) &&
            (ninfv >= 1)) {
            nr = (int) pinfv;
            nc = (int) ninfv;
        }
    }
    if ((nrhs > 3) &&
        (mxGetClassID(prhs[3]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
        pinfv = *((const double*) mxGetData(prhs[3]));
        if ((!mxIsInf(pinfv)) &&
            (!mxIsNaN(pinfv)) &&
            (pinfv >= 0))
            skip = (int) pinfv;
    }

    /* we need to get dims ! */
    if (nc == 0) {

        /* get dims */
        ne = *idim;

        /* for uint8 */
        if (mxGetClassID(*prhs) == mxUINT8_CLASS) {
            str = (const unsigned char*) mxGetData(*prhs);

            /* skip parameter */
            if (skip > 0) {
                str = &str[skip];
                ne -= skip;
            }

            /* skip initial spaces */
            while (ne > 0 && isspacev[*str]) { --ne; ++str; }

            /* wait until end of line */
            while (ne > 0 && !islinebreak[*str]) {

                /* number ? */
                if (isfirstdigit[*str]) {

                    /* one more column */
                    ++nc;

                    /* while number */
                    while (isanydigit[*str]) { --ne; ++str; }
                } else if (!isanyspace[*str])
                    mexErrMsgTxt("invalid character in input.");

                /* skip any additional spaces */
                while (ne > 0 && isspacev[*str]) { --ne; ++str; }
            }
        } else {
            sstr = (const mxChar*) mxGetData(*prhs);

            /* skip parameter */
            if (skip > 0) {
                sstr = &sstr[skip];
                ne -= skip;
            }

            /* skip initial spaces */
            while (ne > 0 && isspacev[SSTR]) { --ne; ++sstr; }

            /* wait until end of line */
            while (ne > 0 && !islinebreak[SSTR]) {

                /* number ? */
                if (isfirstdigit[SSTR]) {

                    /* one more column */
                    ++nc;

                    /* while number */
                    while (isanydigit[SSTR]) { --ne; ++sstr; }
                } else if (!isanyspace[SSTR])
                    mexErrMsgTxt("invalid character in input.");

                /* skip any additional spaces */
                while (ne > 0 && isspacev[SSTR]) { --ne; ++sstr; }
            }
        }

        /* check number of columns */
        if (nc == 0)
            mexErrMsgTxt("the input must not have empty rows.");

        /* detect number of rows */
        --ne;
        if (mxGetClassID(*prhs) == mxUINT8_CLASS) {
            ++str;
            while (ne > 0 && islinebreak[*str]) { --ne; ++str; }
            while (ne > 0) {

                /* any content certain! */
                ++nr;
                --ne;
                ++str;

                /* skip all non-linebreaks */
                while(ne > 0 && !islinebreak[*str]) { --ne; ++str; }

                /* skip all line-breaks */
                while(ne > 0 && islinebreak[*str]) { --ne; ++str; }
            }
        } else {
            ++sstr;
            while (ne > 0 && islinebreak[SSTR]) { --ne; ++sstr; }
            while (ne > 0) {

                /* any content certain! */
                ++nr;
                --ne;
                ++sstr;

                /* skip all non-linebreaks */
                while(ne > 0 && !islinebreak[SSTR]) { --ne; ++sstr; }

                /* skip all line-breaks */
                while(ne > 0 && islinebreak[SSTR]) { --ne; ++sstr; }
            }
        }
    }

    /* get special values */
    pinfv = mxGetInf();
    ninfv = -pinfv;
    nanv = mxGetNaN();

    /* create output matrix */
    *plhs = mxCreateDoubleMatrix(nr, nc, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("error creating output matrix.");
    table = (double *) mxGetData(*plhs);
    if (table == NULL)
        mexErrMsgTxt("error getting pointer to output matrix.");

    /* create column pointers */
    column = (double **) mxCalloc(nc, sizeof(double*));
    if (column == NULL)
        mexErrMsgTxt("error allocating column pointers.");
    *column = table;
    for (col = 1; col < nc; ++col)
        column[col] = &(column[col-1][nr]);

    /* return additional outputs: remaining chars, number of rows parsed */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt("error creating 2nd output argument.");
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        if (plhs[2] == NULL)
            mexErrMsgTxt("error creating 3rd output argument.");
    }

    /* start from the beginning */
    ne = *idim;

    /* until the end is reached */
    if (mxGetClassID(*prhs) == mxUINT8_CLASS) {
        str = (const unsigned char*) mxGetData(*prhs);
        if (skip > 0) {
            str = &str[skip];
            ne -= skip;
        }

        while (ne > 0 && nr > row) {

            /* skip spaces */
            while (ne > 0 && isspacev[*str]) { --ne; ++str; }

            /* set column to 0 */
            col = 0;

            /* until row full or end of line */
            while (ne > 0 && !islinebreak[*str] && col < nc) {

                /* conversion with adapted strtodmod */
                nen = strtodmod(str, ne, column[col]);
                str = &str[ne - nen];
                ne = nen;

                /* increase storage counters */
                ++(column[col++]);

                /* skip any further white space before the end of the row */
                while (ne > 0 && isspacev[*str]) { --ne; ++str; }
            }

            /* column counter underflow ? */
            while (col < nc) {
                *column[col] = nanv;
                ++(column[col++]);
            }

            /* skip rest of line if any (too many columns!) */
            while (ne > 0 && !islinebreak[*str]) { --ne; ++str; }

            /* increase row counter */
            ++row;

            /* skip line breaks */
            while (ne > 0 && islinebreak[*str]) { --ne; ++str; }
        }
    } else {
        sstr = (const mxChar*) mxGetData(*prhs);
        if (skip > 0) {
            sstr = &sstr[skip];
            ne -= skip;
        }

        while (ne > 0 && nr > row) {

            /* skip spaces */
            while (ne > 0 && isspacev[SSTR]) { --ne; ++sstr; }

            /* set column to 0 */
            col = 0;

            /* until row full or end of line */
            while (ne > 0 && !islinebreak[SSTR] && col < nc) {

                /* conversion with adapted strtodmod */
                nen = strtodmodc(sstr, ne, column[col]);
                sstr = &sstr[ne - nen];
                ne = nen;

                /* increase storage counters */
                ++(column[col++]);

                /* skip any further white space before the end of the row */
                while (ne > 0 && isspacev[SSTR]) { --ne; ++sstr; }
            }

            /* column counter underflow ? */
            while (col < nc) {
                *column[col] = nanv;
                ++(column[col++]);
            }

            /* skip rest of line if any (too many columns!) */
            while (ne > 0 && !islinebreak[SSTR]) { --ne; ++sstr; }

            /* increase row counter */
            ++row;

            /* skip line breaks */
            while (ne > 0 && islinebreak[SSTR]) { --ne; ++sstr; }
        }
    }

    /* number of rows parsed */
    if (nlhs > 2) {
        table = (double *) mxGetData(plhs[2]);
        if (table != NULL)
            *table = (double) row;
    }

    /* additional data expected ? */
    while (row < nr) {
        ++row;
        for (col = 0; col < nc; ++col) {
            *column[col] = nanv;
            ++(column[col]);
        }
    }

    /* free column pointer */
    mxFree(column);

    /* return additional output: number of chars left */
    if (nlhs > 1) {
        table = (double *) mxGetData(plhs[1]);
        if (table != NULL)
            *table = (double) ne;
    }
}
