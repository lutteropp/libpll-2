/*
    Copyright (C) 2015 Tomas Flouri, with changes and additions by Sarah Lutteropp in 2019

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
%{
#include "parse_rnetwork.h"
#include "pll.h"

#define YY_USER_ACTION preaction();

static size_t string_length = 0;

int pll_rnetwork_colstart = 0;
int pll_rnetwork_colend = 0;

static void preaction()
{
  pll_rnetwork_colstart = pll_rnetwork_colend;
  pll_rnetwork_colend = pll_rnetwork_colend + pll_rnetwork_leng - 1;
}

static char * xstrndup(const char * s, size_t len)
{
  char * p = (char *)malloc(len+1);
  if (!p)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }
  strncpy(p,s,len);
  p[len] = 0;
  return p;
}

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)malloc((*dstlen + srclen + 1)*sizeof(char));
  if (!mem)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }
  memcpy(mem,pll_rnetwork_lval.s,*dstlen);
  strncpy(mem+(*dstlen),src,srclen);
  mem[*dstlen+srclen] = 0;
  if (*dstlen)
    free(pll_rnetwork_lval.s);
  pll_rnetwork_lval.s =  mem;
  *dstlen += srclen;
  return pll_rnetwork_lval.s;
}

%}
%option noyywrap
%option prefix="pll_rnetwork_"
%option nounput
%option noinput
%option yylineno
%x apos
%x quot

%%
<quot>{
\\\"              { append(&string_length, "\\\"", 2);                         }
\'                { append(&string_length, "\'", 1);                           }
\"                { BEGIN(INITIAL); return STRING;                             }
}

<apos>{
\\\'              { append(&string_length, "\\\'", 2);                         }
\"                { append(&string_length, "\"", 1);                           }
\'                { BEGIN(INITIAL); return STRING;                             }
}

<apos,quot>{
\\n               { append(&string_length, "\\n", 2);                          }
\\t               { append(&string_length, "\\t", 2);                          }
\\                { append(&string_length, "\\", 1);                           }
\\\\              { append(&string_length, "\\\\", 2);                         }
([^\"\'\\]|\n)+   { append(&string_length, pll_rnetwork_text, pll_rnetwork_leng);    }
}

\:                { return ':';                                              }
\;                { return ';';                                          }
\#                { return '#';                                               }
\)                { return ')';                                               }
\(                { return '(';                                               }
\,                { return ',';                                              }
\"                { string_length = 0; BEGIN(quot);                            }
\'                { string_length = 0; BEGIN(apos);                            }
\n                { pll_rnetwork_colstart = pll_rnetwork_colend = 0;                 }
[-+]?[0-9]+       { pll_rnetwork_lval.d = xstrndup(pll_rnetwork_text,pll_rnetwork_leng);
                    return NUMBER;                                             }
[+-]?(([0-9]+[\.]?[0-9]*)|([0-9]*[\.]?[0-9]+))([eE][+-]?[0-9]+)? {
                    pll_rnetwork_lval.d = xstrndup(pll_rnetwork_text,pll_rnetwork_leng); 
                    return NUMBER;                                             }
[^ \'\",\(\):;#\[\]\t\n\r][^ \t\n\r\)\(\[\]\,:;#]* {
                    pll_rnetwork_lval.s = xstrndup(pll_rnetwork_text,pll_rnetwork_leng);
                    return STRING;                                             }
[ \t\r]           { ;                                                          }
.                 { snprintf(pll_errmsg, 200, 
                             "Syntax error (%c)\n", pll_rnetwork_text[0]);
                    pll_errno = PLL_ERROR_NEWICK_SYNTAX;
                    return PLL_FAILURE;                                        }
%%
