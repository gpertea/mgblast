/*   trna2tbl.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  trna2tbl.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   5/31/05
*
* $Revision: 1.1 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

/*
*  To compile on Linux:
*
*   gcc -o trna2tbl trna2tbl.c -lm
*
*  To compile on Solaris:
*
*   cc -xildoff -o trna2tbl trna2tbl.c -lgen -lm
*
*  To compile on SGI:
*
*   cc -mips1 -o trna2tbl trna2tbl.c -lm -lPW -lsun
*
*  To compile on Darwin:
*
*   cc -pipe -o trna2tbl trna2tbl.c -lc
*
*/

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* convenient defines and typedefs from NCBI toolkit */

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifndef Pointer
typedef void * Pointer;
#endif

#ifndef Char
typedef char Char, * CharPtr;
#endif

#ifndef Bool
typedef unsigned char Bool, * BoolPtr;
#endif

#ifndef TRUE
#define TRUE ((Bool)1)
#endif

#ifndef FALSE
#define FALSE ((Bool)0)
#endif

#ifndef Int2
typedef short Int2, * Int2Ptr;
#endif

#ifndef Int4
typedef long Int4, * Int4Ptr;
#endif

#ifndef MIN
#define MIN(a,b)	((a)>(b)?(b):(a))
#endif

/* useful portable character macros from NCBI toolkit (assumes ASCII) */

#define IS_DIGIT(c)	('0'<=(c) && (c)<='9')
#define IS_UPPER(c)	('A'<=(c) && (c)<='Z')
#define IS_LOWER(c)	('a'<=(c) && (c)<='z')
#define IS_ALPHA(c)	(IS_UPPER(c) || IS_LOWER(c))
#define TO_LOWER(c)	((Char)(IS_UPPER(c) ? (c)+' ' : (c)))
#define TO_UPPER(c)	((Char)(IS_LOWER(c) ? (c)-' ' : (c)))
#define IS_WHITESP(c) (((c) == ' ') || ((c) == '\n') || ((c) == '\r') || ((c) == '\t'))
#define IS_ALPHANUM(c) (IS_ALPHA(c) || IS_DIGIT(c))
#define IS_PRINT(c)	(' '<=(c) && (c)<='~')

#define MAX_FIELDS  9

static void RunTrnaScan (void)

{
  CharPtr   aa;
  CharPtr   beg;
  Char      buf [256];
  Char      cmmd [256];
  CharPtr   end;
  CharPtr   field [MAX_FIELDS];
  CharPtr   id;
  Int2      idNotSent = TRUE;
  Int2      inBody = FALSE;
  CharPtr   intronBeg;
  CharPtr   intronEnd;
  long int  intronStart;
  long int  intronStop;
  Int2      numFields = 0;
  CharPtr   ptr;
  CharPtr   speed;
  long int  start;
  long int  stop;
  Char      str [80];

/* line by line processing of tRNAscan-SE output table */

  while (fgets (buf, sizeof (buf), stdin) != NULL) {

    if (inBody) {
      memset (field, 0, sizeof (field));

/*
*  parse tab-delimited output line into array of fields, avoiding use of
*  strtok so that empty columns (adjacent tabs) are properly assigned to
*  field array
*/

      ptr = buf;
      for (numFields = 0; numFields < MAX_FIELDS && ptr != NULL; numFields++) {
        field [numFields] = ptr;
        ptr = strchr (ptr, '\t');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
        }
      }

/* interested in ID, start, stop, amino acid, and intron start and stop */

      id = field [0];
      beg = field [2];
      end = field [3];
      aa = field [4];
      intronBeg = field [6];
      intronEnd = field [7];

      if (numFields > 7 &&
          sscanf (beg, "%ld", &start) == 1 &&
          sscanf (end, "%ld", &stop) == 1 &&
          sscanf (intronBeg, "%ld", &intronStart) == 1 &&
          sscanf (intronEnd, "%ld", &intronStop) == 1) {

/* first line of output gives SeqId from FASTA definition line */

        if (idNotSent) {
          fprintf (stdout, ">Features %s\n", id);
          fflush (stdout);
          idNotSent = FALSE;
        }

/* first line of feature has start (tab) stop (tab) feature key */
/* multiple intervals would have lines of start (tab) stop */

        if (intronStart == 0 && intronStop == 0) {
          fprintf (stdout, "%ld\t%ld\ttRNA\n", (long) start, (long) stop);
          fflush (stdout);
        } else {
          fprintf (stdout, "%ld\t%ld\ttRNA\n", (long) start, (long) (intronStart - 1));
          fprintf (stdout, "%ld\t%ld\t\n", (long) (intronStop + 1), (long) stop);
          fflush (stdout);
        }

/* qualifier lines are (tab) (tab) (tab) qualifier key (tab) value */

        if (strstr (aa, "Pseudo") != NULL) {
          fprintf (stdout, "\t\t\tnote\ttRNA-Pseudo\n");
          fflush (stdout);
        } else {
          fprintf (stdout, "\t\t\tproduct\t%s\n", aa);
          fflush (stdout);
        }

/* dash (formerly empty) gene qualifier to suppress /gene (e.g., if tRNA is in an intron) */

        fprintf (stdout, "\t\t\tgene\t-\n");
        fflush (stdout);
      }
    }

/* detect last line of table header, ignoring everything before data section */

    if (strstr (buf, "-----") != NULL) {
      inBody = TRUE;
    }
  }

  if (idNotSent) {
    fprintf (stdout, ">Message\ntRNAscan-SE found no tRNA genes in this sequence\n");
    fflush (stdout);
  }
}

main (int argc, char *argv[])

{
  RunTrnaScan ();
  return 0;
}

