/*   asn2gnb1.c
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
* File Name:  asn2gnb1.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* $Id: asn2gnb1.c,v 1.100 2006/04/13 20:12:51 kans Exp $
*
* Version Creation Date:   10/21/98
*
* $Revision: 1.100 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objpubme.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbi.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

/* utility functions */

NLM_EXTERN ValNodePtr ValNodeCopyStrToHead (ValNodePtr PNTR head, Int2 choice, CharPtr str)

{
  ValNodePtr newnode;

  if (head == NULL || str == NULL) return NULL;

  newnode = ValNodeNew (NULL);
  if (newnode == NULL) return NULL;

  newnode->choice = (Uint1) choice;
  newnode->data.ptrvalue = StringSave (str);

  newnode->next = *head;
  *head = newnode;

  return newnode;
}

/* the val node strings mechanism will be replaced by a more efficient method later  */

NLM_EXTERN CharPtr MergeFFValNodeStrs (
  ValNodePtr list
)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  if (list == NULL) return NULL;

  for (vnp = list, len = 0; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    len += StringLen (str);
  }
  if (len == 0) return NULL;

  ptr = MemNew (sizeof (Char) * (len + 2));
  if (ptr == NULL) return NULL;

  for (vnp = list, tmp = ptr; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    tmp = StringMove (tmp, str);
  }

  return ptr;
}


NLM_EXTERN void AddValNodeString (
  ValNodePtr PNTR head,
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix
)

{
  Char     buf [256];
  CharPtr  freeme = NULL;
  size_t   len;
  CharPtr  newstr;
  CharPtr  strptr;

  len = StringLen (prefix) + StringLen (string) + StringLen (suffix);
  if (len == 0) return;

  if (len < sizeof (buf)) {

    /* if new string fits in stack buffer, no need to allocate */

    MemSet ((Pointer) buf, 0, sizeof (buf));
    newstr = buf;

  } else {

    /* new string bigger than stack buffer, so allocate sufficient string */

    newstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
    if (newstr == NULL) return;

    /* allocated string will be freed at end of function */

    freeme = newstr;
  }

  strptr = newstr;

  if (prefix != NULL) {
    strptr = StringMove (strptr, prefix);
  }

  if (string != NULL) {
    strptr = StringMove (strptr, string);
  }

  if (suffix != NULL) {
    strptr = StringMove (strptr, suffix);
  }

  /* currently just makes a valnode list, to be enhanced later */

  ValNodeCopyStr (head, 0, newstr);

  /* if large string was allocated, free it now */

  if (freeme != NULL) {
    MemFree (freeme);
  }
}


NLM_EXTERN void FFAddString_NoRedund (
  StringItemPtr unique,
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix
)
{
  CharPtr    str = string;
  Int4       foundPos = 0;
  Boolean    wholeWord = FALSE;

  if ( StringHasNoText(prefix)  &&
       StringHasNoText(string)  &&
       StringHasNoText(suffix)  ) return;

  if (StringNICmp (string, "tRNA-", 5) == 0) {
    str = string+5;
    }

  while ( foundPos >= 0 && !wholeWord ) {
    foundPos = FFStringSearch(unique, str, foundPos);
    if ( foundPos >= 0 ) {
      wholeWord = IsWholeWordSubstr(unique, foundPos, str);
      foundPos += StringLen(str);
    }
  }

  if ( foundPos < 0 || !wholeWord ) {
      FFAddTextToString(unique, prefix, string, suffix, FALSE, FALSE, TILDE_IGNORE);
  }
}



/* s_AddPeriodToEnd () -- Adds a '.' to the end of a given string if */
/*                        there is not already one there.            */
/*                                                                   */
/*                        Note that this adds one character to the   */
/*                        length of the string, leading to a         */
/*                        memory overrun if space was not previously */
/*                        allocated for this.                        */

NLM_EXTERN void s_AddPeriodToEnd (CharPtr someString)
{
  Int4  len;

  if (StringHasNoText (someString)) return;
  len = StringLen (someString);
  if (len < 1) return;
  if (someString[len-1] != '.')
    {
      someString[len] = '.';
      someString[len+1] = '\0';
    }
}

/* s_RemovePeriodFromEnd () -- If the last character in a given      */
/*                             string is a '.', removes it.          */

NLM_EXTERN Boolean s_RemovePeriodFromEnd (CharPtr someString)
{
  Int4  len;

  if (StringHasNoText (someString)) return FALSE;
  len = StringLen (someString);
  if (len < 1) return FALSE;
  if (someString[len-1] == '.') {
    someString[len-1] = '\0';
    return TRUE;
  }
  return FALSE;
}

/**/
/*   isEllipsis () - Determines if a string ends in an ellipses */
/**/

NLM_EXTERN Boolean IsEllipsis (
  CharPtr str
)

{
  size_t   len;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;
  len = StringLen (str);
  if (len < 3) return FALSE;
  ptr = str + len - 3;
  return (Boolean) (ptr [0] == '.' && ptr [1] == '.' && ptr [2] == '.');
}

NLM_EXTERN void A2GBSeqLocReplaceID (
  SeqLocPtr newloc,
  SeqLocPtr ajpslp
)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (ajpslp);
  if (bsp == NULL) return;
  sip = SeqIdFindBest (bsp->id, 0);
  SeqLocReplaceID (newloc, sip);
}

NLM_EXTERN CharPtr asn2gb_PrintDate (
  DatePtr dp
)

{
  Char    buf [30];
  size_t  len;

  if (dp == NULL) return NULL;

  if (DatePrint (dp, buf)) {
    if (StringICmp (buf, "Not given") != 0) {
      len = StringLen (buf);
      if (len > 0) {
        if (buf [len - 1] == '\n') {
          if (buf [len - 2] == '.') {
            buf [len - 2] = '\0';
          } else {
            buf [len - 1] = '\0';
          }
        }
      }
      return StringSave (buf);
    }
  }

  return NULL;
}

static CharPtr month_names [] = {
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
  "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
  "??"
};

NLM_EXTERN CharPtr DateToFF (
  CharPtr buf,
  DatePtr dp,
  Boolean citSub
)

{
  Int2  day;
  Int2  month;
  Int2  year;

  if (buf != NULL) {
    *buf = '\0';
  }
  if (dp == NULL) return NULL;

  if (dp->data [0] == 0) {

    StringCpy (buf, dp->str);

  } else if (dp->data [0] == 1) {

    year = 1900 + (Int2) dp->data [1];
    month = (Int2) dp->data [2];
    day = (Int2) dp->data [3];

    if (citSub) {
      if (month < 1 || month > 12) {
        month = 13;
      }
      if (day < 1 || day > 31) {
        day = 0;
      }
    } else {
      if (month < 1 || month > 12) {
        month = 1;
      }
      if (day < 1 || day > 31) {
        day = 1;
      }
    }

    if (day < 1) {
      sprintf (buf, "\?\?-%s-%ld",
               month_names [month-1], (long) year);
    } else if (day < 10) {
      sprintf (buf, "0%ld-%s-%ld",
               (long) day, month_names [month-1], (long) year);
    } else {
      sprintf(buf, "%ld-%s-%ld",
               (long) day, month_names [month-1], (long) year);
    }
  }

  return buf;
}


NLM_EXTERN StringItemPtr FFGetString (IntAsn2gbJobPtr ajp)

{
  StringItemPtr  sip;

  if (ajp == NULL) return NULL;
  if (ajp->pool != NULL) {
    sip = ajp->pool;
    ajp->pool = sip->next;
    sip->next = NULL;
    MemSet ((Pointer) sip, 0, sizeof (StringItem));
  } else {
    sip = (StringItemPtr) MemNew (sizeof (StringItem));
    if (sip == NULL) return NULL;
  }
  sip->curr = sip;
  sip->iajp = ajp;
  sip->pos = 0;
  return sip;
}

NLM_EXTERN void FFRecycleString (IntAsn2gbJobPtr ajp, StringItemPtr ffstring)

{
  StringItemPtr  nxt;

  if (ajp == NULL || ffstring == NULL) return;
  if ( ffstring->pos == -1 ) return;
  
  nxt = ffstring;
  nxt->pos = -1;
  while (nxt->next != NULL) {
    nxt->pos = -1;
    nxt = nxt->next;
  }
  nxt->next = ajp->pool;
  ajp->pool = ffstring;

  ffstring->curr = NULL;
}

NLM_EXTERN void FFAddOneChar (
  StringItemPtr sip, 
  Char ch,
  Boolean convertQuotes
)
{
  StringItemPtr current = sip->curr;

  if ( current->pos == STRING_BUF_LEN ) {
    current->next = FFGetString(sip->iajp);
    current = current->next;
    current->pos = 0;
    sip->curr = current;
  }

  if ( convertQuotes && ch == '\"' ) {
    ch = '\'';
  }
  current->buf[current->pos] = ch;
  current->pos++;
}

NLM_EXTERN void FFAddNewLine(StringItemPtr ffstring) {
  FFAddOneChar(ffstring, '\n', FALSE);
}

NLM_EXTERN void FFAddNChar (
  StringItemPtr sip, 
  Char ch,
  Int4 n,
  Boolean convertQuotes
)
{
  Int4 i;

  for ( i = 0; i < n; ++i ) {
    FFAddOneChar(sip, ch, convertQuotes);
  }
}
  

NLM_EXTERN void FFExpandTildes (StringItemPtr sip, CharPtr PNTR cpp) {
  Char replace = **cpp;

  if ( **cpp == '~' ) {
    if ( *((*cpp) + 1) == '~' ) {     /* "~~" -> '~' */
      replace = '~';
      (*cpp)++;
    } else {
      replace = '\n';
    }
  } 

  FFAddOneChar(sip, replace, FALSE);
}


NLM_EXTERN void FFReplaceTildesWithSpaces (StringItemPtr ffstring, CharPtr PNTR cpp) {
  Char replace = **cpp, lookahead;
  CharPtr cptr = *cpp;
  
  if ( *cptr == '`' ) {
    FFAddOneChar(ffstring, replace, FALSE);
    return;
  }

  replace = ' ';
  lookahead = *(cptr + 1);

  if ( IS_DIGIT(lookahead) ) {
    replace = '~';
  }
  else {
    if ( (lookahead == ' ') || (lookahead == '(') ) {
      if ( IS_DIGIT(*(cptr + 2)) ) {
        replace = '~';
      }
    }
  }

  FFAddOneChar(ffstring, replace, FALSE);
}

NLM_EXTERN void FFOldExpand (StringItemPtr sip, CharPtr PNTR cpp) {
  /* "~" -> "\n", "~~" or "~~ ~~" -> "\n\n" */ 
  CharPtr cp = *cpp;
  Char current = *cp;
  Char next = *(cp + 1);
  
  /* handle "'~" */
  if ( current == '`' ) {
    if ( next != '~' ) {
        FFAddOneChar(sip, current, FALSE);
    } else {
        FFAddOneChar(sip, '~', FALSE);
        (*cpp)++;
    }
    return;
  }

  /* handle "~", "~~" or "~~ ~~" */
  FFAddOneChar(sip, '\n', FALSE);
  if ( next == '~' ) {
    FFAddOneChar(sip, '\n', FALSE);
    cp++;
    *cpp = cp;
    cp++;
    if ( *cp == ' ' ) {
      cp++;
      if ( *cp == '~' ) {
        cp++;
        if ( *cp == '~' ) { /* saw "~~ ~~" */
          *cpp = cp;
        }
      }
    }
  }
}

NLM_EXTERN void AddCommentStringWithTildes (StringItemPtr ffstring, CharPtr string)
{
/* One "~" is a  new line, "~~" or "~~ ~~" means 2 returns */       

    /* Int2  i; */

    while (*string != '\0') {
        if (*string == '`' && *(string+1) == '~') {
            FFAddOneChar(ffstring, '~', FALSE);
            string += 2;
        } else if (*string == '~') {
            FFAddOneChar(ffstring, '\n', FALSE);
            string++;
            if (*string == '~') {
                /*
                for (i = 0; i < 12; i++) {
                    FFAddOneChar(ffstring, ' ', FALSE);
                }
                */
                FFAddOneChar(ffstring, '\n', FALSE);
                string++;
            if (*string == ' ' && *(string+1) == '~' && *(string+2) == '~') {
                    string += 3;
            }
          }
        } else if (*string == '\"') {   
            *string = '\'';
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        } else {  
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        }
    }
}    /* AddCommentStringWithTildes */


NLM_EXTERN void AddStringWithTildes (StringItemPtr ffstring, CharPtr string)
{
/* One "~" is a  new line, "~~" or "~~ ~~" means 2 returns */       

    while (*string != '\0') {
        if (*string == '`' && *(string+1) == '~') {
            FFAddOneChar(ffstring, '~', FALSE);
            string += 2;
        } else if (*string == '~') {
            FFAddOneChar(ffstring, '\n', FALSE);
            string++;
            if (*string == '~') {
                FFAddOneChar(ffstring, '\n', FALSE);
                string++;
        if (*string == ' ' && *(string+1) == '~' && *(string+2) == '~') {
                    string += 3;
        }
      }
        } else if (*string == '\"') {   
            *string = '\'';
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        } else {  
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        }
    }
}    /* AddStringWithTildes */


NLM_EXTERN void FFProcessTildes (StringItemPtr sip, CharPtr PNTR cpp, Int2 tildeAction) {
    
  switch (tildeAction) {

  case TILDE_EXPAND :
      FFExpandTildes(sip, cpp);
      break;

  case TILDE_OLD_EXPAND :
      FFOldExpand(sip, cpp);
      break;

  case TILDE_TO_SPACES :
      FFReplaceTildesWithSpaces (sip, cpp);
      break;

  case TILDE_IGNORE:
  default:
      FFAddOneChar(sip, **cpp, FALSE);
      break;
  }
}

NLM_EXTERN void FFAddPeriod (StringItemPtr sip) {
  Int4 i;
  Char ch  = '\0';
  StringItemPtr riter = sip->curr, prev;
  IntAsn2gbJobPtr ajp;

  if ( sip == NULL ) return;
  ajp = (IntAsn2gbJobPtr)sip->iajp;
  if ( ajp == NULL ) return;

  for ( i = riter->pos - 1; i >= 0; --i ) {
    ch = riter->buf[i];

    if ( (ch == ' ') || (ch == '\t')  || (ch == '~')  || (ch == '.') || (ch == '\n')) {
      riter->pos--;
      
      if ( i < 0 && riter != sip ) {
        for ( prev = sip; prev->next != NULL; prev = prev->next ) {
          if ( prev->next == riter ) {
            i = prev->pos - 1;
            FFRecycleString(ajp, riter);
            riter = prev;
            riter->next = NULL;
              sip->curr = riter;
            break;
          }
        }
      }

    } else {
      break;
    }
  }

  if (ch != '.') {
    FFAddOneChar(sip, '.', FALSE);
  }
}

NLM_EXTERN void FFAddOneString (
  StringItemPtr sip, 
  CharPtr string,
  Boolean addPeriod, 
  Boolean convertQuotes,
  Int2 tildeAction
)
{
  CharPtr strp = string;

  if ( string == NULL ) return;
  
  while ( *strp != '\0' ) {
    if ( (*strp == '`') || (*strp == '~') ) {
      FFProcessTildes(sip, &strp, tildeAction);
    } else {
      FFAddOneChar(sip, *strp, convertQuotes);
    }
    strp++;
  }

  if ( addPeriod ) {
    FFAddPeriod(sip);
  }
}

NLM_EXTERN void FFCatenateSubString (
  StringItemPtr dest,
  StringItemPtr start_sip, Int4 start_pos,
  StringItemPtr end_sip, Int4 end_pos
)
{
  Int4 max_i, min_i, i;
  StringItemPtr current;
  Boolean in_url = FALSE;
  IntAsn2gbJobPtr ajp = (IntAsn2gbJobPtr)dest->iajp;

  if ( GetWWW(ajp) ) {
    for ( current = start_sip, i = start_pos;
    current != NULL; 
    current = current->next ) {
      if ( current == start_sip ) {
        min_i = start_pos;
      } else {
        min_i = 0;
      }
      
      if ( current == end_sip ) {
        max_i = end_pos;
      } else {
        max_i = current->pos;
      }
      
      for ( i = min_i; i < max_i; ++i ) {
        if ( current->buf[i] == '<' ) {
          if ( !FFIsStartOfLink(current, i) ) {
            FFAddOneString(dest, "&lt;", FALSE, FALSE, TILDE_IGNORE);
            continue;
          } else {
            in_url = TRUE;
          }
        }
        if ( current->buf[i] == '>' ) {
          if ( !in_url ) {
            FFAddOneString(dest, "&gt;", FALSE, FALSE, TILDE_IGNORE);
            continue;
          } else {
            in_url = FALSE;
          }
        } 

        FFAddOneChar(dest, current->buf[i], FALSE);
      }

      if ( current == end_sip ) break;
    }
  } else {
    for ( current = start_sip, i = start_pos;
    current != NULL; 
    current = current->next ) {
      if ( current == start_sip ) {
        min_i = start_pos;
      } else {
        min_i = 0;
      }
      
      if ( current == end_sip ) {
        max_i = end_pos;
      } else {
        max_i = current->pos;
      }
      
      for ( i = min_i; i < max_i; ++i ) {
        FFAddOneChar(dest, current->buf[i], FALSE);
      }
      
      if ( current == end_sip ) break;
    }
  }
}


NLM_EXTERN CharPtr FFToCharPtr (StringItemPtr sip) {
  Int4 size = 0;
  StringItemPtr iter;
  CharPtr result, temp;

  for ( iter = sip; iter != NULL; iter = iter->next ) {
    size += iter->pos;
  }

  result = (CharPtr)MemNew(size + 2);
  temp = result;

  for ( iter = sip; iter != NULL; iter = iter->next ) {
    MemCpy( temp, iter->buf, iter->pos );
    temp += iter->pos;
  }

  *temp = '\0';

  return result;
}



/* word wrap functions */

NLM_EXTERN void FFSkipLink (StringItemPtr PNTR iterp, Int4Ptr ip) {
  StringItemPtr iter = *iterp;
  Int4 i = *ip;

  while ( (iter != NULL) && (iter->buf[i] != '>') ) {
    ++i;

    if ( i == iter->pos ) {
      iter = iter->next;
      i = 0;
    }
  }
  ++i;
  if ( i == iter->pos && iter->next != NULL ) {
    iter = iter->next;
    i = 0;
  }

  *iterp = iter;
  *ip = i;
}

NLM_EXTERN Boolean FFIsStartOfLink (StringItemPtr iter, Int4 pos)  {
  static CharPtr start_link = "<A HREF";
  static CharPtr end_link = "</A>";
  Int4 start_len = StringLen(start_link);
  Int4 end_len = StringLen(end_link);
  Char temp[10];
  Int4 i;

  if ( iter == NULL || pos >= iter->pos ) return FALSE;
  if ( iter->buf[pos] != '<' ) return FALSE;

  MemSet(temp, 0, sizeof(temp));
  for ( i = 0; i < start_len && iter != NULL; ++i ) {
    if ( pos + i < iter->pos ) {
      temp[i] = iter->buf[pos+i];
      if ( i == end_len - 1 ) {
        if ( StringNICmp(temp, end_link, end_len) == 0 ) {
          return TRUE;
        }
      }
    } else {
      iter = iter->next;
      pos = -i;
      --i;
    }
  }

  if ( i == start_len ) {
    if ( StringNICmp(temp, start_link, start_len) == 0 ) {
        return TRUE;
    }
  }

  return FALSE;
}


NLM_EXTERN void FFSavePosition(StringItemPtr ffstring, StringItemPtr PNTR bufptr, Int4 PNTR posptr) {
  *bufptr = ffstring->curr;
  *posptr = ffstring->curr->pos;
}


NLM_EXTERN void FFTrim (
    StringItemPtr ffstring,
    StringItemPtr line_start,
    Int4 line_pos,
    Int4 line_prefix_len
)
{
  StringItemPtr riter, iter;
  Int4 i;
  IntAsn2gbJobPtr ajp = (IntAsn2gbJobPtr)ffstring->iajp;

  for ( i = 0; i < line_prefix_len; ++i ) {
    ++line_pos;
    if ( line_pos == STRING_BUF_LEN ) {
      line_pos = 0;
      line_start= line_start->next;
    }
  }

  riter = ffstring->curr;
  while ( riter != NULL ) {
    for ( i = riter->pos - 1;
          (i >= 0) && !(riter == line_start && i <= line_pos);
          --i ) {
      if ( !IS_WHITESP(riter->buf[i]) || (riter->buf[i] == '\n') ) {
        break;
      }
    }
    if ( i < 0 ) {
      i = STRING_BUF_LEN - 1;
      for ( iter = ffstring; iter != NULL; iter = iter->next ) {
        if ( iter->next == riter ) {
          break;
        }
      }
      if ( iter == NULL ){
        ffstring->pos = 0;
        break;
      } else {
        
        riter = iter;
        ffstring->curr = riter;
      }
    } else {
      riter->pos = i + 1;
      FFRecycleString(ajp, riter->next);
      riter->next = NULL;
      break;
    }
  }
}



/* A line is wrapped when the visble text in th eline exceeds the line size. */
/* Visible text is text that is not an HTML hyper-link.                      */
/* A line may be broken in one of the following characters:                  */
/* space, comma and dash                                                     */
/* the oredr of search is first spaces, then commas and then dashes.         */
/* We nee to take into account the possiblity that a 'new-line' character    */
/* already exists in the line, in such case we break at the 'new-line'       */
/* spaces, dashes and new-lines will be broken at that character wheras for  */
/* commas we break at the character following the comma.                     */

NLM_EXTERN void FFCalculateLineBreak (
  StringItemPtr PNTR break_sip, Int4 PNTR break_pos,
  Int4 init_indent, Int4 visible
)
{
  StringItemPtr iter, prev;
  Int4 i,
       done = FALSE,
       copied = 0, 
       start = *break_pos;
  Char ch;
  Boolean found_comma = FALSE, found_dash = FALSE;
  /* each candidate is a pair of buffer and position withingh this buffer */
  StringItemPtr candidate_sip_space = NULL,
                candidate_sip_comma = NULL,
                candidate_sip_dash  = NULL;
  Int4          candidate_int_space = -1,  
                candidate_int_comma = -1,
                candidate_int_dash  = -1;
  

  iter = *break_sip;
  prev = iter;

  /* skip the first 'init_indent' characters of the line */
  while ( iter != NULL && !done ) {
    for ( i = start; i < iter->pos && init_indent > 0; ++i ) {
      if ( iter->buf[i] == '\n' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
        done = TRUE;
        break;
      }
      if ( FFIsStartOfLink(iter, i) ) {
        FFSkipLink(&iter, &i);
        --i;
        continue;
      }

      --init_indent;
      ++copied;
    }
    if ( init_indent > 0 ) {
      start = 0;
      iter = iter->next;
    } else {
      break;
    }
  }
  start = i;

  while ( iter != NULL && !done ) {
    for ( i = start; i < iter->pos; ++i ) {
      if ( found_comma ) {
        candidate_sip_comma = iter;
        candidate_int_comma = i;
        found_comma = FALSE;
      }
      if ( found_dash ) {
        candidate_sip_dash = iter;
        candidate_int_dash = i;
        found_dash= FALSE;
      }

      ch = iter->buf[i];
      if ( ch == '\n' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
        done = TRUE;
        break;
      } else if ( ch == ' ' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
      } else if ( ch == ',' ) {
        found_comma = TRUE;
      } else if ( ch == '-' ) {
        found_dash = TRUE;
        /*candidate_sip_dash = iter;
        candidate_int_dash = i;*/
      }

      if ( FFIsStartOfLink(iter, i) ) {
        FFSkipLink(&iter, &i);
        --i;
        continue;
      }

      ++copied;
      if ( copied >= visible ) {
        if ( (candidate_sip_space == NULL) && (candidate_int_space == -1) &&
             (candidate_sip_comma == NULL) && (candidate_int_comma == -1) &&
             (candidate_sip_dash == NULL)  && (candidate_int_dash == -1)  ) {
          candidate_sip_space = iter;
          candidate_int_space = i;
        }
        done = TRUE;
        break;
      }      
    }
    start = 0;
    if ( !done ) {
      prev = iter;
      iter = iter->next;
    }
  }
  
  /* the order in which we examine the various candidate breaks is important */
  if ( iter == NULL && !done) { /* reached the end */
    *break_sip = prev;
    *break_pos = prev->pos;
  } else {
    if( candidate_sip_space != NULL ) {
        *break_sip = candidate_sip_space;
        *break_pos = candidate_int_space;
    } else if( candidate_sip_comma != NULL ) {
        *break_sip = candidate_sip_comma;
      *break_pos = candidate_int_comma;
    } else if( candidate_sip_dash != NULL ) {
      *break_sip = candidate_sip_dash;
      *break_pos = candidate_int_dash;
    }
  }
}

/*
 * Scans the given buffer froma given scan position, for the next occurrence of 
 * the indicated character. The search breaks when the character is found, or the
 * supplied break position is reached.
 * On exit, the scan position will either be on the character found, or at the 
 * given break position.
 *
 *  *p_line_sip:    in: points to the buffer where scan should start
 *                  out: points to the buffer where the scan ended
 *  *p_line_pos:    in: points to the position in *p_line_sip where the scan should
 *                    start
 *                  out: points to the position in *p_line_sip where the scan ended.
 *  break_sip:      points to buffer where the scan should stop
 *  break_pos:      position in *break_sip where the scan should stop
 *  c:              the character we are looking for
 */
NLM_EXTERN Boolean FFFindSingleChar(
  StringItemPtr* p_line_sip,
  Int4* p_line_pos,
  StringItemPtr break_sip,
  Int4 break_pos,
  char c )
{
  while( *p_line_pos >= (*p_line_sip)->pos) {
    *p_line_pos -= (*p_line_sip)->pos;
    (*p_line_sip) = (*p_line_sip)->next;
    if ( *p_line_sip == NULL ) {
      return FALSE;
    }
  }
  while (*p_line_sip != break_sip){
    while (*p_line_pos <(*p_line_sip)->pos) {
      if ((*p_line_sip)->buf[ *p_line_pos ] == c)
        return TRUE;
      else
        ++(*p_line_pos);
    }
    *p_line_pos = 0;
    *p_line_sip = (*p_line_sip)->next;
  }
  while (*p_line_pos < break_pos){
    if ( (*p_line_sip)->buf[ *p_line_pos ] == c )
      return TRUE;
    else
      ++(*p_line_pos);
  }
  return FALSE;
}

/*
 * Returns the number of bytes remaining in the buffer chain, starting from the 
 * given buffer and a read mark inside it.
 *
 *  sip:            points to the buffer where the string starts,
 *  cur_pos:        read mark in the buffer
 */
NLM_EXTERN Int4 FFRemainingLength(
  StringItemPtr sip,
  Int4 cur_pos )
{
  return FFLength(sip)-cur_pos;
}

/*
 * Scans the given line for the next opening tag of an HTML hyperlink. Ajusts 
 * the line position to immediately after the opening tag (if such a tag is 
 * found) or the the end of the line (if no such tag is found).
 * If a character buffer is supplied, this function will copy any opening tag
 * it finds into that buffer.
 *
 *  *p_line_sip:    in: points to the string buffer where the scan should start
 *                  out: points to the string buffer where the scan ended
 *  *p_line_pos:    in: position in **p_start_sip where the scan should start
 *                  out: position in **p_start_sip where the scan ended
 *  break_sip:      buffer that contain the line break
 *  break_pos:      position in break_sip that represents the line break
 *  buf_open_link:  character buffer to hold a copy of the opening link found
 *                  (or =0 if this information is not required).
 */
NLM_EXTERN Boolean FFExtractNextOpenLink(
  StringItemPtr* p_line_sip,
  Int4* p_line_pos,
  StringItemPtr break_sip,
  Int4 break_pos,
  char* buf_open_link )
{
  int i;

  const char* buf_markup_open = "<A HREF";
  const int markup_size = strlen(buf_markup_open);

   while ((*p_line_sip != break_sip) || (*p_line_pos < break_pos)) {

    if (FFFindSingleChar(p_line_sip, p_line_pos, break_sip, break_pos, '<' )) {
      
      if (FFRemainingLength(*p_line_sip, *p_line_pos) < markup_size) {
        *p_line_sip = break_sip;
        *p_line_pos = break_pos;
        return FALSE;
      }
      for ( i=0; i < markup_size; ++i ) {
        if (buf_markup_open[i] != toupper( FFCharAt( *p_line_sip, (*p_line_pos)+i )))
          break;
      }
      if ( i == markup_size ) {
        if (buf_open_link != 0) {
          char ch = FFCharAt( *p_line_sip, *p_line_pos );
          for (i=0; '>' != (buf_open_link[i] = FFCharAt( *p_line_sip, *p_line_pos )); ++i)
            ++(*p_line_pos);
          buf_open_link[i+1] = 0;
        } else {
          *p_line_pos += markup_size;
        }

        return TRUE;
      } else {
        ++(*p_line_pos);
      }
    }
  }      
  return FALSE;
}

/*
 * Scans the given line for the next closing tag of an HTML hyperlink. Ajusts 
 * the line position to immediately after the closing tag (if such a tag is 
 * found) or the the end of the line (if no such tag is found).
 *
 *  *p_line_sip:    in: points to the string buffer where the scan should start
 *                  out: points to the string buffer where the scan ended
 *  *p_line_pos:    in: position in **p_start_sip where the scan should start
 *                  out: position in **p_start_sip where the scan ended
 *  break_sip:      buffer that contain the line break
 *  break_pos:      position in break_sip that represents the line break
 */
NLM_EXTERN Boolean FFExtractNextCloseLink(
  StringItemPtr* p_line_sip,
  Int4* p_line_pos,
  StringItemPtr break_sip,
  Int4 break_pos )
{
  int i;

  const char* buf_close_link = "</A>";
  const int markup_close_size = strlen(buf_close_link);

   while ((*p_line_sip != break_sip) || (*p_line_pos < break_pos)) {

    if (FFFindSingleChar(p_line_sip, p_line_pos, break_sip, break_pos, '<' )) {
      if (FFRemainingLength(*p_line_sip, *p_line_pos) < markup_close_size) {
        *p_line_sip = break_sip;
        *p_line_pos = break_pos;
        return FALSE;
      }
      for ( i=0; i < markup_close_size; ++i ) {
        if (buf_close_link[i] != toupper(FFCharAt( *p_line_sip, (*p_line_pos)+i)))
          break;
      }
      if (i == markup_close_size) {
          (*p_line_pos) += markup_close_size;
        return TRUE;
      } else {
        ++(*p_line_pos);
      }
    }
  }      
  return FALSE;
}

/*  
 * Checks a given line whether its end falls between the opening and the closing
 * tag of an HTML link.
 *
 *  start_sip:      string buffer where the given line starts,
 *  start_pos:      position in start_sip where the given line starts,
 *  break_sip:      string buffer where the given line ends,
 *  break_pos:      position in break_pos where the given line ends,
 *  buf_link_open:  optional buffer where the open tag of the split link will be
 *                  written to. Leave =0 if you don't need this.
 */
NLM_EXTERN Boolean FFLineBreakSplitsHtmlLink( 
  StringItemPtr start_sip, 
  Int4 start_pos, 
  StringItemPtr break_sip, 
  Int4 break_pos,
  char* buf_link_open ) 
{
  int open_count = 0;

  StringItemPtr cur_iter=0;
  int cur_pos=0;

  if ( ! GetWWW((IntAsn2gbJobPtr)start_sip->iajp) )
      return FALSE;
  
  cur_iter = start_sip;
  cur_pos = start_pos;

  while ((cur_iter != break_sip) || (cur_pos < break_pos)) {
    switch(open_count) {
    case 0:
      if (FFExtractNextOpenLink(&cur_iter, &cur_pos, break_sip, break_pos, buf_link_open ))
        ++open_count;
      break;
    case 1:
      if (FFExtractNextCloseLink(&cur_iter, &cur_pos, break_sip, break_pos ))
        --open_count;
      break;
    default:
      break;
    }
  }
  return (open_count);
} /*FFLineBreakSplitsHtmlLink*/

NLM_EXTERN void FFLineWrap (
  StringItemPtr dest, 
  StringItemPtr src, 
  Int4 init_indent,
  Int4 cont_indent, 
  Int4 line_max,
  CharPtr eb_line_prefix
)
{
  /* line break candidate is a pair <StringItemPtr, position> */
  StringItemPtr break_sip = src;
  Int4          break_pos = 0;
  StringItemPtr line_start = NULL;
  Int4          line_pos = 0;
  Int4          i, line_prefix_len = 0;
  StringItemPtr iter;
  Boolean       cont = FALSE;

  FFSavePosition(dest, &line_start, &line_pos);

  for ( iter = src; iter != NULL; iter = iter->next ) {
    for ( i = 0; i < iter->pos; ) {

      Boolean linebreak_splits_link = FALSE;
      char buf_split_link_open[ 1024 ];
      const char* buf_split_link_close = "</a>";
      MemSet( (void*)buf_split_link_open, 0, sizeof(buf_split_link_open) );

      break_pos = i;
      break_sip = iter;

      FFCalculateLineBreak(
        &break_sip, &break_pos, init_indent, line_max - line_prefix_len + 1);
      linebreak_splits_link = 
        FFLineBreakSplitsHtmlLink(iter, i, break_sip, break_pos, 
          buf_split_link_open);
      FFCatenateSubString(dest, iter, i, break_sip, break_pos);
      FFTrim(dest, line_start, line_pos, cont_indent);
      if ( linebreak_splits_link ) {
        FFAddOneString( dest, 
          (char*)buf_split_link_close, FALSE, FALSE, TILDE_IGNORE );
      }
      FFAddOneChar(dest, '\n', FALSE);
      
      FFSavePosition(dest, &line_start, &line_pos);

      /* for EMBL 'XX' lines */
      if (eb_line_prefix != NULL) {
        cont = FALSE;
        if (break_pos > 1) {
          if (break_sip->buf[break_pos-1] == 'X' && break_sip->buf[break_pos-2] == 'X') {
            if ((break_pos == 2) || (break_sip->buf[break_pos-3] == '\n')) {
              ++break_pos;
              cont = TRUE;
            }
          }
        } else if (break_pos == 1) {
          if (break_sip->buf[0] == 'X' && iter->buf[iter->pos-1] == 'X') {
            if ((iter->pos > 1)  &&  iter->buf[iter->pos-2] == '\n') {
              ++break_pos;
              cont = TRUE;
            }
          }
        }
      }

      i = break_pos;
      iter = break_sip;

      if (cont) continue;

      if ( IS_WHITESP(iter->buf[i]) ) {
        i++;
      }
      if ( iter != src->curr || i < iter->pos ) {
        if ( eb_line_prefix != NULL ) {
          FFAddOneString(dest, eb_line_prefix, FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddNChar(dest, ' ', cont_indent - StringLen(eb_line_prefix), FALSE);
        if ( linebreak_splits_link ) {
          FFAddOneString( dest, buf_split_link_open, FALSE, FALSE, TILDE_IGNORE );
        }
        init_indent = 0;
        line_prefix_len = cont_indent;
        /*FFSkipGarbage(&iter, &i);*/
      }
    }
  }
}

/* === */

NLM_EXTERN void FFStartPrint (
  StringItemPtr sip,
  FmtType format,
  Int4 gb_init_indent,
  Int4 gb_cont_indent,
  CharPtr gb_label,
  Int4 gb_tab_to,
  Int4 eb_init_indent,
  Int4 eb_cont_indent,
  CharPtr eb_line_prefix,
  Boolean eb_print_xx 
)

{
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    FFAddNChar(sip, ' ', gb_init_indent, FALSE);
    FFAddOneString(sip, gb_label, FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(sip, ' ', gb_tab_to - gb_init_indent - StringLen(gb_label), FALSE);
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    if ( eb_print_xx ) {
      FFAddOneString(sip, "XX\n", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString(sip, eb_line_prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(sip, ' ', eb_init_indent - StringLen(eb_line_prefix), FALSE);
  }
}

NLM_EXTERN void FFAddTextToString (
  StringItemPtr ffstring, 
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix,
  Boolean addPeriod,
  Boolean convertQuotes,
  Int2 tildeAction
)

{
  FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, string, FALSE, convertQuotes, tildeAction);
  FFAddOneString (ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);

  if ( addPeriod ) {
    FFAddPeriod(ffstring);
  }
}
   
NLM_EXTERN CharPtr FFEndPrint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  FmtType format,
  Int2 gb_init_indent,
  Int2 gb_cont_indent,
  Int2 eb_init_indent,
  Int2 eb_cont_indent,
  CharPtr eb_line_prefix
)
{
  StringItemPtr temp = FFGetString(ajp);
  CharPtr result;

  if ( (ffstring == NULL) || (ajp == NULL) ) return NULL;

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    FFLineWrap(temp, ffstring, gb_init_indent, gb_cont_indent, ASN2FF_GB_MAX, NULL);
  } else {
    FFLineWrap(temp, ffstring, eb_init_indent, eb_cont_indent, ASN2FF_EMBL_MAX, eb_line_prefix);
  }
  result = FFToCharPtr(temp);
  FFRecycleString(ajp, temp);
  return result;
}

NLM_EXTERN Uint4 FFLength(StringItemPtr ffstring) {
  Uint4 len = 0;
  StringItemPtr current;

  for ( current = ffstring; current != NULL; current = current->next ) {
    len += current->pos;
  }

  return len;
}


NLM_EXTERN Char FFCharAt(StringItemPtr ffstring, Uint4 pos) {
  Uint4 inbufpos = pos % STRING_BUF_LEN;
  Uint4 count = 0;
  StringItemPtr current = NULL;

  inbufpos = pos % STRING_BUF_LEN;
  
  for ( current = ffstring; current != NULL; current = current->next ) {
    count += current->pos;
    if ( count > pos ) break;
  }

  if ( current != NULL && inbufpos <= pos )  {
    return current->buf[inbufpos];
  }

  return '\0';
}


NLM_EXTERN Char FFFindChar (
  StringItemPtr ffstring,   /* StringItem to search in */
  StringItemPtr start_buf,  /* the position of the last char searched for (buffer) */
  Uint4 start_pos,          /* the position of the last char searched for (pos) */
  Uint4 old_pos,         /* the global position searched for */
  Uint4 new_pos             /* new search position */
)
{
  Uint4 delta;
  Uint4 count;
  StringItemPtr current = NULL;

  Char result = '\0';

  if ( new_pos == old_pos ) {
    result = start_buf->buf[start_pos];
  } 

  if ( new_pos > old_pos ) {
    delta = new_pos - old_pos;
    current = start_buf;
    count = current->pos - start_pos - 1;
    current = current->next;
    
    while ( delta > count && current != NULL ) {
      current = current->next;
      count += current->pos;
    }
    
    if ( current != NULL  )  {
      result = current->buf[new_pos % STRING_BUF_LEN];
    }
    
  } else /* new_pos < old_pos */ {
    delta = old_pos - new_pos;
    if ( old_pos % STRING_BUF_LEN >= delta ) {
      result = start_buf->buf[new_pos % STRING_BUF_LEN];
    } else {
      result = FFCharAt(ffstring, new_pos);
    }
  }

  return result;
}

NLM_EXTERN Boolean FFEmpty(StringItemPtr ffstring) {
  if ( ffstring != NULL && ffstring->pos != 0 ) {
    return FALSE;
  }
  return TRUE;
}

/*
 * Compute the right-most position in the pattern at which character a occurs,
 * for each character a in the alphabet (assumed ASCII-ISO 8859-1)
 * 
 * The result is returned in the supplied vector.
 */
static void ComputeLastOccurrence(const CharPtr pattern, Uint4 last_occurrence[])
{
    Uint4 i;
    Uint4 pat_len;

    /* Initilalize vector */
    for ( i = 0; i < 256; ++i ) {
        last_occurrence[i] = 0;
    }

    /* compute right-most occurrence */
    pat_len = StringLen(pattern);
    for ( i = 0; i < pat_len; ++i ) {
        last_occurrence[(Uint1)pattern[i]] = i;
    }
}

static void ComputePrefix(const CharPtr pattern, Uint4 longest_prefix[])
{
    Uint4 pat_len = StringLen(pattern);
    Uint4 k, q;

    longest_prefix[0] = 0;

    k = 0;
    for ( q = 1; q < pat_len; ++q ) {
        while ( k > 0 && pattern[k] != pattern[q] ) {
            k = longest_prefix[k - 1];
        }
        if ( pattern[k] == pattern[q] ) {
            ++k;
        }
        longest_prefix[q] = k;
    }
}


static void ComputeGoodSuffix(const CharPtr pattern, Uint4 good_suffix[])
{
    Uint4 pat_len = StringLen(pattern);
    Uint4Ptr longest_prefix, reverse_longest_prefix;
    CharPtr reverse_pattern;
    Uint4 i, j;

    /* allocate memory */
    longest_prefix = MemNew(pat_len * sizeof(Uint4));
    reverse_longest_prefix = MemNew(pat_len * sizeof(Uint4));
    reverse_pattern = MemNew((pat_len + 1) * sizeof(Char));

    if ( longest_prefix == NULL  ||
         reverse_longest_prefix == NULL  ||
         reverse_pattern == NULL ) {
      MemFree(longest_prefix);
      MemFree(reverse_longest_prefix);
      MemFree(reverse_pattern);
      return;
    }

    /* compute reverse pattern */
    for ( i = 0; i < pat_len; ++i ) {
      reverse_pattern[pat_len - i] = pattern[i];
    }

    ComputePrefix(pattern, longest_prefix);
    ComputePrefix(reverse_pattern, reverse_longest_prefix);

    for ( j = 0; j < pat_len; ++j) {
        good_suffix[j] = pat_len - longest_prefix[pat_len-1];
    }

    for ( i = 0; i < pat_len; ++i ) {
        j = pat_len - reverse_longest_prefix[i] - 1;
        if ( good_suffix[j] > i - reverse_longest_prefix[i] + 1) {
            good_suffix[j] = i - reverse_longest_prefix[i] + 1;
        }
    }

    MemFree(longest_prefix);
    MemFree(reverse_longest_prefix);
    MemFree(reverse_pattern);
}


/*
 * searches for a pattern in a StringItem.
 * Using the Boyer-Moore algorithm for the search.
 */
NLM_EXTERN Int4 FFStringSearch (
  StringItemPtr text,
  const CharPtr pattern,
  Uint4 position )
{
  Uint4 text_len = FFLength(text);
  Uint4 pat_len = StringLen(pattern);
  Uint4 last_occurrence[256];
  Uint4Ptr good_suffix;
  Uint4 shift;
  Int4 j;

  if ( pat_len == 0 ) return 0;
  if ( text_len == 0 || pat_len > text_len - position ) return -1;
  
  good_suffix = (Uint4Ptr)MemNew(pat_len * sizeof(Int4));
  if ( good_suffix == NULL ) return -1;

  ComputeLastOccurrence(pattern, last_occurrence);
  ComputeGoodSuffix(pattern, good_suffix);

  shift = position;
  while ( shift <= text_len - pat_len ) {
    j = pat_len - 1;
    while( j >= 0 && pattern[j] == FFCharAt(text,shift + j) ) {
      --j;
    }
    if ( j == -1 ) {
      MemFree (good_suffix);
      return shift;
    } else {
        shift += MAX( (Int4)good_suffix[(int) j],
              (Int4)(j - last_occurrence[(int) FFCharAt(text,shift + j)]));
    }
  }
  MemFree (good_suffix);

  return -1;
}


/*                                                                   */
/* IsWholeWordSubstr () -- Determines if a substring that is         */
/*                         contained in another string is a whole    */
/*                         word or phrase -- i.e. is it both         */
/*                         preceded and followed by white space.     */
/*                                                                   */

NLM_EXTERN Boolean IsWholeWordSubstr (
  StringItemPtr searchStr,
  Uint4 foundPos,
  CharPtr subStr
)
{
    Boolean left, right;
    Char ch;


    /* check on the left only if there is a character there */
    if (foundPos > 0) {
        ch = FFCharAt(searchStr, foundPos - 1);
        left = IS_WHITESP(ch) || ispunct(ch);
    } else {
        left = TRUE;
    }

    foundPos += StringLen(subStr);
  if ( foundPos == FFLength(searchStr) ) {
    right = TRUE;
  } else {
    ch = FFCharAt(searchStr, foundPos);
      right = IS_WHITESP(ch) || ispunct(ch);
  }

    return left; /* see comment above */
  /* return left && right;  this is how it should be!*/
}


/* functions to record sections or blocks in linked lists */

NLM_EXTERN BaseBlockPtr Asn2gbAddBlock (
  Asn2gbWorkPtr awp,
  BlockType blocktype,
  size_t size
)

{
  BaseBlockPtr  bbp;
  ValNodePtr    vnp;

  if (awp == NULL || size < 1) return NULL;

  bbp = (BaseBlockPtr) MemNew (size);
  if (bbp == NULL) return NULL;
  bbp->blocktype = blocktype;
  bbp->section = awp->currsection;

  vnp = ValNodeAddPointer (&(awp->lastblock), 0, bbp);
  if (vnp == NULL) return bbp;

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }

  return bbp;
}


/*--------------------------------------------------------*/
/*                                                        */
/*  s_LocusGetBaseName() -                                */
/*                                                        */
/*--------------------------------------------------------*/

static Boolean s_LocusGetBaseName (BioseqPtr parent, BioseqPtr segment, CharPtr baseName)
{
  Char          parentName[SEQID_MAX_LEN];
  Char          segName[SEQID_MAX_LEN];
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;
  Char          prefix[5];
  Char          bufTmp[SEQID_MAX_LEN];
  Int2          deleteChars;
  Int2          newLength;
  Int2          i;
  Uint2         segNameLen;

  /* Get the parent Sequence ID */

  parentName [0] = '\0';
  sip = NULL;
  for (sip = parent->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
    if (sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) break;
  }

  if (sip != NULL) {
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip != NULL && (! StringHasNoText (tsip->name))) {
      StringNCpy_0 (parentName, tsip->name, sizeof (parentName));
    }
  }

  if (StringHasNoText (parentName)) {
    StringNCpy_0 (parentName, baseName, sizeof (parentName));
  }

  /* Get segment id */

  segName [0] = '\0';
  segNameLen = 0;
  sip = NULL;
  for (sip = segment->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
    if (sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) break;
    }

  if (sip != NULL) {
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip != NULL && (! StringHasNoText (tsip->name))) {
      StringNCpy_0 (segName, tsip->name, sizeof (segName));
      segNameLen = StringLen(segName);
    }
  }

  /* If there's no "SEG_" prefix, then */
  /* just use the parent ID.           */

  StringNCpy_0 (prefix,parentName,sizeof (prefix));
  prefix[4] = '\0';
  if (StringCmp(prefix,"SEG_") != 0)
    {
      StringCpy(baseName,parentName);
      return FALSE;
    }

  /* Otherwise, eliminate the "SEG_" ... */

  StringCpy(bufTmp, &parentName[4]);
  StringCpy(parentName,bufTmp);

  /* ... And calculate a base name */

  if (segNameLen > 0 &&
      (segName[segNameLen-1] == '1') &&
      (StringLen(parentName) == segNameLen) &&
      (parentName[segNameLen-1] == segName[segNameLen-1]))
    {
      deleteChars = 1;
      for (i = segNameLen-2; i >= 0; i--)
    if (parentName[i] == '0')
      deleteChars++;
    else
      break;
      newLength = segNameLen - deleteChars;
      StringNCpy (parentName,segName,newLength); /* not StringNCpy_0 */
      parentName[newLength] = '\0';
    }

  /* Return the base name in the basename parameter */

  StringCpy(baseName,parentName);
  return TRUE;
}

/* ********************************************************************** */

  static Uint1 fasta_order [NUM_SEQID] = {
    33, /* 0 = not set */
    20, /* 1 = local Object-id */
    15, /* 2 = gibbsq */
    16, /* 3 = gibbmt */
    30, /* 4 = giim Giimport-id */
    10, /* 5 = genbank */
    10, /* 6 = embl */
    10, /* 7 = pir */
    10, /* 8 = swissprot */
    15, /* 9 = patent */
    20, /* 10 = other TextSeqId */
    20, /* 11 = general Dbtag */
    255, /* 12 = gi */
    10, /* 13 = ddbj */
    10, /* 14 = prf */
    12, /* 15 = pdb */
    10, /* 16 = tpg */
    10, /* 17 = tpe */
    10, /* 18 = tpd */
    10  /* 19 = gpp */
  };

/* DoOneSection builds a single report for one bioseq or segment */

static Asn2gbSectPtr Asn2gbAddSection (
  Asn2gbWorkPtr awp
)

{
  Asn2gbSectPtr  asp;
  ValNodePtr     vnp;

  if (awp == NULL) return NULL;

  asp = (Asn2gbSectPtr) MemNew (sizeof (IntAsn2gbSect));
  if (asp == NULL) return NULL;

  vnp = ValNodeAddPointer (&(awp->lastsection), 0, asp);
  if (vnp == NULL) return asp;

  awp->lastsection = vnp;
  if (awp->sectionList == NULL) {
    awp->sectionList = vnp;
  }

  return asp;
}

NLM_EXTERN Boolean DeltaLitOnly (
  BioseqPtr bsp
)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

NLM_EXTERN Boolean SegHasParts (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return FALSE;
  sep = bsp->seqentry;
  if (sep == NULL) return FALSE;
  sep = sep->next;
  if (sep == NULL || (! IS_Bioseq_set (sep))) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) return TRUE;
  return FALSE;
}

NLM_EXTERN void DoOneSection (
  BioseqPtr target,
  BioseqPtr parent,
  BioseqPtr bsp,
  BioseqPtr refs,
  SeqLocPtr slp,
  Uint2 seg,
  Int4 from,
  Int4 to,
  Boolean contig,
  Boolean onePartOfSeg,
  Asn2gbWorkPtr awp
)

{
  size_t               acclen;
  Asn2gbFormatPtr      afp;
  IntAsn2gbJobPtr      ajp;
  Asn2gbSectPtr        asp;
  SeqMgrBioseqContext  bcontext;
  BaseBlockPtr         PNTR blockArray;
  SeqMgrDescContext    dcontext;
  Boolean              hasRefs;
  Int4                 i;
  IntAsn2gbSectPtr     iasp;
  Boolean              isRefSeq = FALSE;
  MolInfoPtr           mip;
  Boolean              nsgenome = FALSE;
  Int4                 numBlocks;
  Int4                 numsegs = 0;
  SeqDescrPtr          sdp;
  SeqIdPtr             sip;
  TextSeqIdPtr         tsip;
  ValNodePtr           vnp;
  Boolean              wgsmaster = FALSE;
  Boolean              wgstech = FALSE;
  Boolean              willshowwgs = FALSE;
  Boolean              willshowgenome = FALSE;
  Boolean              willshowcontig = FALSE;
  Boolean              willshowsequence = FALSE;

  if (target == NULL || parent == NULL || bsp == NULL || awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (awp->mode == RELEASE_MODE && awp->style == CONTIG_STYLE) {
    if (bsp->repr == Seq_repr_seg) {
    } else if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
    } else return;
  }

  if (ajp->flags.suppressLocalID) {
    sip = SeqIdSelect (bsp->id, fasta_order, NUM_SEQID);
    if (sip == NULL || sip->choice == SEQID_LOCAL) return;
  }

  if (seg == 0) {
    awp->basename[0] = '\0';
  } else if (seg == 1) {
    s_LocusGetBaseName (parent, bsp, awp->basename);
  }

  asp = Asn2gbAddSection (awp);
  if (asp == NULL) return;

  afp = awp->afp;
  if (afp != NULL) {
    afp->asp = asp;
  }

  numsegs = awp->partcount;
  if (numsegs == 0 && SeqMgrGetBioseqContext (parent, &bcontext)) {
    numsegs = bcontext.numsegs;
  }

  /* set working data fields */

  awp->asp = asp;

  awp->target = target;
  awp->parent = parent;
  awp->bsp = bsp;
  awp->refs = refs;
  awp->slp = slp;
  (awp->sectionCount)++;
  awp->currGi = 0;
  awp->seg = seg;
  awp->numsegs = numsegs;
  awp->from = from;
  awp->to = to;
  awp->contig = contig;

  awp->firstfeat = TRUE;
  awp->featseen = FALSE;
  awp->featjustseen = FALSE;
  awp->wgsaccnlist = NULL;

  /* initialize empty blockList for this section */

  awp->blockList = NULL;
  awp->lastblock = NULL;

  /* and store section data into section fields */

  asp->target = target;
  asp->bsp = bsp;
  asp->slp = slp;
  asp->seg = seg;
  asp->numsegs = numsegs;
  asp->from = from;
  asp->to = to;

  iasp = (IntAsn2gbSectPtr) asp;

  asp->blockArray = NULL;
  asp->numBlocks = 0;

  /* WGS master and NS_ virtual records treated differently */

  if (bsp->repr == Seq_repr_virtual) {

    /* check for certain ID types */

    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GENBANK ||
          sip->choice == SEQID_EMBL ||
          sip->choice == SEQID_DDBJ) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          acclen = StringLen (tsip->accession);
          if (acclen == 12) {
            if (StringCmp (tsip->accession + 6, "000000") == 0) {
              wgsmaster = TRUE;
            }
          } else if (acclen == 13) {
            if (StringCmp (tsip->accession + 6, "0000000") == 0) {
              wgsmaster = TRUE;
            }
          }
        }
     } else if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            wgsmaster = TRUE;
          } else if (StringNICmp (tsip->accession, "NS_", 3) == 0) {
            nsgenome = TRUE;
          } else if (StringNICmp (tsip->accession, "NZ_", 3) == 0) {
            if (StringLen (tsip->accession) == 15) {
              if (StringCmp (tsip->accession + 9, "000000") == 0) {
                wgsmaster = TRUE;
              }
            }
          }
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL && mip->tech == MI_TECH_wgs) {
        wgstech = TRUE;
      }
    }
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      isRefSeq = TRUE;
    } else if (sip->choice == SEQID_GI) {
      awp->currGi = (Int4) sip->data.intvalue;
    }
  }

  /* start exploring and populating paragraphs */

  if (awp->format == FTABLE_FMT) {
    AddFeatHeaderBlock (awp);
    if (awp->showFtableRefs) {
      AddReferenceBlock (awp, isRefSeq);
    }
    if (! awp->hideSources) {
      AddSourceFeatBlock (awp);
    }
    AddFeatureBlock (awp);

  } else {

    if (wgsmaster && wgstech) {
      willshowwgs = TRUE;
    } else if (nsgenome) {
      willshowgenome = TRUE;
    } else if (contig) {
      willshowcontig = TRUE;
      if (awp->showContigAndSeq) {
        if (! awp->hideSequence) {
          willshowsequence = TRUE;
        }
      }
    } else {
      if (awp->showContigAndSeq) {
        if (bsp->repr == Seq_repr_seg && (! SegHasParts (bsp))) {
          willshowcontig = TRUE;
        } else if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
          willshowcontig = TRUE;
        }
      }
      if (! awp->hideSequence) {
        willshowsequence = TRUE;
      }
    }

    AddLocusBlock (awp, willshowwgs, willshowgenome, willshowcontig, willshowsequence);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

      AddDeflineBlock (awp);
      AddAccessionBlock (awp);

      if (ISA_aa (bsp->mol)) {
        /*
        AddPidBlock (awp);
        */
      }

      AddVersionBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddProjectBlock (awp);
      }

      if (ISA_aa (bsp->mol)) {
        AddDbsourceBlock (awp);
      }

    } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

      AddAccessionBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddVersionBlock (awp);
      }

      if (ISA_aa (bsp->mol)) {
        /* AddPidBlock (awp); */
        /* AddDbsourceBlock (awp); */
      }

      AddDateBlock (awp);

      AddDeflineBlock (awp);
    }

    AddKeywordsBlock (awp);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
      AddSegmentBlock (awp, onePartOfSeg);
    }

    AddSourceBlock (awp);
    AddOrganismBlock (awp);

    if (! awp->hidePubs) {

      /* !!! RELEASE_MODE should check return value of AddReferenceBlock !!! */

      hasRefs = AddReferenceBlock (awp, isRefSeq);
      if (! hasRefs) {
        if (ajp->flags.needAtLeastOneRef) {
          /* RefSeq does not require a publication */
          if (! isRefSeq) {
            awp->failed = TRUE;
          }
        }
      }
    }

    AddCommentBlock (awp);
    AddPrimaryBlock (awp);

    AddFeatHeaderBlock (awp);
    if (! awp->hideSources) {
      AddSourceFeatBlock (awp);
    }

    if (wgsmaster && wgstech) {

      AddWGSBlock (awp);

    } else if (nsgenome) {

      AddGenomeBlock (awp);

    } else if (contig) {

      if (awp->showconfeats) {
        AddFeatureBlock (awp);
      } else if (awp->smartconfeats && bsp->length <= 1000000) {
        AddFeatureBlock (awp);
      }
      AddContigBlock (awp);

      if (awp->showContigAndSeq) {
        if (ISA_na (bsp->mol) && ajp->gbseq == NULL) {
          if (awp->showBaseCount) {
            AddBasecountBlock (awp);
          }
        }
        AddOriginBlock (awp);

        if (! awp->hideSequence) {
          AddSequenceBlock (awp);
        }
      }

    } else {

      AddFeatureBlock (awp);

      if (awp->showContigAndSeq) {
        if (bsp->repr == Seq_repr_seg && (! SegHasParts (bsp))) {
          AddContigBlock (awp);
        } else if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
          AddContigBlock (awp);
        }
      }

      if (ISA_na (bsp->mol) && ajp->gbseq == NULL) {
        if (awp->showBaseCount) {
          AddBasecountBlock (awp );
        }
      }
      AddOriginBlock (awp);

      if (! awp->hideSequence) {
        AddSequenceBlock (awp);
      }
    }

    AddSlashBlock (awp);
  }

  /* allocate block array for this section */

  numBlocks = ValNodeLen (awp->blockList);
  asp->numBlocks = numBlocks;

  if (numBlocks > 0) {
    blockArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numBlocks + 1));
    asp->blockArray = blockArray;

    if (blockArray != NULL) {
      for (vnp = awp->blockList, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        blockArray [i] = (BaseBlockPtr) vnp->data.ptrvalue;
      }
    }
  }

  /* free blockList, but leave data, now pointed to by blockArray elements */

  awp->blockList = ValNodeFree (awp->blockList);
  awp->lastblock = NULL;

  (awp->currsection)++;
}

/* ********************************************************************** */

/*
the following functions handle various kinds of input, all calling
DoOneSection once for each component that gets its own report
*/

static Boolean LIBCALLBACK Asn2Seg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
  BioseqPtr      parent;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  parent = context->parent;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* may remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    parent = bsp;
    from = 0;
    to = bsp->length - 1;
  }

  if (bsp->repr != Seq_repr_virtual) {
    (awp->seg)++;
    DoOneSection (bsp, parent, bsp, bsp, /* slp */ NULL, awp->seg, from, to, FALSE, FALSE, awp);
  }

  BioseqUnlock (bsp);

  return TRUE;
}

static Int4 CountRealParts (
  SeqLocPtr slp_head
)

{
  SeqIdPtr   id;
  Int4       numparts;
  BioseqPtr  part;
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  numparts = 0;
  for (slp = (SeqLocPtr) slp_head; slp != NULL; slp = slp->next) {
    sip = SeqLocId (slp);
    if (sip == NULL) continue;
    if (sip->choice == SEQID_GI) {
      part = BioseqFind (sip);
      if (part == NULL) continue;
      for (id = part->id; id != NULL; id = id->next) {
        if (id->choice == SEQID_GIBBSQ ||
            id->choice == SEQID_GIBBMT ||
            id->choice == SEQID_GIIM) break;
      }
      if (id != NULL && part->repr == Seq_repr_virtual) continue;
    }
    numparts++;
  }
  return numparts;
}

typedef struct findseg {
  BioseqPtr  bsp;
  Uint2      seg;
} FindSeg, PNTR FindSegPtr;

static Boolean LIBCALLBACK FindSegForPart (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  FindSegPtr  fsp;
  BioseqPtr   bsp = NULL;
  SeqLocPtr   loc;
  SeqIdPtr    sip;

  if (slp == NULL || context == NULL) return TRUE;
  fsp = (FindSegPtr) context->userdata;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFind (sip);
  if (bsp == NULL) return TRUE;

  if (bsp->repr != Seq_repr_virtual) {
    (fsp->seg)++;
  }

  if (bsp != fsp->bsp) return TRUE;

  return FALSE;
}

NLM_EXTERN void DoOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  IntAsn2gbJobPtr       ajp;
  Asn2gbWorkPtr         awp;
  BioseqSetPtr          bssp;
  SeqMgrSegmentContext  context;
  Boolean               contig = FALSE;
  Int4                  from;
  FindSeg               fs;
  SeqEntryPtr           oldscope;
  BioseqPtr             parent;
  Boolean               segmented = FALSE;
  SeqEntryPtr           sep;
  Int4                  to;

  if (bsp == NULL) return;
  awp = (Asn2gbWorkPtr) userdata;
  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  /* return if molecule not right for format */

  if (ISA_na (bsp->mol)) {
    if (ajp->format == GENPEPT_FMT || ajp->format == EMBLPEPT_FMT) return;

    /* only do mRNA feature tables in GPS if targeted to a specific mRNA */

    if (ajp->format == FTABLE_FMT && ajp->skipMrnas) {
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bsp->idx.parentptr;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
          if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bsp->idx.parentptr;
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
              return;
            }
          }
        }
      }
    }

  } else if (ISA_aa (bsp->mol)) {
    if (ajp->format == GENBANK_FMT || ajp->format == EMBL_FMT) return;

    /* only do protein feature tables if targeted to a specific protein */

    if (ajp->format == FTABLE_FMT && ajp->skipProts) return;
  }

  if (awp->style == SEGMENT_STYLE) {
    segmented = TRUE;
  }
  if (awp->style == CONTIG_STYLE) {
    contig = TRUE;
  }
  /* Never do segmented style in FTABLE format */
  if (awp->format == FTABLE_FMT) {
      segmented = FALSE;
      contig = FALSE;
  }

  awp->partcount = 0;

  if (bsp->repr == Seq_repr_seg && awp->style == NORMAL_STYLE) {

    /* if bsp followed by parts set, then do not default to contig style */

    if (SegHasParts (bsp)) {
      segmented = TRUE;
      contig = FALSE;

      if (bsp->seq_ext_type == 1) {

        /* count only non-virtual parts */

        sep = GetTopSeqEntryForEntityID (awp->entityID);
        oldscope = SeqEntrySetScope (sep);
        awp->partcount = CountRealParts ((SeqLocPtr) bsp->seq_ext);
        SeqEntrySetScope (oldscope);
      }
    } else {
      segmented = FALSE;
      contig = TRUE;
    }
  }
  if (bsp->repr == Seq_repr_delta && awp->style == NORMAL_STYLE) {
    if (! DeltaLitOnly (bsp)) {
      contig = TRUE;
    }
  }

  if (bsp->repr == Seq_repr_seg) {

    /* this is a segmented bioseq */

    if (segmented) {

      /* show all segments individually */

      awp->seg = 0;
      SeqMgrExploreSegments (bsp, (Pointer) awp, Asn2Seg);

    } else {

      /* show as single bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (parent, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }

  } else if (bsp->repr == Seq_repr_raw ||
             bsp->repr == Seq_repr_const ||
             bsp->repr == Seq_repr_delta ||
             bsp->repr == Seq_repr_virtual) {

    parent = SeqMgrGetParentOfPart (bsp, &context);
    if (parent != NULL) {

      /* this is a part of an indexed segmented bioseq */

      from = context.cumOffset;
      to = from + context.to - context.from;

      s_LocusGetBaseName (parent, bsp, awp->basename);

      fs.bsp = bsp;
      fs.seg = 0;
      SeqMgrExploreSegments (parent, (Pointer) &fs, FindSegForPart);
      awp->showAllFeats = TRUE;

      DoOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, fs.seg, from, to, contig, TRUE, awp);

    } else {

      /* this is a regular non-segmented bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }
  }
}

static void DoBioseqSetList (
  SeqEntryPtr seq_set,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (seq_set == NULL || awp == NULL) return;

  /* iterate rather than recurse unless multiple nested sets > nuc-prot */

  for (sep = seq_set; sep != NULL; sep = sep->next) {

    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) continue;

      if (bssp->_class == BioseqseqSet_class_genbank ||
          bssp->_class == BioseqseqSet_class_mut_set ||
          bssp->_class == BioseqseqSet_class_pop_set ||
          bssp->_class == BioseqseqSet_class_phy_set ||
          bssp->_class == BioseqseqSet_class_eco_set ||
          bssp->_class == BioseqseqSet_class_wgs_set ||
          bssp->_class == BioseqseqSet_class_gen_prod_set) {

        /* if popset within genbank set, for example, recurse */

        DoBioseqSetList (bssp->seq_set, awp);

        continue;
      }
    }

    /* at most nuc-prot set, so do main bioseqs that fit the format */

    VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, DoOneBioseq);
  }
}

static void DoOneBioseqSet (
  SeqEntryPtr sep,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL || awp == NULL) return;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;

    if (bssp->_class == BioseqseqSet_class_genbank ||
        bssp->_class == BioseqseqSet_class_mut_set ||
        bssp->_class == BioseqseqSet_class_pop_set ||
        bssp->_class == BioseqseqSet_class_phy_set ||
        bssp->_class == BioseqseqSet_class_eco_set ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_gen_prod_set) {

      /* this is a pop/phy/mut/eco set, catenate separate reports */

      DoBioseqSetList (bssp->seq_set, awp);

      return;
    }
  }

  /* at most nuc-prot set, so do main bioseqs that fit the format */

  VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, DoOneBioseq);
}

/* ********************************************************************** */

static void RecordOneSection (
  Asn2gbWorkPtr awp,
  BioseqPtr bsp,
  SeqIdPtr sip
)

{
  IntAsn2gbJobPtr  ajp;
  ValNodePtr       vnp;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) break;
    }
  }

  if (sip == NULL) return;

  if (sip->choice == SEQID_GI) {
    vnp = ValNodeAddInt (&(ajp->gitail), 0, (Int4) sip->data.intvalue);
    if (ajp->gihead == NULL) {
      ajp->gihead = vnp;
    }
    ajp->gitail = vnp;
  }

  (awp->sectionMax)++;
}

static void CountOneSection (
  BioseqPtr target,
  BioseqPtr parent,
  BioseqPtr bsp,
  BioseqPtr refs,
  SeqLocPtr slp,
  Uint2 seg,
  Int4 from,
  Int4 to,
  Boolean contig,
  Boolean onePartOfSeg,
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr  ajp;
  SeqIdPtr         sip;

  if (target == NULL || parent == NULL || bsp == NULL || awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (awp->mode == RELEASE_MODE && awp->style == CONTIG_STYLE) {
    if (bsp->repr == Seq_repr_seg) {
    } else if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
    } else return;
  }

  if (ajp->flags.suppressLocalID) {
    sip = SeqIdSelect (bsp->id, fasta_order, NUM_SEQID);
    if (sip == NULL || sip->choice == SEQID_LOCAL) return;
  }

  RecordOneSection (awp, bsp, NULL);
}


static Boolean LIBCALLBACK Count2Seg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Int4           from;
  SeqLocPtr      loc;
  BioseqPtr      parent;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  parent = context->parent;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFindCore (sip);
  if (bsp != NULL && bsp->repr == Seq_repr_virtual) return TRUE;

  RecordOneSection (awp, NULL, sip);

  return TRUE;
}

static void CountOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  IntAsn2gbJobPtr       ajp;
  Asn2gbWorkPtr         awp;
  BioseqSetPtr          bssp;
  SeqMgrSegmentContext  context;
  Boolean               contig = FALSE;
  Int4                  from;
  BioseqPtr             parent;
  Boolean               segmented = FALSE;
  Int4                  to;

  if (bsp == NULL) return;
  awp = (Asn2gbWorkPtr) userdata;
  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (ISA_na (bsp->mol)) {
    if (ajp->format == GENPEPT_FMT || ajp->format == EMBLPEPT_FMT) return;

    if (ajp->format == FTABLE_FMT && ajp->skipMrnas) {
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bsp->idx.parentptr;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
          if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bsp->idx.parentptr;
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
              return;
            }
          }
        }
      }
    }

  } else if (ISA_aa (bsp->mol)) {
    if (ajp->format == GENBANK_FMT || ajp->format == EMBL_FMT) return;

    if (ajp->format == FTABLE_FMT && ajp->skipProts) return;
  }

  if (awp->style == SEGMENT_STYLE) {
    segmented = TRUE;
  }
  if (awp->style == CONTIG_STYLE) {
    contig = TRUE;
  }
  if (awp->format == FTABLE_FMT) {
      segmented = FALSE;
      contig = FALSE;
  }

  if (bsp->repr == Seq_repr_seg && awp->style == NORMAL_STYLE) {

    if (SegHasParts (bsp)) {
      segmented = TRUE;
      contig = FALSE;
    } else {
      segmented = FALSE;
      contig = TRUE;
    }
  }
  if (bsp->repr == Seq_repr_delta && awp->style == NORMAL_STYLE) {
    if (! DeltaLitOnly (bsp)) {
      contig = TRUE;
    }
  }

  if (bsp->repr == Seq_repr_seg) {

    if (segmented) {

      SeqMgrExploreSegments (bsp, (Pointer) awp, Count2Seg);

    } else {

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      CountOneSection (parent, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }

  } else if (bsp->repr == Seq_repr_raw ||
             bsp->repr == Seq_repr_const ||
             bsp->repr == Seq_repr_delta ||
             bsp->repr == Seq_repr_virtual) {

    parent = SeqMgrGetParentOfPart (bsp, &context);
    if (parent != NULL) {

      from = context.cumOffset;
      to = from + context.to - context.from;

      CountOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, TRUE, awp);

    } else {

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      CountOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }
  }
}



static void CountBioseqSetList (
  SeqEntryPtr seq_set,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (seq_set == NULL || awp == NULL) return;

  for (sep = seq_set; sep != NULL; sep = sep->next) {

    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) continue;

      if (bssp->_class == BioseqseqSet_class_genbank ||
          bssp->_class == BioseqseqSet_class_mut_set ||
          bssp->_class == BioseqseqSet_class_pop_set ||
          bssp->_class == BioseqseqSet_class_phy_set ||
          bssp->_class == BioseqseqSet_class_eco_set ||
          bssp->_class == BioseqseqSet_class_wgs_set ||
          bssp->_class == BioseqseqSet_class_gen_prod_set) {

        CountBioseqSetList (bssp->seq_set, awp);

        continue;
      }
    }

    VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, CountOneBioseq);
  }
}

static void CountOneBioseqSet (
  SeqEntryPtr sep,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL || awp == NULL) return;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;

    if (bssp->_class == BioseqseqSet_class_genbank ||
        bssp->_class == BioseqseqSet_class_mut_set ||
        bssp->_class == BioseqseqSet_class_pop_set ||
        bssp->_class == BioseqseqSet_class_phy_set ||
        bssp->_class == BioseqseqSet_class_eco_set ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_gen_prod_set) {

      CountBioseqSetList (bssp->seq_set, awp);

      return;
    }
  }

  VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, CountOneBioseq);
}

/* ********************************************************************** */

/* public functions */

static int LIBCALLBACK SortParagraphByIDProc (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  BaseBlockPtr  bbp1, bbp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  bbp1 = *((BaseBlockPtr PNTR) vp1);
  bbp2 = *((BaseBlockPtr PNTR) vp2);
  if (bbp1 == NULL || bbp2 == NULL) return 0;

  if (bbp1->entityID > bbp2->entityID) return 1;
  if (bbp1->entityID < bbp2->entityID) return -1;

  if (bbp1->itemtype > bbp2->itemtype) return 1;
  if (bbp1->itemtype < bbp2->itemtype) return -1;

  if (bbp1->itemID > bbp2->itemID) return 1;
  if (bbp1->itemID < bbp2->itemID) return -1;

  if (bbp1->paragraph > bbp2->paragraph) return 1;
  if (bbp1->paragraph < bbp2->paragraph) return -1;

  return 0;
}

static void IsBspRefseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BoolPtr   has_refseqP;
  SeqIdPtr  sip;

  if (bsp == NULL || userdata == NULL) return;
  has_refseqP = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      *has_refseqP = TRUE;
    }
  }
}

static Boolean IsSepRefseq (
  SeqEntryPtr sep
)

{
  Boolean  is_refseq = FALSE;

  if (sep == NULL) return FALSE;
  VisitBioseqsInSep (sep, (Pointer) &is_refseq, IsBspRefseq);
  return is_refseq;
}

typedef struct modeflags {
  Boolean  flags [29];
} ModeFlags, PNTR ModeFlagsPtr;

static ModeFlags flagTable [] = {

  /* RELEASE_MODE */
  {TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  TRUE},

  /* ENTREZ_MODE */
  {FALSE, TRUE,  TRUE,  TRUE,  TRUE,
   FALSE, TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  FALSE, TRUE,  TRUE,
   TRUE,  TRUE,  FALSE, FALSE, TRUE,
   TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
   TRUE,  TRUE,  TRUE,  FALSE},

  /* SEQUIN_MODE */
  {FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, TRUE,  FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, TRUE,  FALSE, FALSE},

  /* DUMP_MODE */
  {FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE}
};

static void SetFlagsFromMode (
  IntAsn2gbJobPtr ajp,
  ModType mode
)

{
  BoolPtr       bp;
  ModeFlagsPtr  mfp;
  SeqEntryPtr   sep;

  if (ajp == NULL) return;
  if (! (mode >= RELEASE_MODE && mode <= DUMP_MODE)) {
    mode = DUMP_MODE;
  }
  mfp = &(flagTable [(int) (mode - 1)]);
  bp = &(mfp->flags [0]);

  ajp->flags.suppressLocalID = *(bp++);
  ajp->flags.validateFeats = *(bp++);
  ajp->flags.ignorePatPubs = *(bp++);
  ajp->flags.dropShortAA = *(bp++);
  ajp->flags.avoidLocusColl = *(bp++);

  ajp->flags.iupacaaOnly = *(bp++);
  ajp->flags.dropBadCitGens = *(bp++);
  ajp->flags.noAffilOnUnpub = *(bp++);
  ajp->flags.dropIllegalQuals = *(bp++);
  ajp->flags.checkQualSyntax = *(bp++);

  ajp->flags.needRequiredQuals = *(bp++);
  ajp->flags.needOrganismQual = *(bp++);
  ajp->flags.needAtLeastOneRef = *(bp++);
  ajp->flags.citArtIsoJta = *(bp++);
  ajp->flags.dropBadDbxref = *(bp++);

  ajp->flags.useEmblMolType = *(bp++);
  ajp->flags.hideBankItComment = *(bp++);
  ajp->flags.checkCDSproductID = *(bp++);
  ajp->flags.suppressSegLoc = *(bp++);
  ajp->flags.srcQualsToNote = *(bp)++;

  ajp->flags.hideEmptySource = *(bp++);
  ajp->flags.goQualsToNote = *(bp++);
  ajp->flags.geneSynsToNote = *(bp++);
  ajp->flags.refSeqQualsToNote = *(bp++);
  ajp->flags.selenocysteineToNote = *(bp++);

  ajp->flags.pyrrolysineToNote = *(bp++);
  ajp->flags.extraProductsToNote = *(bp++);
  ajp->flags.codonRecognizedToNote = *(bp++);
  ajp->flags.forGbRelease = *(bp++);

  /* unapproved qualifiers suppressed for flatfile, okay for GBSeq XML */

  if (ajp->gbseq == NULL) {

    /* collaboration unapproved source quals on their own line only in indexer Sequin - relaxed */

    /*
    if (GetAppProperty ("InternalNcbiSequin") == NULL) {

      ajp->flags.srcQualsToNote = TRUE;
    }
    */

    sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
    if (IsSepRefseq (sep)) {

      ajp->flags.srcQualsToNote = FALSE;
      ajp->flags.goQualsToNote = FALSE;
      ajp->flags.refSeqQualsToNote = FALSE;

      /* selenocysteine always a separate qualifier for RefSeq */

      ajp->flags.selenocysteineToNote = FALSE;
      ajp->flags.pyrrolysineToNote = FALSE;

    } else {

      /* collaboration unapproved Gene Ontology quals on their own line only for RefSeq */

      /* ajp->flags.goQualsToNote = TRUE; */
      ajp->flags.geneSynsToNote = TRUE;
    }

  } else {

    sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
    if (IsSepRefseq (sep)) {

      ajp->flags.srcQualsToNote = FALSE;

      /* selenocysteine always a separate qualifier for RefSeq */

      ajp->flags.selenocysteineToNote = FALSE;
      ajp->flags.pyrrolysineToNote = FALSE;

    }
  }

  if (ajp->refseqConventions) {
    ajp->flags.srcQualsToNote = FALSE;
    ajp->flags.goQualsToNote = FALSE;
    ajp->flags.refSeqQualsToNote = FALSE;
  }
}

static void CheckVersionWithGi (BioseqPtr bsp, Pointer userdata)

{
  Boolean       hasGi = FALSE;
  BoolPtr       missingVersion;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;
  Boolean       zeroVersion = FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_TPG:
      case SEQID_TPE:
      case SEQID_TPD:
      case SEQID_GENBANK:
      case SEQID_EMBL:
      case SEQID_DDBJ:
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->version == 0) {
          zeroVersion = TRUE;
        }
        break;
      case SEQID_GI :
        hasGi = TRUE;
        break;
      default :
        break;
    }
  }
  if (hasGi && zeroVersion) {
    missingVersion = (BoolPtr) userdata;
    *missingVersion = TRUE;
  }
}


typedef struct lookforids {
  Boolean isG;
  Boolean isGED;
  Boolean isNTorNW;
  Boolean isNC;
  Boolean isRefSeq;
  Boolean isGeneral;
  Boolean isTPA;
  Boolean isTPG;
  Boolean isNuc;
  Boolean isProt;
  Boolean isLocal;
  Boolean isNonLocal;
} LookForIDs, PNTR LookForIDsPtr;

static void LookForSeqIDs (BioseqPtr bsp, Pointer userdata)

{
  DbtagPtr       dbt;
  LookForIDsPtr  lfip;
  SeqIdPtr       sip;
  TextSeqIdPtr   tsip;

  lfip = (LookForIDsPtr) userdata;
  if (ISA_na (bsp->mol)) {
    lfip->isNuc = TRUE;
  }
  if (ISA_aa (bsp->mol)) {
    lfip->isProt = TRUE;
  }
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GENBANK :
        lfip->isG = TRUE;
        /* and fall through to EMBL and DDBJ */
      case SEQID_EMBL :
      case SEQID_DDBJ :
        lfip->isGED = TRUE;
        lfip->isNonLocal = TRUE;
        break;
      case SEQID_TPG :
        lfip->isTPG = TRUE;
        /* and fall through to TPE and TPD */
      case SEQID_TPE :
      case SEQID_TPD :
        lfip->isTPA = TRUE;
        lfip->isNonLocal = TRUE;
        break;
      case SEQID_OTHER :
        lfip->isRefSeq = TRUE;
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL) {
          if (StringNCmp (tsip->accession, "NC_", 3) == 0) {
            lfip->isNC = TRUE;
          } else if (StringNCmp (tsip->accession, "NT_", 3) == 0) {
            lfip->isNTorNW = TRUE;
          } else if (StringNCmp (tsip->accession, "NW_", 3) == 0) {
            lfip->isNTorNW = TRUE;
          }
        }
        lfip->isNonLocal = TRUE;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "TMSMART") != 0 && StringICmp (dbt->db, "BankIt") != 0) {
            lfip->isGeneral = TRUE;
            lfip->isNonLocal = TRUE;
          }
        }
        break;
      case SEQID_LOCAL :
        lfip->isLocal = TRUE;
        break;
      default :
        lfip->isNonLocal = TRUE;
        break;
    }
  }
}

static void LookForGEDetc (
  SeqEntryPtr topsep,
  BoolPtr isG,
  BoolPtr isGED,
  BoolPtr isNTorNW,
  BoolPtr isNC,
  BoolPtr isRefSeq,
  BoolPtr isGeneral,
  BoolPtr isTPA,
  BoolPtr isTPG,
  BoolPtr isNuc,
  BoolPtr isProt,
  BoolPtr isOnlyLocal
)

{
  LookForIDs  lfi;

  MemSet ((Pointer) &lfi, 0, sizeof (LookForIDs));
  VisitBioseqsInSep (topsep, (Pointer) &lfi, LookForSeqIDs);
  *isG = lfi.isG;
  *isGED = lfi.isGED;
  *isNTorNW = lfi.isNTorNW;
  *isNC = lfi.isNC;
  *isRefSeq = lfi.isRefSeq;
  *isGeneral = lfi.isGeneral;
  *isTPA = lfi.isTPA;
  *isTPG = lfi.isTPG;
  *isNuc = lfi.isNuc;
  *isProt = lfi.isProt;
  if (lfi.isLocal && (! lfi.isNonLocal)) {
    *isOnlyLocal = TRUE;
  } else {
    *isOnlyLocal = FALSE;
  }
}

static void MakeGapFeats (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Char             buf [32];
  Int4             currpos = 0;
  BioseqPtr        fakebsp = NULL;
  IntFuzzPtr       fuzz;
  ValNodePtr PNTR  gapvnp;
  ImpFeatPtr       ifp;
  SeqLitPtr        litp;
  SeqAnnotPtr      sap = NULL;
  SeqFeatPtr       sfp;
  SeqIdPtr         sip;
  SeqLocPtr        slp;
  ValNodePtr       vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return;
  gapvnp = (ValNodePtr PNTR) userdata;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  /* suppress on far delta contigs for now */
  if (! DeltaLitOnly (bsp)) return;

  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      if (slp == NULL) continue;
      currpos += SeqLocLen (slp);
    }
    if (vnp->choice == 2) {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL) continue;
      if (litp->seq_data == NULL && litp->length > 0) {
        if (fakebsp == NULL) {
          /* to be freed with MemFree, not BioseqFree */
          fakebsp = MemNew (sizeof (Bioseq));
          if (fakebsp == NULL) return;
          sap = SeqAnnotNew ();
          if (sap == NULL) return;
          sap->type = 1;
          fakebsp->annot = sap;
          ValNodeAddPointer (gapvnp, 0, (Pointer) fakebsp);
        }
        ifp = ImpFeatNew ();
        if (ifp == NULL) continue;
        ifp->key = StringSave ("gap");
        sfp = SeqFeatNew ();
        if (sfp == NULL) continue;
        sfp->data.choice = SEQFEAT_IMP;
        sfp->data.value.ptrvalue = (Pointer) ifp;
        sfp->next = (SeqFeatPtr) sap->data;
        sap->data = (Pointer) sfp;
        fuzz = litp->fuzz;
        if (fuzz != NULL && fuzz->choice == 4 && fuzz->a == 0) {
          AddQualifierToFeature (sfp, "estimated_length", "unknown");
          sfp->location = AddIntervalToLocation (NULL, sip, currpos, currpos + litp->length - 1, FALSE, FALSE);
        } else {
          sprintf (buf, "%ld", (long) litp->length);
          AddQualifierToFeature (sfp, "estimated_length", buf);
          sfp->location = AddIntervalToLocation (NULL, sip, currpos, currpos + litp->length - 1, FALSE, FALSE);
        }
      }
      currpos += litp->length;
    }
  }
}

static void LookForFeatFetchPolicy (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BoolPtr        forceOnlyNearFeatsP;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  forceOnlyNearFeatsP = (BoolPtr) userdata;
  if (forceOnlyNearFeatsP == NULL) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "FeatureFetchPolicy") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
    if (StringCmp (oip->str, "Policy") == 0) {
      if (StringICmp ((CharPtr) ufp->data.ptrvalue, "OnlyNearFeatures") == 0) {
        *forceOnlyNearFeatsP = TRUE;
      }
    }
  }
}

static CharPtr bad_html_strings [] = {
  "<script", "<object", "<applet", "<embed", "<form", "javascript:", NULL
};

static CharPtr defHead = "\
<html>\n\
<head>\n\
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=us-ascii\" />\
<title>GenBank entry</title>\n\
</head>\n\
<body>\n\
<hr />\n\
<pre>";

static CharPtr defTail = "\
</pre>\n\
<hr />\n\
</body>\n\
</html>\n";

#define FAR_TRANS_MASK (SHOW_FAR_TRANSLATION | TRANSLATE_IF_NO_PRODUCT | ALWAYS_TRANSLATE_CDS)
#define FEAT_FETCH_MASK (ONLY_NEAR_FEATURES | FAR_FEATURES_SUPPRESS | NEAR_FEATURES_SUPPRESS)
#define HTML_XML_ASN_MASK (CREATE_HTML_FLATFILE | CREATE_XML_GBSEQ_FILE | CREATE_ASN_GBSEQ_FILE)
#define PUBLICATION_MASK (HIDE_GENE_RIFS | ONLY_GENE_RIFS | ONLY_REVIEW_PUBS | NEWEST_PUBS | OLDEST_PUBS | HIDE_ALL_PUBS)

static Asn2gbJobPtr asn2gnbk_setup_ex (
  BioseqPtr bsp,
  BioseqSetPtr bssp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  Boolean stream,
  FILE *fp,
  AsnIoPtr aip,
  AsnTypePtr atp
)

{
  Asn2gbFormat     af;
  IntAsn2gbJobPtr  ajp = NULL;
  Asn2gbSectPtr    asp;
  Asn2gbWork       aw;
  BaseBlockPtr     bbp;
  BaseBlockPtr     PNTR blockArray;
  Uint2            eID = 0;
  Uint2            entityID = 0;
  CharPtr          ffhead = NULL;
  CharPtr          fftail = NULL;
  Asn2gbWriteFunc  ffwrite = NULL;
  Boolean          forceOnlyNearFeats = FALSE;
  ValNodePtr       gapvnp = NULL;
  GBSeqPtr         gbseq = NULL;
  Int4             i;
  IndxPtr          index = NULL;
  Boolean          isG;
  Boolean          isGED;
  Boolean          isGeneral;
  Boolean          isNTorNW;
  Boolean          isNC;
  Boolean          isNuc;
  Boolean          isOnlyLocal;
  Boolean          isProt;
  Boolean          isRefSeq;
  Boolean          isTPA;
  Boolean          isTPG;
  Int4             j;
  Int4             k;
  SeqLocPtr        loc = NULL;
  Boolean          lockFarComp;
  Boolean          lockFarLocs;
  Boolean          lockFarProd;
  Boolean          lookupFarComp;
  Boolean          lookupFarHist;
  Boolean          lookupFarInf;
  Boolean          lookupFarLocs;
  Boolean          lookupFarOthers;
  Boolean          lookupFarProd;
  Boolean          missingVersion;
  Int4             nextGi = 0;
  Boolean          noLeft;
  Boolean          noRight;
  Int4             numBlocks;
  Int4             numSections;
  SeqEntryPtr      oldscope;
  ObjMgrDataPtr    omdp;
  Int4             numParagraphs;
  BaseBlockPtr     PNTR paragraphArray;
  BaseBlockPtr     PNTR paragraphByIDs;
  BioseqPtr        parent = NULL;
  Int4             prevGi = 0;
  Int2             q;
  Pointer          remotedata = NULL;
  Asn2gbFreeFunc   remotefree = NULL;
  Asn2gbLockFunc   remotelock = NULL;
  ValNodePtr       remotevnp = NULL;
  Asn2gbSectPtr    PNTR sectionArray;
  SubmitBlockPtr   sbp;
  SeqEntryPtr      sep;
  SeqIntPtr        sintp;
  SeqIdPtr         sip;
  Boolean          skipMrnas = FALSE;
  Boolean          skipProts = FALSE;
  SeqSubmitPtr     ssp;
  BioseqSetPtr     topbssp;
  Pointer          userdata = NULL;
  ValNodePtr       vnp;
  Boolean          is_html = FALSE;

  if (format == 0) {
    format = GENBANK_FMT;
  }
  if (mode == 0) {
    mode = SEQUIN_MODE;
  }
  if (style == 0) {
    style = NORMAL_STYLE;
  }

  if (extra != NULL) {
    ffwrite = extra->ffwrite;
    ffhead = extra->ffhead;
    fftail = extra->fftail;
    index = extra->index;
    gbseq = extra->gbseq;
    userdata = extra->userdata;
    remotelock = extra->remotelock;
    remotefree = extra->remotefree;
    remotedata = extra->remotedata;
    prevGi = extra->prevGi;
    nextGi = extra->nextGi;
  }

  if (slp != NULL) {
    sip = SeqLocId (slp);
    bsp = BioseqFind (sip);
    if (bsp == NULL) {
      bsp = BioseqFindFromSeqLoc (slp);
    }
    if (bsp == NULL) return NULL;

    /* if location is on part of segmented set, need to map to segmented bioseq */

    if (slp->choice == SEQLOC_WHOLE) {
      /* Entrez server may pass in whole location on part instead of part bioseq */
      slp = NULL;
    } else if (sip == NULL) {
      parent = bsp;
    } else {

      /* SeqMgrGetParentOfPart depends upon feature indexing */

      eID = ObjMgrGetEntityIDForPointer (bsp);
      if (SeqMgrFeaturesAreIndexed (eID) == 0) {
        SeqMgrIndexFeatures (eID, NULL);
      }

      parent = SeqMgrGetParentOfPart (bsp, NULL);
    }
    if (parent != NULL) {
      CheckSeqLocForPartial (slp, &noLeft, &noRight);
      loc = SeqLocMergeEx (parent, slp, NULL, FALSE, TRUE, FALSE, FALSE);
      slp = loc;
      FreeAllFuzz (slp);
      SetSeqLocPartial (slp, noLeft, noRight);
    }

    /* if location is whole, generate normal bioseq report */

    if (slp == NULL) {
      /* reality check in case SeqLocMergeEx fails and sets slp to NULL, or if was cleared above */
    } else if (slp->choice == SEQLOC_WHOLE) {
      slp = NULL;
      SeqLocFree (loc);
      loc = NULL;
    } else if (slp->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL &&
          sintp->from == 0 &&
          sintp->to == bsp->length - 1 &&
          sintp->strand == Seq_strand_plus) {
        slp = NULL;
        SeqLocFree (loc);
        loc = NULL;
      }
    }
  }

  if (bsp != NULL) {
    bssp = NULL;
    entityID = ObjMgrGetEntityIDForPointer (bsp);
  } else if (bssp != NULL) {
    entityID = ObjMgrGetEntityIDForPointer (bssp);
    if (format == FTABLE_FMT) {
      skipProts = TRUE;
      skipMrnas = TRUE;
    }
  }
  if ((Boolean) ((custom & SHOW_PROT_FTABLE) != 0)) {
    skipProts = FALSE;
    skipMrnas = FALSE;
  }

  if (entityID == 0) return NULL;

  sep = GetTopSeqEntryForEntityID (entityID);

  LookForGEDetc (sep, &isG, &isGED, &isNTorNW, &isNC, &isRefSeq,
                 &isGeneral, &isTPA, &isTPG, &isNuc, &isProt, &isOnlyLocal);

  if (mode == RELEASE_MODE) {
    missingVersion = FALSE;
    VisitBioseqsInSep (sep, (Pointer) &missingVersion, CheckVersionWithGi);
    if (missingVersion) return NULL;
  }

  ajp = (IntAsn2gbJobPtr) MemNew (sizeof (IntAsn2gbJob));
  if (ajp == NULL) return NULL;

  VisitDescriptorsInSep (sep, (Pointer) &forceOnlyNearFeats, LookForFeatFetchPolicy);

  gapvnp = NULL;
  if (format != FTABLE_FMT) {
    if (isGED /* was isG */ || isTPG || isOnlyLocal || isRefSeq || (isGeneral && (! isGED))) {
      if ((Boolean) ((custom & HIDE_GAP_FEATS) == 0)) {
        VisitBioseqsInSep (sep, (Pointer) &gapvnp, MakeGapFeats);
      }
    }
  }
  ajp->gapvnp = gapvnp;

  remotevnp = NULL;
  ajp->remotelock = remotelock;
  ajp->remotefree = remotefree;
  ajp->remotedata = remotedata;
  if (remotelock != NULL && bsp != NULL) {
    sip = SeqIdFindBest (bsp->id, SEQID_GI);
    if (sip != NULL) {
      remotevnp = remotelock (sip, remotedata);
    }
  }
  ajp->remotevnp = remotevnp;

  if (gapvnp != NULL || remotevnp != NULL) {
    /* if both gapvnp and remotevnp, link together so everything is indexed */
    if (gapvnp != NULL) {
      ValNodeLink(&gapvnp, remotevnp);
    } else {
      gapvnp = remotevnp;
    }
    SeqMgrClearFeatureIndexes (entityID, NULL);
    SeqMgrIndexFeaturesExEx (entityID, NULL, FALSE, FALSE, gapvnp);
    gapvnp->next = NULL;
  }

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  is_html = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_HTML_FLATFILE);
  if (is_html) {
    InitWWW(ajp);
  }

  ajp->ajp.entityID = entityID;
  ajp->ajp.bsp = bsp;
  ajp->ajp.bssp = bssp;

  if (loc != NULL) {
    ajp->ajp.slp = loc;
  } else if (slp != NULL) {
    ajp->ajp.slp = AsnIoMemCopy ((Pointer) slp,
                                 (AsnReadFunc) SeqLocAsnRead,
                                 (AsnWriteFunc) SeqLocAsnWrite);
  } else {
    ajp->ajp.slp = NULL;
  }

  /* if location specified, normal defaults to master style */

  if (ajp->ajp.slp != NULL && style == NORMAL_STYLE) {
    style = MASTER_STYLE;
  }

  ajp->format = format;
  ajp->mode = mode; /* for showing new qualifiers before quarantine ends */

  ajp->index = index;
  ajp->gbseq = gbseq; /* gbseq output can relax srcQualsToNote or goQualsToNote strictness */
  ajp->aip = aip;
  ajp->atp = atp;

  ajp->refseqConventions = (Boolean) ((flags & REFSEQ_CONVENTIONS) != 0);

  SetFlagsFromMode (ajp, mode);

  lockFarComp = (Boolean) ((locks & LOCK_FAR_COMPONENTS) != 0);
  lockFarLocs = (Boolean) ((locks & LOCK_FAR_LOCATIONS) != 0);
  lockFarProd = (Boolean) ((locks & LOCK_FAR_PRODUCTS) != 0);

  if (lockFarComp || lockFarLocs || lockFarProd) {

    /* lock all bioseqs in advance, including remote genome components */

    if (ajp->ajp.slp != NULL && lockFarComp) {
      ajp->lockedBspList = LockFarComponentsEx (sep, FALSE, lockFarLocs, lockFarProd, ajp->ajp.slp);
    } else {
      ajp->lockedBspList = LockFarComponentsEx (sep, lockFarComp, lockFarLocs, lockFarProd, NULL);
    }
  }

  lookupFarComp = (Boolean) ((locks & LOOKUP_FAR_COMPONENTS) != 0);
  lookupFarLocs = (Boolean) ((locks & LOOKUP_FAR_LOCATIONS) != 0);
  lookupFarProd = (Boolean) ((locks & LOOKUP_FAR_PRODUCTS) != 0);
  lookupFarHist = (Boolean) ((locks & LOOKUP_FAR_HISTORY) != 0);
  lookupFarInf = (Boolean) ((locks & LOOKUP_FAR_INFERENCE) != 0);
  lookupFarOthers = (Boolean) ((locks & LOOKUP_FAR_OTHERS) != 0);

  if (lookupFarComp || lookupFarLocs || lookupFarProd || lookupFarHist || lookupFarInf || lookupFarOthers) {

    /* lookukp all far SeqIDs in advance */

    LookupFarSeqIDs (sep, lookupFarComp, lookupFarLocs, lookupFarProd, FALSE, lookupFarHist, lookupFarInf, lookupFarOthers);
  }

  ajp->showFarTransl = (Boolean) ((flags & FAR_TRANS_MASK) == SHOW_FAR_TRANSLATION);
  ajp->transIfNoProd = (Boolean) ((flags & FAR_TRANS_MASK) == TRANSLATE_IF_NO_PRODUCT);
  ajp->alwaysTranslCds = (Boolean) ((flags & FAR_TRANS_MASK) == ALWAYS_TRANSLATE_CDS);
  if (ajp->transIfNoProd || ajp->alwaysTranslCds) {
    ajp->showFarTransl = TRUE;
  }

  ajp->masterStyle = (Boolean) (style == MASTER_STYLE);

  ajp->showTranscript = (Boolean) ((flags & SHOW_TRANCRIPTION) != 0);
  ajp->showPeptide = (Boolean) ((flags & SHOW_PEPTIDE) != 0);

  ajp->produceInsdSeq = (Boolean) ((flags & PRODUCE_OLD_GBSEQ) == 0);

  if (stream && (format == GENBANK_FMT || format == GENPEPT_FMT)) {
    ajp->specialGapFormat = (Boolean) ((flags & SPECIAL_GAP_DISPLAY) != 0);
    if (is_html && mode == ENTREZ_MODE) {
      ajp->specialGapFormat = TRUE;
    }
    if ((custom & EXPANDED_GAP_DISPLAY) != 0) {
      ajp->specialGapFormat = FALSE;
    }
  } else {
    ajp->specialGapFormat = FALSE;
  }
  ajp->seqGapCurrLen = 0;

  ajp->gihead = NULL;
  ajp->gitail = NULL;

  ajp->hideGoTerms = (Boolean) ((custom & HIDE_GO_TERMS) != 0);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ajp->newSourceOrg = TRUE;
  }

  ajp->relModeError = FALSE;
  ajp->skipProts = skipProts;
  ajp->skipMrnas = skipMrnas;

  MemSet ((Pointer) (&aw), 0, sizeof (Asn2gbWork));
  aw.ajp = ajp;
  aw.entityID = entityID;

  aw.sectionList = NULL;
  aw.lastsection = NULL;

  aw.currsection = 0;
  aw.showAllFeats = FALSE;

  aw.showconfeats = (Boolean) ((flags & SHOW_CONTIG_FEATURES) != 0);
  aw.showconsource = (Boolean) ((flags & SHOW_CONTIG_SOURCES) != 0);

  aw.format = format;
  aw.mode = mode;
  aw.style = style;

  /* sectionCount used for hyperlinks */

  aw.sectionCount = 0;
  aw.sectionMax = 0;
  aw.gilistpos = NULL;

  aw.currGi = 0;
  aw.prevGi = prevGi;
  aw.nextGi = nextGi;

  /* internal format pointer if writing at time of creation */

  if (stream) {
    MemSet ((Pointer) &af, 0, sizeof (Asn2gbFormat));
    af.ajp = ajp;
    af.asp = NULL;
    af.qvp = NULL;
    af.format = format;
    af.ffwrite = ffwrite;
    af.userdata = userdata;
    af.fp = fp;
    af.aip = aip;
    af.atp = atp;

    aw.afp = &af;
  }

  /* special types of records override feature fetching and contig display parameters */

  if (mode == ENTREZ_MODE) {
    if (! aw.showconfeats) {
      aw.smartconfeats = TRUE;  /* features suppressed if CONTIG style and length > 1 MB */
      aw.showconfeats = FALSE;
      aw.showconsource = FALSE;
    }
  }

  aw.onlyNearFeats = FALSE;
  aw.farFeatsSuppress = FALSE;
  aw.nearFeatsSuppress = FALSE;

  if (isNC) {

    if ((Boolean) ((flags & FEAT_FETCH_MASK) == ONLY_NEAR_FEATURES)) {
      aw.onlyNearFeats = TRUE;
    } else if (forceOnlyNearFeats) {
      aw.onlyNearFeats = TRUE;
    } else {
      aw.nearFeatsSuppress = TRUE;
    }

  } else if (isNTorNW || isTPA) {

    aw.onlyNearFeats = TRUE;

  } else if (isGED) {

    if ((Boolean) ((flags & FEAT_FETCH_MASK) == ONLY_NEAR_FEATURES)) {
      aw.onlyNearFeats = TRUE;
    } else if (forceOnlyNearFeats) {
      aw.onlyNearFeats = TRUE;
    } else {
      aw.nearFeatsSuppress = TRUE;
    }
    ajp->showFarTransl = TRUE;

  } else if (forceOnlyNearFeats) {

    aw.onlyNearFeats = TRUE;
    
  } else {

    aw.onlyNearFeats = (Boolean) ((flags & FEAT_FETCH_MASK) == ONLY_NEAR_FEATURES);
    aw.farFeatsSuppress = (Boolean) ((flags & FEAT_FETCH_MASK) == FAR_FEATURES_SUPPRESS);
    aw.nearFeatsSuppress = (Boolean) ((flags & FEAT_FETCH_MASK) == NEAR_FEATURES_SUPPRESS);
  }

  /* continue setting flags */

  aw.hideImpFeats = (Boolean) ((custom & HIDE_IMP_FEATS) != 0);
  aw.hideRemImpFeats = (Boolean) ((custom & HIDE_REM_IMP_FEATS) != 0);

  aw.hideSnpFeats = (Boolean) ((custom & HIDE_SNP_FEATS) != 0);
  aw.hideExonFeats = (Boolean) ((custom & HIDE_EXON_FEATS) != 0);
  aw.hideIntronFeats = (Boolean) ((custom & HIDE_INTRON_FEATS) != 0);
  aw.hideMiscFeats = (Boolean) ((custom & HIDE_MISC_FEATS) != 0);

  aw.hideCddFeats = (Boolean) ((custom & HIDE_CDD_FEATS) != 0);
  aw.hideCdsProdFeats = (Boolean) ((custom & HIDE_CDS_PROD_FEATS) != 0);

  aw.hideGeneRIFs = (Boolean) ((custom & PUBLICATION_MASK) == HIDE_GENE_RIFS);
  aw.onlyGeneRIFs = (Boolean) ((custom & PUBLICATION_MASK) == ONLY_GENE_RIFS);
  aw.onlyReviewPubs = (Boolean) ((custom & PUBLICATION_MASK) == ONLY_REVIEW_PUBS);
  aw.newestPubs = (Boolean) ((custom & PUBLICATION_MASK) == NEWEST_PUBS);
  aw.oldestPubs = (Boolean) ((custom & PUBLICATION_MASK) == OLDEST_PUBS);
  aw.hidePubs = (Boolean) ((custom & PUBLICATION_MASK) == HIDE_ALL_PUBS);

  aw.showFtableRefs = (Boolean) ((custom & SHOW_FTABLE_REFS) != 0);
  aw.hideSources = (Boolean) ((custom & HIDE_SOURCE_FEATS) != 0);
  aw.hideGaps = (Boolean) ((custom & HIDE_GAP_FEATS) != 0);
  aw.hideSequence = (Boolean) ((custom & HIDE_SEQUENCE) != 0);

  aw.isGPS = FALSE;
  if (sep != NULL && IS_Bioseq_set (sep)) {
    topbssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (topbssp != NULL && topbssp->_class == BioseqseqSet_class_gen_prod_set) {
      aw.isGPS = TRUE;
      aw.copyGpsCdsUp = (Boolean) ((flags & COPY_GPS_CDS_UP) != 0);
      aw.copyGpsGeneDown = (Boolean) ((flags & COPY_GPS_GENE_DOWN) != 0);
    }
  }

  aw.showContigAndSeq = (Boolean) ((flags & SHOW_CONTIG_AND_SEQ) != 0);
  /*
  if (style != MASTER_STYLE && style != SEGMENT_STYLE) {
    aw.showContigAndSeq = FALSE;
  }
  */

  aw.newLocusLine = TRUE;
  aw.showBaseCount = FALSE;

  if ((Boolean) ((flags & DDBJ_VARIANT_FORMAT) != 0)) {
    aw.citSubsFirst = TRUE;
    aw.hideGeneFeats = TRUE;
    aw.newLocusLine = FALSE;
    aw.showBaseCount = TRUE;
    ajp->newSourceOrg = FALSE;
  }
  if (mode == SEQUIN_MODE || mode == DUMP_MODE) {
    aw.showBaseCount = TRUE;
  }
  aw.forcePrimaryBlock = (Boolean) ((flags & FORCE_PRIMARY_BLOCK) != 0);

  aw.localFeatCount = VisitFeaturesInSep (sep, NULL, NULL);

  aw.hup = FALSE;
  aw.ssp = NULL;

  aw.failed = FALSE;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      aw.ssp = ssp;
      sbp = ssp->sub;
      if (sbp != NULL) {
        aw.hup = sbp->hup;
      }
    }
  }

  ajp->bad_html_fsa = TextFsaNew ();

  for (q = 0; bad_html_strings [q] != NULL; q++) {
    TextFsaAdd (ajp->bad_html_fsa, bad_html_strings [q]);
  }

  oldscope = SeqEntrySetScope (sep);

  if (stream) {
    /* send optional head string */

    is_html = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_HTML_FLATFILE);
    if (ffhead == NULL && is_html) {
      ffhead = defHead;
    }
    if (ffhead != NULL) {
      if (fp != NULL) {
        fprintf (fp, ffhead);
      }
    }
    if (ffwrite != NULL) {
      ffwrite (ffhead, userdata, HEAD_BLOCK);
    }
    if (is_html) {
      DoQuickLinkFormat (aw.afp, "<div class=\"sequence\">");
    }
  }

  /* if Web Entrez, set awp->sectionMax to decide when Next hyperlink is needed */

  if (is_html && mode == ENTREZ_MODE && stream &&
      (format == GENBANK_FMT || format == GENPEPT_FMT)) {
     /* add dummy node as prev id for first section */
     ajp->gihead = ValNodeAddInt (&(ajp->gitail), 0, (Int4) 0);
     ajp->gitail = ajp->gihead;
    if (bssp != NULL) {
      CountOneBioseqSet (SeqMgrGetSeqEntryForData (bssp), &aw);
    } else {
      CountOneBioseq (bsp, &aw);
    }
  }

  if (bssp != NULL) {

    /* handle all components of a pop/phy/mut/eco set */

    sep = SeqMgrGetSeqEntryForData (bssp);
    DoOneBioseqSet (sep, &aw);

  } else {

    /* handle single bioseq, which may be segmented or a local part */

    DoOneBioseq (bsp, &aw);
  }

  if (stream) {
    if (is_html) {
      DoQuickLinkFormat (aw.afp, "</div>");
    }

    /* send optional tail string */

    if (fftail == NULL && is_html) {
      fftail = defTail;
    }
    if (fftail != NULL) {
      if (fp != NULL) {
        fprintf (fp, fftail);
      }
    }
    if (ffwrite != NULL) {
      ffwrite (fftail, userdata, TAIL_BLOCK);
    }
  }

  SeqEntrySetScope (oldscope);

  /* check for failure to populate anything */

  numSections = ValNodeLen (aw.sectionList);
  ajp->ajp.numSections = numSections;

  if (numSections == 0) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  /* allocate section array for this job - needed for memory cleanup even if streamed */

  sectionArray = (Asn2gbSectPtr PNTR) MemNew (sizeof (Asn2gbSectPtr) * (numSections + 1));
  ajp->ajp.sectionArray = sectionArray;

  if (sectionArray == NULL) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  /* fill in section and paragraph arrays */

  numParagraphs = 0;
  for (vnp = aw.sectionList, i = 0; vnp != NULL && i < numSections; vnp = vnp->next, i++) {
    asp = (Asn2gbSectPtr) vnp->data.ptrvalue;
    sectionArray [i] = asp;
    if (asp != NULL) {
      numParagraphs += asp->numBlocks;
    }
  }

  /* allocate paragraph array pointing to all blocks in all sections */

  ajp->ajp.numParagraphs = numParagraphs;
  if (numParagraphs == 0) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  paragraphArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->ajp.paragraphArray = paragraphArray;

  paragraphByIDs = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->ajp.paragraphByIDs = paragraphByIDs;

  if (paragraphArray == NULL || paragraphByIDs == NULL) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  k = 0;
  for (i = 0; i < numSections; i++) {
    asp = sectionArray [i];
    if (asp != NULL) {

      numBlocks = asp->numBlocks;
      blockArray = asp->blockArray;
      if (blockArray != NULL) {

        for (j = 0; j < numBlocks; j++) {
          bbp = blockArray [j];

          paragraphArray [k] = bbp;
          paragraphByIDs [k] = bbp;
          bbp->paragraph = k;
          k++;
        }
      }
    }
  }

  /* sort paragraphByIDs array by entityID/itemtype/itemID/paragraph */

  HeapSort (paragraphByIDs, (size_t) numParagraphs, sizeof (BaseBlockPtr), SortParagraphByIDProc);

  /* free sectionList, but leave data, now pointed to by sectionArray elements */

  ValNodeFree (aw.sectionList);

  /* check for failure to to make legal flatfile */

  if (ajp->flags.needAtLeastOneRef && aw.failed) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  return (Asn2gbJobPtr) ajp;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_setup (
  BioseqPtr bsp,
  BioseqSetPtr bssp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra
)

{
  return asn2gnbk_setup_ex (bsp, bssp, slp, format, mode, style,
                            flags, locks, custom, extra,
                            FALSE, NULL, NULL, NULL);
}

/* ********************************************************************** */

/* format functions allocate printable string for given paragraph */

NLM_EXTERN CharPtr DefaultFormatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  if (afp == NULL || bbp == NULL) return NULL;

  /* default format function assumes string pre-allocated by add block function */

  return StringSaveNoNull (bbp->string);
}

typedef CharPtr (*FormatProc) (Asn2gbFormatPtr afp, BaseBlockPtr bbp);

static FormatProc asn2gnbk_fmt_functions [28] = {
  NULL,
  NULL,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  FormatSourceBlock,
  FormatOrganismBlock,
  FormatReferenceBlock,
  DefaultFormatBlock,
  FormatCommentBlock,
  FormatFeatHeaderBlock,
  FormatSourceFeatBlock,
  FormatFeatureBlock,
  FormatBasecountBlock,
  DefaultFormatBlock,
  FormatSequenceBlock,
  FormatContigBlock,
  DefaultFormatBlock,
  DefaultFormatBlock,
  FormatSlashBlock,
  NULL
};

NLM_EXTERN void PrintFtableIntervals (
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqLocPtr location,
  CharPtr label
)

{
  IntFuzzPtr  ifp;
  Boolean     partial5;
  Boolean     partial3;
  SeqLocPtr   slp;
  SeqPntPtr   spp;
  Int4        start;
  Int4        stop;
  Char        str [64];
  Char        str1 [32];
  Char        str2 [32];

  if (head == NULL || target == NULL || location == NULL || label == NULL) return;

  if (location->choice == SEQLOC_PNT) {
    spp = (SeqPntPtr) location->data.ptrvalue;
    if (spp != NULL) {
      ifp = spp->fuzz;
      if (ifp != NULL && ifp->choice == 4 && ifp->a == 3) {
        sprintf (str, "%ld^\t%ld\t%s\n", (long) (spp->point + 1),
                 (long) (spp->point + 2), label);
        ValNodeCopyStr (head, 0, str);
        return;
      }
    }
  }

  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return;

  start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  if (partial5) {
    sprintf (str1, "<%ld", (long) start);
  } else {
    sprintf (str1, "%ld", (long) start);
  }
  if (partial3) {
    sprintf (str2, ">%ld", (long) stop);
  } else {
    sprintf (str2, "%ld", (long) stop);
  }
  sprintf (str, "%s\t%s\t%s\n", str1, str2, label);
  ValNodeCopyStr (head, 0, str);

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (partial5) {
      sprintf (str1, "<%ld", (long) start);
    } else {
      sprintf (str1, "%ld", (long) start);
    }
    if (partial3) {
      sprintf (str2, ">%ld", (long) stop);
    } else {
      sprintf (str2, "%ld", (long) stop);
    }
    if (start != 0 && stop != 0) {
      sprintf (str, "%s\t%s\n", str1, str2);
      ValNodeCopyStr (head, 0, str);
    }
  }
}

static CharPtr goQualList [] = {
  "", "go_process", "go_component", "go_function", NULL
};

static void PrintFTUserFld (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr     entry;
  CharPtr          evidence = NULL;
  Char             gid [32];
  CharPtr          goid = NULL;
  ValNodePtr PNTR  head;
  Int2             i;
  Int2             j;
  size_t           len;
  ObjectIdPtr      oip;
  Int4             pmid = 0;
  CharPtr          str;
  CharPtr          textstr = NULL;
  Char             tmp [16];

  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; goQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, goQualType [i]) == 0) break;
  }
  if (goQualType [i] == NULL) return;

  /* loop to allow multiple entries for each type of GO term */
  for (entry = ufp->data.ptrvalue; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
  
    pmid = 0;
    goid = NULL;
    evidence = NULL;
    textstr = NULL;

    for (ufp = (UserFieldPtr) entry->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; goFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, goFieldType [j]) == 0) break;
      }
      if (goFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }
    /* if (StringHasNoText (textstr)) break; */
  
    str = (CharPtr) MemNew (StringLen (textstr) + StringLen (goid) + StringLen (evidence) + 40);
    if (str == NULL) return;
    StringCpy (str, "\t\t\t");
    StringCat (str, goQualList [i]);
    StringCat (str, "\t");
    StringCat (str, textstr);
    if (! StringHasNoText (goid)) {
      StringCat (str, "|");
      StringCat (str, goid);
    } else {
      StringCat (str, "|");
    }
    if (pmid != 0) {
      sprintf (tmp, "|%ld", (long) pmid);
      StringCat (str, tmp);
    } else {
      StringCat (str, "|");
    }
    if (! StringHasNoText (evidence)) {
      StringCat (str, "|");
      StringCat (str, evidence);
    }
    len = StringLen (str);
    while (len > 0 && str [len - 1] == '|') {
      str [len - 1] = '\0';
      len--;
    }
  
    head = (ValNodePtr PNTR) userdata;
    StringCat (str, "\n");
    ValNodeCopyStr (head, 0, str);
  }
}

static void PrintFTUserObj (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "GeneOntology") != 0) return;
  VisitUserFieldsInUop (uop, userdata, PrintFTUserFld);
}

static void PrintFTCodeBreak (
  ValNodePtr PNTR head,
  CodeBreakPtr cbp,
  BioseqPtr target
)

{
  Char             buf [80];
  Choice           cbaa;
  IntAsn2gbJob     iaj;
  CharPtr          ptr;
  Uint1            residue;
  SeqCodeTablePtr  sctp;
  Uint1            seqcode;
  SeqLocPtr        slp;
  CharPtr          str;

  seqcode = 0;
  sctp = NULL;
  cbaa = cbp->aa;
  switch (cbaa.choice) {
    case 1 :
      seqcode = Seq_code_ncbieaa;
      break;
    case 2 :
      seqcode = Seq_code_ncbi8aa;
      break;
    case 3 :
      seqcode = Seq_code_ncbistdaa;
      break;
    default :
      break;
  }
  if (seqcode == 0) return;
  sctp = SeqCodeTableFind (seqcode);
  if (sctp == NULL) return;

  MemSet ((Pointer) &iaj, 0, sizeof (IntAsn2gbJob));
  iaj.flags.iupacaaOnly = FALSE;
  iaj.relModeError = FALSE;

  slp = SeqLocFindNext (cbp->loc, NULL);
  while (slp != NULL) {
    str = FFFlatLoc (&iaj, target, slp, FALSE);
    if (str != NULL) {
      residue = cbaa.value.intvalue;
      ptr = Get3LetterSymbol (&iaj, seqcode, sctp, residue);
      if (ptr == NULL) {
        ptr = "OTHER";
      }
      sprintf (buf, "\t\t\ttransl_except\t(pos:%s,aa:%s)\n", str, ptr);
      ValNodeCopyStr (head, 0, buf);
      MemFree (str);
    }
    slp = SeqLocFindNext (cbp->loc, slp);
  }
}

static Boolean SeqIdWriteForTable (SeqIdPtr sip, CharPtr buf, size_t buflen, IntAsn2gbJobPtr ajp, Boolean giOK)

{
  SeqIdPtr  accn = NULL, local = NULL, patent = NULL,
            pdb = NULL, general = NULL, gi = NULL;
  DbtagPtr  dbt;
  Char      id [41], str [64];
  Int2      numids;
  CharPtr   prefix = NULL;

  if (sip == NULL || buf == NULL || ajp == NULL) return FALSE;

  while (sip != NULL) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        local = sip;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_PIR :
      case SEQID_SWISSPROT :
      case SEQID_DDBJ :
      case SEQID_PRF :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
      case SEQID_OTHER :
      case SEQID_GPIPE :
        accn = sip;
        break;
      case SEQID_PATENT :
        patent = sip;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "TMSMART") != 0 && StringICmp (dbt->db, "BankIt") != 0) {
            general = sip;
          }
        }
        break;
      case SEQID_PDB :
        pdb = sip;
        break;
      case SEQID_GI :
        gi = sip;
        break;
      default :
        break;
    }
    sip = sip->next;
  }

  str [0] = '\0';
  numids = 0;

  if (accn != NULL) {
    if (SeqIdWrite (accn, id, PRINTID_FASTA_SHORT, sizeof (id) - 1) != NULL) {
      StringCat (str, prefix);
      StringCat (str, id);
      prefix = "|";
      numids++;
    }
  }

  if (general != NULL) {
    if (SeqIdWrite (general, id, PRINTID_FASTA_SHORT, sizeof (id) - 1) != NULL) {
      StringCat (str, prefix);
      StringCat (str, id);
      prefix = "|";
      numids++;
    }
  }

  if (local != NULL && (! ajp->flags.suppressLocalID) && numids == 0) {
    if (SeqIdWrite (local, id, PRINTID_FASTA_SHORT, sizeof (id) - 1) != NULL) {
      StringCat (str, prefix);
      StringCat (str, id);
      prefix = "|";
      numids++;
    }
  }

  if (gi != NULL && giOK && numids == 0) {
    if (SeqIdWrite (accn, id, PRINTID_FASTA_SHORT, sizeof (id) - 1) != NULL) {
      StringCat (str, prefix);
      StringCat (str, id);
      prefix = "|";
      numids++;
    }
  }

  StringNCpy_0 (buf, str, buflen);
  if (StringHasNoText (buf)) return FALSE;

  return TRUE;
}

static void PrintBioSourceFtableEntry (
  ValNodePtr PNTR head,
  BioSourcePtr    biop
)
{
  OrgModPtr    mod;
  SubSourcePtr ssp;
  Char         str [256];

  if (head == NULL || biop == NULL) return;

  if (biop->org != NULL && ! StringHasNoText (biop->org->taxname))
  {
    sprintf (str, "\t\t\torganism\t%s\n", biop->org->taxname);
    ValNodeCopyStr (head, 0, str);
  }

  /* add OrgMods */
  if (biop->org != NULL && biop->org->orgname != NULL)
  {
    for (mod = biop->org->orgname->mod;
         mod != NULL;
         mod = mod->next)
    {
      switch (mod->subtype)
      {
        case ORGMOD_strain :
          sprintf (str, "\t\t\tstrain\t");
          break;
        case ORGMOD_substrain :
          sprintf (str, "\t\t\tsubstrain\t");
          break;
        case ORGMOD_type :
          sprintf (str, "\t\t\ttype\t");
          break;
        case ORGMOD_subtype :
          sprintf (str, "\t\t\tsubtype\t");
          break;
        case ORGMOD_variety :
          sprintf (str, "\t\t\tvariety\t");
          break;
        case ORGMOD_serotype :
          sprintf (str, "\t\t\tserotype\t");
          break;
        case ORGMOD_serogroup :
          sprintf (str, "\t\t\tserogroup\t");
          break;
        case ORGMOD_serovar :
          sprintf (str, "\t\t\tserovar\t");
          break;
        case ORGMOD_cultivar :
          sprintf (str, "\t\t\tcultivar\t");
          break;
        case ORGMOD_pathovar :
          sprintf (str, "\t\t\tpathovar\t");
          break;
        case ORGMOD_chemovar :
          sprintf (str, "\t\t\tchemovar\t");
          break;
        case ORGMOD_biovar :
          sprintf (str, "\t\t\tbiovar\t");
          break;
        case ORGMOD_biotype :
          sprintf (str, "\t\t\tbiotype\t");
          break;
        case ORGMOD_group :
          sprintf (str, "\t\t\tgroup\t");
          break;
        case ORGMOD_subgroup :
          sprintf (str, "\t\t\tsubgroup\t");
          break;
        case ORGMOD_isolate :
          sprintf (str, "\t\t\tisolate\t");
          break;
        case ORGMOD_common :
          sprintf (str, "\t\t\tcommon\t");
          break;
        case ORGMOD_acronym :
          sprintf (str, "\t\t\tacronym\t");
          break;
        case ORGMOD_dosage :
          sprintf (str, "\t\t\tdosage\t");
          break;
        case ORGMOD_nat_host :
          sprintf (str, "\t\t\tnat_host\t");
          break;
        case ORGMOD_sub_species :
          sprintf (str, "\t\t\tsub_species\t");
          break;
        case ORGMOD_specimen_voucher :
          sprintf (str, "\t\t\tspecimen_voucher\t");
          break;
        case ORGMOD_authority :
          sprintf (str, "\t\t\tauthority\t");
          break;
        case ORGMOD_forma :
          sprintf (str, "\t\t\tforma\t");
          break;
        case ORGMOD_forma_specialis :
          sprintf (str, "\t\t\tforma_specialis\t");
          break;
        case ORGMOD_ecotype :
          sprintf (str, "\t\t\tecotype\t");
          break;
        case ORGMOD_synonym :
          sprintf (str, "\t\t\tsynonym\t");
          break;
        case ORGMOD_anamorph :
          sprintf (str, "\t\t\tanamorph\t");
          break;
        case ORGMOD_teleomorph :
          sprintf (str, "\t\t\tteleomorph\t");
          break;
        case ORGMOD_breed :
          sprintf (str, "\t\t\tbreed\t");
          break;
        case ORGMOD_gb_acronym :
          sprintf (str, "\t\t\tgb_acronym\t");
          break;
        case ORGMOD_gb_anamorph :
          sprintf (str, "\t\t\tgb_anamorph\t");
          break;
        case ORGMOD_old_lineage :
          sprintf (str, "\t\t\told_lineage\t");
          break;
        case ORGMOD_old_name :
          sprintf (str, "\t\t\told_name\t");
          break;
        case ORGMOD_other :
          sprintf (str, "\t\t\tnote\t");
          break;
        default :
          str [0] = 0;
      }
      if ( str [0] == 0) continue;
      if (! StringHasNoText (mod->subname))
      {
        StringNCat (str, mod->subname, sizeof (str) - StringLen (str) - 2);
        str [sizeof (str) - 2] = 0;
      }
      StringCat (str, "\n");
      ValNodeCopyStr (head, 0, str);
    }
  }
   
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next)
  {
    switch (ssp->subtype)
    {
      case SUBSRC_chromosome :
        sprintf (str, "\t\t\tchromosome\t");
        break;
      case SUBSRC_map :
        sprintf (str, "\t\t\tmap\t");
        break;
      case SUBSRC_clone :
        sprintf (str, "\t\t\tclone\t");
        break;
      case SUBSRC_haplotype :
        sprintf (str, "\t\t\thaplotype\t");
        break;
      case SUBSRC_genotype :
        sprintf (str, "\t\t\tgenotype\t");
        break;
      case SUBSRC_sex :
        sprintf (str, "\t\t\tsex\t");
        break;
      case SUBSRC_cell_line :
        sprintf (str, "\t\t\tcell_line\t");
        break;
      case SUBSRC_cell_type :
        sprintf (str, "\t\t\tcell_type\t");
        break;
      case SUBSRC_tissue_type :
        sprintf (str, "\t\t\ttissue_type\t");
        break;
      case SUBSRC_clone_lib :
        sprintf (str, "\t\t\tclone_lib\t");
        break;
      case SUBSRC_dev_stage :
        sprintf (str, "\t\t\tdev_stage\t");
        break;
      case SUBSRC_frequency :
        sprintf (str, "\t\t\tfrequency\t");
        break;
      case SUBSRC_germline :
        sprintf (str, "\t\t\tgermline\t");
        break;
      case SUBSRC_rearranged :
        sprintf (str, "\t\t\trearranged\t");
        break;
      case SUBSRC_lab_host :
        sprintf (str, "\t\t\tlab_host\t");
        break;
      case SUBSRC_pop_variant :
        sprintf (str, "\t\t\tpop_variant\t");
        break;
      case SUBSRC_tissue_lib :
        sprintf (str, "\t\t\ttissue_lib\t");
        break;
      case SUBSRC_plasmid_name :
        sprintf (str, "\t\t\tplasmid_name\t");
        break;
      case SUBSRC_transposon_name :
        sprintf (str, "\t\t\ttransposon_name\t");
        break;
      case SUBSRC_insertion_seq_name :
        sprintf (str, "\t\t\tinsertion_seq_name\t");
        break;
      case SUBSRC_plastid_name :
        sprintf (str, "\t\t\tplastid_name\t");
        break;
      case SUBSRC_country :
        sprintf (str, "\t\t\tcountry\t");
        break;
      case SUBSRC_segment :
        sprintf (str, "\t\t\tsegment\t");
        break;
      case SUBSRC_endogenous_virus_name :
        sprintf (str, "\t\t\tendogenous_virus_name\t");
        break;
      case SUBSRC_transgenic :
        sprintf (str, "\t\t\ttransgenic\t");
        break;
      case SUBSRC_environmental_sample :
        sprintf (str, "\t\t\tenvironmental_sample\t");
        break;
      case SUBSRC_isolation_source :
        sprintf (str, "\t\t\tisolation_source\t");
        break;
      case SUBSRC_lat_lon :
        sprintf (str, "\t\t\tlat_lon\t");
        break;
      case SUBSRC_collection_date :
        sprintf (str, "\t\t\tcollection_date\t");
        break;
      case SUBSRC_collected_by :
        sprintf (str, "\t\t\tcollected_by\t");
        break;
      case SUBSRC_identified_by :
        sprintf (str, "\t\t\tidentified_by\t");
        break;
      case SUBSRC_fwd_primer_seq :
        sprintf (str, "\t\t\tfwd_pcr_primer_seq\t");
        break;
      case SUBSRC_rev_primer_seq :
        sprintf (str, "\t\t\trev_pcr_primer_seq\t");
        break;
      case SUBSRC_fwd_primer_name :
        sprintf (str, "\t\t\tfwd_pcr_primer_name\t");
        break;
      case SUBSRC_rev_primer_name :
        sprintf (str, "\t\t\trev_pcr_primer_name\t");
        break;
      case SUBSRC_other :
          sprintf (str, "\t\t\tnote\t");
          break;
      default :
        str [0] = 0;
    }
    if ( str [0] == 0) continue;
    if (! StringHasNoText (ssp->name))
    {
      StringNCat (str, ssp->name, sizeof (str) - StringLen (str) - 2);
      str [sizeof (str) - 2] = 0;
    }
    StringCat (str, "\n");
    ValNodeCopyStr (head, 0, str);
  }  
}

NLM_EXTERN void PrintFtableLocAndQuals (
  IntAsn2gbJobPtr ajp,
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  CharPtr            aa;
  Int2               bondidx;
  BioseqSetPtr       bssp;
  CodeBreakPtr       cbp;
  BioseqPtr          cdna;
  SeqFeatPtr         cds;
  CdRegionPtr        crp;
  SeqMgrDescContext  dcontext;
  DbtagPtr           dbt;
  SeqMgrFeatContext  fcontext;
  GBQualPtr          gbq;
  ValNodePtr         geneorprotdb;
  GeneRefPtr         grp;
  Boolean            is_gps_genomic = FALSE;
  CharPtr            label;
  MolInfoPtr         mip;
  SeqLocPtr          newloc;
  Char               numbuf [32];
  Int2               numcodons;
  ObjectIdPtr        oip;
  BioseqPtr          prod;
  SeqFeatPtr         prot;
  ProtRefPtr         prp;
  Boolean            pseudo;
  RnaRefPtr          rrp;
  SeqDescrPtr        sdp;
  Int4               sec_str;
  SeqIdPtr           sip;
  SeqIdPtr           sip2;
  Int2               siteidx;
  SeqLocPtr          slp;
  Char               str [256];
  Char               tmp [300];
  CharPtr            tmpx;
  tRNAPtr            trp;
  ValNodePtr         vnp;

  if (head == NULL || target == NULL || sfp == NULL || context == NULL) return;
  /* label = (CharPtr) FeatDefTypeLabel (sfp); */
  label = FindKeyFromFeatDefType (sfp->idx.subtype, FALSE);
  if (StringCmp (label, "Gene") == 0) {
    label = "gene";
  }
  else if (StringCmp (label, "Src") == 0) {
    label = "source";
  }
  if (StringHasNoText (label)) {
    label = "???";
  }

  /* check if genomic sequence in genomic product set */

  if (target->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) target->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      sdp = SeqMgrGetNextDescriptor (target, NULL, Seq_descr_molinfo, &dcontext);
      if (sdp != NULL) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip != NULL && mip->biomol == MOLECULE_TYPE_GENOMIC) {
          is_gps_genomic = TRUE;
        }
      }
    }
  }

  PrintFtableIntervals (head, target, sfp->location, label);

  geneorprotdb = NULL;
  pseudo = sfp->pseudo;

  switch (context->seqfeattype) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        geneorprotdb = grp->db;
        pseudo |= grp->pseudo;

        StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
        if (! StringHasNoText (str)) {
          sprintf (tmp, "\t\t\tgene\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tgene_syn\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
        if (! StringHasNoText (grp->desc)) {
          sprintf (tmp, "\t\t\tgene_desc\t%s\n", grp->desc);
          ValNodeCopyStr (head, 0, tmp);
        }
        if (! StringHasNoText (grp->maploc)) {
          sprintf (tmp, "\t\t\tmap\t%s\n", grp->maploc);
          ValNodeCopyStr (head, 0, tmp);
        }
        if (! StringHasNoText (grp->locus_tag)) {
          sprintf (tmp, "\t\t\tlocus_tag\t%s\n", grp->locus_tag);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
      break;
    case SEQFEAT_CDREGION :
      prod = BioseqFind (SeqLocId (sfp->product));
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      if (prot != NULL) {
        prp = (ProtRefPtr) prot->data.value.ptrvalue;
        if (prp != NULL) {
          geneorprotdb = prp->db;
          if (prp->name != NULL) {
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
          }
          if (prp->desc != NULL) {
            StringNCpy_0 (str, prp->desc, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tprot_desc\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tfunction\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tEC_number\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
        StringNCpy_0 (str, prot->comment, sizeof (str));
        if (! StringHasNoText (str)) {
          sprintf (tmp, "\t\t\tprot_note\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        if (crp->frame > 1 && crp->frame <= 3) {
          sprintf (tmp, "\t\t\tcodon_start\t%d\n", (int) crp->frame);
          ValNodeCopyStr (head, 0, tmp);
        }
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          PrintFTCodeBreak (head, cbp, target);
        }
      }
      if (prod != NULL) {
        if (SeqIdWriteForTable (prod->id, str, sizeof (str), ajp, FALSE)) {
          sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
        if (is_gps_genomic) {
          cds = SeqMgrGetCDSgivenProduct (prod, NULL);
          if (cds != NULL) {
            cdna = BioseqFindFromSeqLoc (cds->location);
            if (cdna != NULL) {
              if (SeqIdWriteForTable (cdna->id, str, sizeof (str), ajp, FALSE)) {
                sprintf (tmp, "\t\t\ttranscript_id\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
          }
        }
      } else if (sfp->product != NULL) {
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          if (sip->choice == SEQID_GI) {
            sip2 = GetSeqIdForGI (sip->data.intvalue);
            if (sip2 != NULL) {
              sip = sip2;
            }
          }
          if (SeqIdWriteForTable (sip, str, sizeof (str), ajp, TRUE)) {
            sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      prod = BioseqFind (SeqLocId (sfp->product));
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL) {
        switch (rrp->ext.choice) {
          case 1 :
            StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tproduct\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
            break;
          case 2 :
            trp = rrp->ext.value.ptrvalue;
            if (trp != NULL) {
              FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
              numcodons = ComposeCodonsRecognizedString (trp, numbuf, sizeof (numbuf));
              if (numcodons > 0 && StringDoesHaveText (numbuf)) {
                sprintf (tmp, "\t\t\tcodon_recognized\t%s\n", numbuf);
                ValNodeCopyStr (head, 0, tmp);
              }
              slp = trp->anticodon;
              newloc = NULL;
              if (slp != NULL && ajp->ajp.slp != NULL) {
                sip = SeqIdParse ("lcl|dummy");
                newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle);
                SeqIdFree (sip);
                slp = newloc;
                if (newloc != NULL) {
                  A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
                }
              }
              aa = str;
              if (StringNICmp (aa, "tRNA-", 5) == 0) {
                aa += 5;
              }
              if (slp != NULL && StringDoesHaveText (aa)) {
                tmpx = FFFlatLoc (ajp, target, slp, ajp->masterStyle);
                if (tmpx != NULL) {
                  sprintf (tmp, "\t\t\tanticodon\t(pos:%s,aa:%s)\n", tmpx, aa);
                  ValNodeCopyStr (head, 0, tmp);
                }
                MemFree (tmpx);
              }
              if (newloc != NULL) {
                SeqLocFree (newloc);
              }
            }
            break;
          default :
            break;
        }
      }
      if (prod != NULL) {
        if (SeqIdWriteForTable (prod->id, str, sizeof (str), ajp, FALSE)) {
          sprintf (tmp, "\t\t\ttranscript_id\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
        if (is_gps_genomic) {
          cds = SeqMgrGetNextFeature (prod, NULL, SEQFEAT_CDREGION, 0, &fcontext);
          if (cds != NULL && SeqMgrGetNextFeature (prod, cds, SEQFEAT_CDREGION, 0, &fcontext) == NULL) {
            prod = BioseqFindFromSeqLoc (cds->product);
            if (prod != NULL) {
              if (SeqIdWriteForTable (prod->id, str, sizeof (str), ajp, FALSE)) {
                sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
          }
        }
      } else if (sfp->product != NULL) {
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          if (sip->choice == SEQID_GI) {
            sip2 = GetSeqIdForGI (sip->data.intvalue);
            if (sip2 != NULL) {
              sip = sip2;
            }
          }
          if (SeqIdWriteForTable (sip, str, sizeof (str), ajp, TRUE)) {
            sprintf (tmp, "\t\t\ttranscript_id\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp != NULL) {
        if (prp->name != NULL) {
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tproduct\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
        if (prp->desc != NULL) {
          StringNCpy_0 (str, prp->desc, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tprot_desc\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
        for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tfunction\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
        for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tEC_number\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
      }
      StringNCpy_0 (str, sfp->comment, sizeof (str));
      if (! StringHasNoText (str)) {
        sprintf (tmp, "\t\t\tprot_note\t%s\n", str);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_REGION :
      StringNCpy_0 (str, (CharPtr) sfp->data.value.ptrvalue, sizeof (str));
      if (! StringHasNoText (str)) {
        sprintf (tmp, "\t\t\tregion\t%s\n", str);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_BOND :
      bondidx = (Int2) sfp->data.value.intvalue;
      if (bondidx == 255) {
        bondidx = 5;
      }
      if (bondidx > 0 && bondidx < 6) {
        sprintf (tmp, "\t\t\tbond_type\t%s\n", bondList [bondidx]);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_SITE :
      siteidx = (Int2) sfp->data.value.intvalue;
      if (siteidx == 255) {
        siteidx = 26;
      }
      if (siteidx > 0 && siteidx < 27) {
        sprintf (tmp, "\t\t\tsite_type\t%s\n", siteList [siteidx]);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_PSEC_STR :
      sec_str = (Int2) sfp->data.value.intvalue;
      if (sec_str > 0 && sec_str <= 3) {
        sprintf (tmp, "\t\t\tsec_str_type\t%s\n", secStrText [sec_str]);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_HET :
      StringNCpy_0 (str, (CharPtr) sfp->data.value.ptrvalue, sizeof (str));
      if (! StringHasNoText (str)) {
        sprintf (tmp, "\t\t\theterogen\t%s\n", str);
        ValNodeCopyStr (head, 0, tmp);
      }
      break;
    case SEQFEAT_BIOSRC :
      PrintBioSourceFtableEntry (head, sfp->data.value.ptrvalue);
      break;
    default :
      break;
  }
  if (pseudo) {
    ValNodeCopyStr (head, 0, "\t\t\tpseudo\n");
  }
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (SeqMgrGeneIsSuppressed (grp)) {
      ValNodeCopyStr (head, 0, "\t\t\tgene\t-\n");
    } else {
      if (StringDoesHaveText (grp->locus)) {
        sprintf (tmp, "\t\t\tgene\t%s\n", grp->locus);
        ValNodeCopyStr (head, 0, tmp);
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        sprintf (tmp, "\t\t\tlocus_tag\t%s\n", grp->locus_tag);
        ValNodeCopyStr (head, 0, tmp);
      }
    }
  }
  if (! StringHasNoText (sfp->comment)) {
    ValNodeCopyStr (head, 0, "\t\t\tnote\t");
    ValNodeCopyStr (head, 0, sfp->comment);
    ValNodeCopyStr (head, 0, "\n");
  }
  switch (sfp->exp_ev) {
    case 1 :
      ValNodeCopyStr (head, 0, "\t\t\tevidence\texperimental\n");
      break;
    case 2 :
      ValNodeCopyStr (head, 0, "\t\t\tevidence\tnot_experimental\n");
      break;
    default :
      break;
  }
  if (! StringHasNoText (sfp->except_text)) {
    ValNodeCopyStr (head, 0, "\t\t\texception\t");
    ValNodeCopyStr (head, 0, sfp->except_text);
    ValNodeCopyStr (head, 0, "\n");
  } else if (sfp->excpt) {
    ValNodeCopyStr (head, 0, "\t\t\texception\n");
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (! StringHasNoText (gbq->qual)) {
      if (! StringHasNoText (gbq->val)) {
        sprintf (tmp, "\t\t\t%s\t%s\n", gbq->qual, gbq->val);
        ValNodeCopyStr (head, 0, tmp);
      }
    }
  }
  VisitUserObjectsInUop (sfp->ext, (Pointer) head, PrintFTUserObj);
  for (vnp = geneorprotdb; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
          ValNodeCopyStr (head, 0, tmp);
        } else {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
    }
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
          ValNodeCopyStr (head, 0, tmp);
        } else {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
    }
  }
}

static BioseqPtr FindFirstBioseq (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sep == NULL || sep->data.ptrvalue == NULL ||
      /* sep->choice < 0 || */ sep->choice > 2) return NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    return bsp;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    bsp = FindFirstBioseq (sep);
    if (bsp != NULL) return bsp;
  }
  return NULL;
}

static BioseqPtr BioseqLockAndIndexByEntity (Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;

  if (entityID < 1) return NULL;

  sep = SeqMgrGetSeqEntryForEntityID (entityID);
  if (sep == NULL) return NULL;

  bsp = FindFirstBioseq (sep);
  if (bsp == NULL) return NULL;

  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return NULL;

  bsp = BioseqLockById (sip);
  if (bsp == NULL) return NULL;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  return bsp;
}

NLM_EXTERN CharPtr FormatFtableSourceFeatBlock (
  BaseBlockPtr bbp,
  BioseqPtr target
)

{
  SeqFeatPtr        sfp;
  SeqDescPtr        sdp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioSourcePtr      biop;
  ValNodePtr        head;
  IntSrcBlockPtr    isp;
  CharPtr           str;

  if (bbp == NULL) return NULL;

  isp = (IntSrcBlockPtr) bbp;
  head = NULL;
  biop = NULL;

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID,
                                      0, NULL, &dcontext);
    if (sdp == NULL) return NULL;
    biop = sdp->data.ptrvalue;
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp == NULL) return NULL;
    biop = sfp->data.value.ptrvalue;
  }
  if (biop == NULL) return NULL;
  PrintFtableIntervals (&head, target, isp->loc, "source");
  PrintBioSourceFtableEntry (&head, biop);

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return str;
}

NLM_EXTERN void DoImmediateFormat (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  BlockType        blocktype;
  BioseqPtr        bsp;
  FormatProc       fmt;
  Boolean          is_www;
  size_t           max;
  SeqEntryPtr      oldscope;
  QualValPtr       qv = NULL;
  SeqEntryPtr      sep;
  CharPtr          str = NULL;

  if (afp == NULL || bbp == NULL) return;
  ajp = afp->ajp;
  if (ajp == NULL) return;
  is_www = GetWWW (ajp);

  blocktype = bbp->blocktype;
  if (blocktype < LOCUS_BLOCK || blocktype > SLASH_BLOCK) return;
  fmt = asn2gnbk_fmt_functions [(int) blocktype];
  if (fmt == NULL) return;

  max = (size_t) (MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR));
  qv = MemNew (sizeof (QualVal) * (max + 5));
  if (qv == NULL) return;

  sep = GetTopSeqEntryForEntityID (bbp->entityID);

  bsp = BioseqLockAndIndexByEntity (bbp->entityID);
  oldscope = SeqEntrySetScope (sep);

  afp->qvp = qv;
  str = fmt (afp, bbp);
  afp->qvp = NULL;

  SeqEntrySetScope (oldscope);
  BioseqUnlock (bsp);

  if (str != NULL) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "%s", str);
    }
    if (afp->ffwrite != NULL) {
      afp->ffwrite (str, afp->userdata, blocktype);
    }
  } else {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "?\n");
    }
    if (afp->ffwrite != NULL) {
      afp->ffwrite ("?\n", afp->userdata, blocktype);
    }
  }

  MemFree (str);
  MemFree (qv);
}

NLM_EXTERN void DoQuickLinkFormat (
  Asn2gbFormatPtr afp,
  CharPtr str
)

{
  if (afp == NULL || StringHasNoText (str)) return;

  if (str != NULL) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "%s", str);
    }
    if (afp->ffwrite != NULL) {
      afp->ffwrite (str, afp->userdata, (BlockType) 0);
    }
  }
}

NLM_EXTERN CharPtr asn2gnbk_format (
  Asn2gbJobPtr ajp,
  Int4 paragraph
)

{
  Asn2gbFormat     af;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp;
  BlockType        blocktype;
  BioseqPtr        bsp;
  FormatProc       fmt;
  IntAsn2gbJobPtr  iajp;
  size_t           max;
  SeqEntryPtr      oldscope;
  QualValPtr       qv;
  Int4             section;
  SeqEntryPtr      sep;
  CharPtr          str = NULL;

  /* qv must hold MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR) */

  iajp = (IntAsn2gbJobPtr) ajp;
  if (iajp == NULL || ajp->sectionArray == NULL || ajp->paragraphArray == NULL) return NULL;
  if (paragraph < 0 || paragraph >= ajp->numParagraphs) return NULL;

  bbp = ajp->paragraphArray [paragraph];
  if (bbp == NULL) return NULL;

  section = bbp->section;
  if (section < 0 || section >= ajp->numSections) return NULL;

  asp = ajp->sectionArray [section];
  if (asp == NULL) return NULL;

  blocktype = bbp->blocktype;
  if (blocktype < LOCUS_BLOCK || blocktype > SLASH_BLOCK) return NULL;

  max = (size_t) (MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR));
  qv = MemNew (sizeof (QualVal) * (max + 5));
  if (qv == NULL) return NULL;

  MemSet ((Pointer) &af, 0, sizeof (Asn2gbFormat));
  af.ajp = (IntAsn2gbJobPtr) ajp;
  af.asp = asp;
  af.qvp = qv;
  af.format = iajp->format;
  af.aip = iajp->aip;
  af.atp = iajp->atp;

  sep = GetTopSeqEntryForEntityID (bbp->entityID);

  fmt = asn2gnbk_fmt_functions [(int) blocktype];
  if (fmt == NULL) return NULL;

  bsp = BioseqLockAndIndexByEntity (bbp->entityID);
  oldscope = SeqEntrySetScope (sep);

  str = fmt (&af, bbp);

  SeqEntrySetScope (oldscope);
  BioseqUnlock (bsp);

  if (str == NULL) {
    str = StringSave ("???\n");
  }

  MemFree (qv);

  return str;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_cleanup (
  Asn2gbJobPtr ajp
)

{
  Asn2gbSectPtr     asp;
  BaseBlockPtr      bbp;
  BaseBlockPtr      PNTR blockArray;
  BioseqPtr         bsp;
  ValNodePtr        gapvnp;
  Int4              i;
  IntAsn2gbJobPtr   iajp;
  IntAsn2gbSectPtr  iasp;
  IntCdsBlockPtr    icp;
  IntFeatBlockPtr   ifp;
  IntRefBlockPtr    irp;
  IntSrcBlockPtr    isp;
  Int4              j;
  Int4              numBlocks;
  Int4              numSections;
  RefBlockPtr       rbp;
  Asn2gbFreeFunc    remotefree;
  ValNodePtr        remotevnp;
  SeqAnnotPtr       sap;
  SeqAnnotPtr       sapnext;
  Asn2gbSectPtr     PNTR sectionArray;
  StringItemPtr     sip, nxt;
  SeqBlockPtr       sbp;
  ValNodePtr        vnp;

  iajp = (IntAsn2gbJobPtr) ajp;
  if (iajp == NULL) return NULL;

  SeqLocFree (iajp->ajp.slp);

  numSections = ajp->numSections;
  sectionArray = ajp->sectionArray;

  if (sectionArray != NULL) {

    for (i = 0; i < numSections; i++) {
      asp = sectionArray [i];
      if (asp != NULL) {
        iasp = (IntAsn2gbSectPtr) asp;

        numBlocks = asp->numBlocks;
        blockArray = asp->blockArray;
        if (blockArray != NULL) {

          for (j = 0; j < numBlocks; j++) {
            bbp = blockArray [j];
            if (bbp != NULL) {

              MemFree (bbp->string);

              if (bbp->blocktype == REFERENCE_BLOCK) {
                rbp = (RefBlockPtr) bbp;
                MemFree (rbp->uniquestr);
                irp = (IntRefBlockPtr) rbp;
                DateFree (irp->date);
                SeqLocFree (irp->loc);
                MemFree (irp->authstr);
                MemFree (irp->fig);
                MemFree (irp->maploc);

              } else if (bbp->blocktype == SOURCEFEAT_BLOCK) {

                isp = (IntSrcBlockPtr) bbp;
                SeqLocFree (isp->loc);

              } else if (bbp->blocktype == FEATURE_BLOCK) {

                ifp = (IntFeatBlockPtr) bbp;
                if (ifp->isCDS) {
                  icp = (IntCdsBlockPtr) ifp;
                  MemFree (icp->fig);
                  MemFree (icp->maploc);
                }

              } else if (bbp->blocktype == SEQUENCE_BLOCK) {

                sbp = (SeqBlockPtr) bbp;
                MemFree (sbp->bases);
              }

              MemFree (bbp);
            }
          }
        }
        MemFree (asp->blockArray);
        MemFree (asp->referenceArray);
        MemFree (asp);
      }
    }
  }

  MemFree (ajp->sectionArray);
  MemFree (ajp->paragraphArray);
  MemFree (ajp->paragraphByIDs);

  sip = iajp->pool;
  while (sip != NULL) {
    nxt = sip->next;
    MemFree (sip);
    sip = nxt;
  }

  if (iajp->lockedBspList != NULL) {
    UnlockFarComponents (iajp->lockedBspList);
  }

  if (iajp->gapvnp != NULL || iajp->remotevnp != NULL) {
    SeqMgrClearFeatureIndexes (ajp->entityID, NULL);
  }

  if (iajp->gapvnp != NULL) {
    gapvnp = iajp->gapvnp;
    gapvnp->next = NULL;  /* unlink in case remotevnp still linked after gapvnp */
    bsp = (BioseqPtr) gapvnp->data.ptrvalue;
    if (bsp != NULL) {
      sap = bsp->annot;
      while (sap != NULL) {
        sapnext = sap->next;
        SeqAnnotFree (sap);
        sap = sapnext;
      }
    }
    /* frees fake Bioseq that was created by MemNew, not BioseqNew */
    ValNodeFreeData (gapvnp);
  }

  if (iajp->remotevnp != NULL) {
    remotevnp = iajp->remotevnp;
    remotefree = iajp->remotefree;
    if (remotefree != NULL) {
      /* if remotefree exists, it is responsible for all freeing */
      remotefree (remotevnp, iajp->remotedata);
    } else {
      /* otherwise free Bioseqs and ValNode chain ourselves */
      for (vnp = remotevnp; vnp != NULL; vnp = vnp->next) {
        bsp = (BioseqPtr) vnp->data.ptrvalue;
        if (bsp != NULL) {
          BioseqFree (bsp);
        }
      }
      ValNodeFree (remotevnp);
    }
  }

  TextFsaFree (iajp->bad_html_fsa);

  ValNodeFree (iajp->gihead);

  free_buff ();
  FiniWWW (iajp);

  MemFree (iajp);

  return NULL;
}

NLM_EXTERN Boolean SeqEntryToGnbk (
  SeqEntryPtr sep,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  FILE *fp
)

{
  AsnIoPtr           aip = NULL;
  AsnIoPtr           aipfree = NULL;
  Asn2gbJobPtr       ajp;
  AsnTypePtr         atp = NULL;
  BioseqPtr          bsp = NULL;
  BioseqSetPtr       bssp = NULL;
  Boolean            do_gbseq_asn = FALSE;
  Boolean            do_gbseq_xml = FALSE;
  Asn2gbWriteFunc    ffwrite = NULL;
  GBSeqPtr           gbseq = NULL;
  GBSeq              gbsq;
  IntAsn2gbJobPtr    iajp;
  Boolean            rsult = FALSE;
  Int1               type = ASNIO_TEXT_OUT;
  Pointer            userdata = NULL;
  XtraBlock          xtra;
  /*
  BaseBlockPtr       bbp;
  BlockType          block;
  CharPtr            ffhead = NULL;
  CharPtr            fftail = NULL;
  Int4               i;
  Boolean            is_html;
  Int4               numParagraphs;
  BaseBlockPtr PNTR  paragraphArray;
  CharPtr            str;
  */
#ifdef WIN_MAC
#if __profile__
  ValNodePtr         bsplist = NULL;
  Uint2              entityID;
  Boolean            lockFarComp;
  Boolean            lockFarLocs;
  Boolean            lockFarProd;
  Boolean            lookupFarComp;
  Boolean            lookupFarHist;
  Boolean            lookupFarInf;
  Boolean            lookupFarLocs;
  Boolean            lookupFarOthers;
  Boolean            lookupFarProd;
#endif
#endif

  if (extra != NULL) {
    ffwrite = extra->ffwrite;
    /*
    ffhead = extra->ffhead;
    fftail = extra->fftail;
    */
    gbseq = extra->gbseq;
    aip = extra->aip;
    atp = extra->atp;
    userdata = extra->userdata;
  }
  if (fp == NULL && ffwrite == NULL && aip == NULL) return FALSE;
  if (sep == NULL && slp == NULL) return FALSE;
  if (sep != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
    }
  }

#ifdef WIN_MAC
#if __profile__
  /* this allows profiling of just the formatter, without feature indexing, on the Mac */

  if (sep != NULL) {
    entityID = ObjMgrGetEntityIDForPointer (sep->data.ptrvalue);
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }
  }

  lockFarComp = (Boolean) ((locks & LOCK_FAR_COMPONENTS) != 0);
  lockFarLocs = (Boolean) ((locks & LOCK_FAR_LOCATIONS) != 0);
  lockFarProd = (Boolean) ((locks & LOCK_FAR_PRODUCTS) != 0);

  if (lockFarComp || lockFarLocs || lockFarProd) {
    locks = locks ^ (LOCK_FAR_COMPONENTS | LOCK_FAR_LOCATIONS | LOCK_FAR_PRODUCTS);
    if (slp != NULL && lockFarComp) {
      bsplist = LockFarComponentsEx (sep, FALSE, lockFarLocs, lockFarProd, slp);
    } else {
      bsplist = LockFarComponentsEx (sep, lockFarComp, lockFarLocs, lockFarProd, NULL);
    }
  }

  lookupFarComp = (Boolean) ((locks & LOOKUP_FAR_COMPONENTS) != 0);
  lookupFarLocs = (Boolean) ((locks & LOOKUP_FAR_LOCATIONS) != 0);
  lookupFarProd = (Boolean) ((locks & LOOKUP_FAR_PRODUCTS) != 0);
  lookupFarHist = (Boolean) ((locks & LOOKUP_FAR_HISTORY) != 0);
  lookupFarInf = (Boolean) ((locks & LOOKUP_FAR_INFERENCE) != 0);
  lookupFarOthers = (Boolean) ((locks & LOOKUP_FAR_OTHERS) != 0);

  if (lookupFarComp || lookupFarLocs || lookupFarProd || lookupFarHist || lookupFarInf || lookupFarOthers) {
    locks = locks ^ (LOOKUP_FAR_COMPONENTS | LOOKUP_FAR_LOCATIONS | LOOKUP_FAR_PRODUCTS | LOOKUP_FAR_HISTORY | LOOKUP_FAR_INFERENCE | LOOKUP_FAR_OTHERS);
    LookupFarSeqIDs (sep, lookupFarComp, lookupFarLocs, lookupFarProd, FALSE, lookupFarHist, lookupFarInf, lookupFarOthers);
  }

  ProfilerSetStatus (TRUE);
#endif
#endif

  do_gbseq_xml = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_XML_GBSEQ_FILE);
  do_gbseq_asn = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_ASN_GBSEQ_FILE);

  if (do_gbseq_xml || do_gbseq_asn) {
    if (fp != NULL && aip == NULL) {
      if (do_gbseq_xml) {
        type |= ASNIO_XML;
      }
      aip = AsnIoNew (type, fp, NULL, NULL, NULL);
      aipfree = aip;
      fp = NULL;
    }
    if (extra == NULL) {
      MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
      extra = &xtra;
    }
    if (extra->gbseq == NULL) {
      MemSet ((Pointer) &gbsq, 0, sizeof (GBSeq));
      extra->gbseq = &gbsq;
      gbseq = extra->gbseq;
    }
  }

  /* pass TRUE for stream to do immediate write at time of creation for speed */

  ajp = asn2gnbk_setup_ex (bsp, bssp, slp, format, mode, style,
                           flags, locks, custom, extra,
                           TRUE, fp, aip, atp);

  if (ajp != NULL) {
    rsult = TRUE;
    iajp = (IntAsn2gbJobPtr) ajp;

#if 0
    /* if streaming, all output was written in setup function, otherwise output here */

    if (! stream) {

      /* send optional head string */

      is_html = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_HTML_FLATFILE);
      if (ffhead == NULL && is_html) {
        ffhead = defHead;
      }
      if (ffhead != NULL) {
        if (fp != NULL) {
          fprintf (fp, ffhead);
        }
      }
      if (ffwrite != NULL) {
        ffwrite (ffhead, userdata, HEAD_BLOCK);
      }

      /* send each paragraph */

      numParagraphs = ajp->numParagraphs;
      paragraphArray = ajp->paragraphArray;

      for (i = 0; i < numParagraphs; i++) {
        str = asn2gnbk_format (ajp, i);
        block = (BlockType) 0;
        if (paragraphArray != NULL) {
          bbp = paragraphArray [i];
          if (bbp != NULL) {
            block = bbp->blocktype;
          }
        }
        if (str != NULL) {
          if (fp != NULL) {
            fprintf (fp, "%s", str);
          }
          if (ffwrite != NULL) {
            ffwrite (str, userdata, block);
          }
        } else {
          if (fp != NULL) {
            fprintf (fp, "?\n");
          }
          if (ffwrite != NULL) {
            ffwrite ("?\n", userdata, block);
          }
        }

        MemFree (str);
      }

      /* send optional tail string */

      if (fftail == NULL && is_html) {
        fftail = defTail;
      }
      if (fftail != NULL) {
        if (fp != NULL) {
          fprintf (fp, fftail);
        }
      }
      if (ffwrite != NULL) {
        ffwrite (fftail, userdata, TAIL_BLOCK);
      }
    }
#endif

    /* if RELEASE_MODE, warn if unresolved gi numbers, missing translation, etc. */
    
    if (iajp->relModeError && mode == RELEASE_MODE) {
      rsult = FALSE;
    }

    asn2gnbk_cleanup (ajp);
  }

  if (aipfree != NULL) {
    AsnIoFree (aipfree, FALSE);
  }

#ifdef WIN_MAC
#if __profile__
  ProfilerSetStatus (FALSE);

  UnlockFarComponents (bsplist);
#endif
#endif

  return rsult;
}

NLM_EXTERN Boolean BioseqToGnbk (
  BioseqPtr bsp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  FILE *fp
)

{
  SeqEntryPtr  sep = NULL;

  if (bsp == NULL && slp == NULL) return FALSE;
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  }
  return SeqEntryToGnbk (sep, slp, format, mode, style, flags, locks, custom, extra, fp);
}

