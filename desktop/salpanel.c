/*   salpanel.c
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
* File Name:  salpanel.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.81 $
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
#include <saledit.h>
#include <salpanel.h>
#include <salutil.h>
#include <salfiles.h>
#include <salstruc.h>
#include <fstyle.h>
#include <salmedia.h>
#include <sequtil.h>
#include <alignmgr2.h>
#include <edutil.h>
#include <dlogutil.h>
#include <import.h>
#include <seqpanel.h>

#define OBJ_VIRT 254

static Uint1 rectSym [] = {
  0xFE, 0x82, 0x82, 0x82, 0x82, 0x82, 0xFE, 0x00
};
static Uint1 diamondSym [] = {
  0x10, 0x28, 0x44, 0x82, 0x44, 0x28, 0x10, 0x00
};
static Uint1 ovalSym [] = {
  0x38, 0x44, 0x82, 0x82, 0x82, 0x44, 0x38, 0x00
};
static Uint1 leftTriSym [] = {
  0x06, 0x1A, 0x62, 0x82, 0x62, 0x1A, 0x06, 0x00
};
static Uint1 rightTriSym [] = {
  0xC0, 0xB0, 0x8C, 0x82, 0x8C, 0xB0, 0xC0, 0x00
};
static Uint1 upTriSym [] = {
  0x10, 0x28, 0x28, 0x44, 0x44, 0x82, 0xFE, 0x00
};
static Uint1 downTriSym [] = {
  0xFE, 0x82, 0x44, 0x44, 0x28, 0x28, 0x10, 0x00
};
static Uint1 rectFillSym [] = {
  0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0x00
};
static Uint1 diamondFillSym [] = {
  0x10, 0x38, 0x7C, 0xFE, 0x7C, 0x38, 0x10, 0x00
};
static Uint1 ovalFillSym [] = {
  0x38, 0x7C, 0xFE, 0xFE, 0xFE, 0x7C, 0x38, 0x00
};
static Uint1 leftTriFillSym [] = {
  0x06, 0x1E, 0x7E, 0xFE, 0x7E, 0x1E, 0x06, 0x00
};
static Uint1 rightTriFillSym [] = {
  0xC0, 0xF0, 0xFC, 0xFE, 0xFC, 0xF0, 0xC0, 0x00
};
static Uint1 upTriFillSym [] = {
  0x10, 0x38, 0x38, 0x7C, 0x7C, 0xFE, 0xFE, 0x00
};
static Uint1 downTriFillSym [] = {
  0xFE, 0xFE, 0x7C, 0x7C, 0x38, 0x38, 0x10, 0x00
};
static Uint1 rightOvalSym [] = {
  0x18, 0x14, 0x12, 0x12, 0x12, 0x14, 0x18, 0x00
};
static Uint1 leftOvalSym [] = {
  0x38, 0x44, 0x82, 0x82, 0x82, 0x44, 0x38, 0x00
};
static Uint1 rightOvalFillSym [] = {
  0x18, 0x1C, 0x1E, 0x1E, 0x1E, 0x1C, 0x18, 0x00
};
static Uint1 leftOvalFillSym [] = {
  0x30, 0x50, 0x90, 0x90, 0x90, 0x50, 0x30, 0x00
};

/*######################################################################
#
#       functions for setting up the color for different object
#
######################################################################*/
#define RGB_B(x) (Uint1)((x)&255);
#define RGB_G(x) (Uint1)(((x)>>8)&255);
#define RGB_R(x) (Uint1)(((x)>>16)&255);

static Boolean convert_color(Int4 val, Uint1Ptr color)
{

        if(val<0 || color == NULL)
                return FALSE;
        color[0] = RGB_R(val);
        color[1] = RGB_G(val);
        color[2] = RGB_B(val);
        return TRUE;
}

static Uint4 getcolor_fromstyles (Uint2 itemsubtype)
{
  Int4 c_val;
  Uint1 color [3];

  c_val = GetMuskCParam (itemsubtype, MSM_SEGMENT, MSM_COLOR);
  if ( convert_color (c_val, color) )
     return GetColorRGB (color[0], color[1], color[2]);
  return GetColorRGB (0, 0, 0);
}

static Uint4 getcolorforvirtfeat (Uint2 item)
{
  return GetColorRGB ((Uint1)((item % 3) * 40), 
                      (Uint1)(255-(item % 3) * 40), 
                      (Uint1)((item % 3) * 40)); 
}

/**********************************************
***   GetAlignEditData
**********************************************/
extern WindoW getwindow_frompanel (PaneL pnl)
{
/*
#ifdef WIN_MAC
  w = FrontWindow ();
#else
  w = ParentWindow (pnl);
#endif
*/
  return ParentWindow (pnl);
}

extern PaneL GetPanelFromWindow (WindoW w)
{
  SeqEditViewFormPtr wdp;

  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  if (wdp == NULL) { 
    return NULL; 
  }
  if ( wdp->pnl == NULL ) { 
    return NULL;
  }
  return wdp->pnl;
}

extern EditAlignDataPtr GetAlignEditData (WindoW w)
{
  PaneL            pnl;
  EditAlignDataPtr adp;

  pnl = GetPanelFromWindow (w);
  if (pnl == NULL) return NULL;
  GetPanelExtra (pnl, &adp);
  if ( adp == NULL )
     return NULL;
  if ( adp->firstssp == NULL )
     return NULL;
  return adp;
}

extern EditAlignDataPtr GetAlignDataPanel (PaneL pnl)
{
  EditAlignDataPtr adp;

  if ( pnl == NULL )
     return NULL;
  GetPanelExtra (pnl, &adp);
  if ( adp == NULL )
     return NULL;
  if ( adp->firstssp == NULL )
     return NULL;
  return adp;
}

/**********************************************
***   AlignDataSet_Restore
***     Inval Rect
***     RestorePort
**********************************************/
static void get_client_rect (PaneL p, RectPtr prc)
{
  ObjectRect (p, prc);
  InsetRect (prc, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
}

static void get_client_rectxy (PaneL p, RectPtr prc, Int2 x, Int2 y)
{
  ObjectRect (p, prc);
  InsetRect (prc, x, y);
}

/**************************************************************
***
***   
**************************************************************/
static SelStructPtr get_firstline (SelEdStructPtr sesp1, SelStructPtr buffer)
{
  SelEdStructPtr sesp;
  SelStructPtr   buf;

  if (buffer == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in get_firstline [25]");
         return NULL;
  }
  if (sesp1 == NULL) {
         return buffer;
  }
  for (buf=buffer; buf!=NULL; buf=buf->next) 
  {
         sesp = (SelEdStructPtr) buf->region;
         if (is_sameId (sesp1->entityID, sesp1->itemID, sesp1->itemtype, 255, sesp->entityID, sesp->itemID, sesp->itemtype, 255) ) 
            break;
  }
  if (buf == NULL) {
         return buffer;
  }
  return buf;
}

/*********************************************************************
***  FindNextSegment
*********************************************************************/
static SelEdStructPtr FindNextSegment (SelEdStructPtr current)
{
  Uint2  ei, ii;
  if ( current->next == NULL ) return NULL;
  ei = current->entityID;
  ii = current->itemID;
  current = current->next;
  while ( current != NULL )
  {
     if (current->entityID == ei &&  current->itemID == ii) break;
     current = current->next;
  }
  return current;
}

static CharPtr get_substring (CharPtr str, Int4 drw_start, Int4 drw_width)
{
  Int4            width;
  Int4            stringlens;
  CharPtr         strp;

  if (str == NULL ) 
     return NULL; 
  if (str[0]=='\0')
     return NULL;
  stringlens = StringLen (str);
  if ( drw_start >= stringlens ) { 
     return NULL; } 
  strp = str + drw_start;
  stringlens = StringLen (strp);
  if (stringlens == 0) 
     return NULL; 
  width = MIN ((Int4) drw_width, (Int4) stringlens);
  if ( !not_empty_string (strp, width) ) 
     return NULL;
  return strp;
}

static SelStructPtr go_to_next_to_draw (EditAlignDataPtr adp, Boolean next, Int4 offset)
{
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  Int4             from_inbuf;  /* alignment coordinates in buffer */
  Int4             from_inseq;  /* alignment coordinates in buffer */
  Int4             drw_width;   /* length of drw_str */
  Uint2            itemsubtype;
  CharPtr          curstr = NULL;
  Int4             offsettmp;

  if (adp->voffset == 0)
     return adp->buffer;
  if (offset == 0)
     return adp->firstssp;
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int4) adp->visibleWidth, (Int4) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 ) {
         return NULL;
  }
  curvnp = adp->buffer;
  offsettmp = 0;
  while (offsettmp < offset && from_inseq < adp->length)
  {
         curssp = (SelEdStructPtr) curvnp->region;
         itemsubtype = curvnp->itemtype;
         if (itemsubtype == EDITDEF_SCA)      
                offsettmp++;
         else if (itemsubtype == EDITDEF_SCB) 
                offsettmp;
         else if (curssp->itemtype ==OBJ_BIOSEQ) /*&&itemsubtype==FEATDEF_BAD)*/

         {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )   
                   offsettmp++;
         }
         if (offsettmp < offset) {
            if (next)
               curvnp = curvnp->next;
            else
               curvnp = curvnp->prev;
         }
         if (offsettmp < offset && curvnp == NULL) 
         {
            if (next) {
               curvnp = adp->buffer;
               adp->hoffset += adp->visibleWidth;
               from_inbuf += drw_width;
               from_inseq += drw_width;
            } else {
               curvnp = adp->buffertail;
               adp->hoffset -= adp->visibleWidth;
               from_inbuf -= drw_width;
               from_inseq -= drw_width;
            } 
            if (from_inseq >= adp->length) break;
            if (from_inseq <= 0 || from_inbuf <= 0 || adp->hoffset<=0) break;
            drw_width = MIN ((Int4) adp->visibleWidth, 
                (Int4)(adp->bufferlength -(adp->hoffset -adp->bufferstart)));
         }
  }
  if (from_inseq < 0 || from_inbuf < 0 || adp->hoffset<0)  {
     adp->hoffset = adp->bufferstart;
     return adp->buffer;
  }
  if (from_inseq >= adp->length) 
     return NULL;
  return curvnp;
}


static SelStructPtr next_to_draw (EditAlignDataPtr adp, Boolean next)
{
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelStructPtr     ssptmp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  SeqLocPtr        curslp;
  SeqIdPtr         sip;
  CharPtr          curstr = NULL;
  Int4             start, stop;
  Int4             start2, stop2;
  Int4             from_inbuf;  /* alignment coordinates in buffer */
  Int4             from_inseq;  /* alignment coordinates in buffer */
  Int4             drw_width;   /* length of drw_str */
  Int4             chklocp;
  Uint2            itemsubtype;
  SeqAlignPtr      salp = (SeqAlignPtr) adp->sap_align->data;
  Boolean          empty_line;

  if (next) {
     if ( adp->firstssp->next == NULL)
     {
         adp->hoffset += adp->visibleWidth;
         return adp->buffer;
     }
     else {
         ssptmp = adp->firstssp->next;
     }
  }
  if (!next) {
     if ( adp->firstssp->prev == NULL)
     {
         if (adp->hoffset == 0) {
            adp->voffset = 0;
            return adp->buffer;
         }
         adp->hoffset -= adp->visibleWidth;
         ssptmp = adp->buffertail;
     }
     else ssptmp = adp->firstssp->prev;
  }
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int4) adp->visibleWidth, 
                   (Int4) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 ) {
         return NULL;
  }
  empty_line = FALSE;
  curvnp = ssptmp;
  while ( from_inseq < adp->length )
  {
         curssp = (SelEdStructPtr) curvnp->region;
         itemsubtype = curvnp->itemtype;
         if (itemsubtype == EDITDEF_SCA)      
                return curvnp;
         if (itemsubtype == EDITDEF_SCB) 
                return curvnp;
         if (itemsubtype == FEATDEF_TRSL) 
         {
            if (!empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<=stop && start<from_inseq+drw_width && stop>from_inseq) 
                     return curvnp;      
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else if (itemsubtype>=EDITDEF_RF1 && itemsubtype<=EDITDEF_RF6) 
         {
            if (curssp->region !=NULL && !empty_line) {
                curslp = (SeqLocPtr) curssp->region;
                if ( SeqLocStop (curslp) > from_inseq && curssp->data != NULL) 
                   return curvnp; 
            }
         }
         else if (curssp->itemtype ==OBJ_BIOSEQ) /*&&itemsubtype==FEATDEF_BAD)*/

         {
            if (curssp->data !=NULL) {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )   
                   if (adp->draw_emptyline || (!adp->draw_emptyline && !stringhasnochar(curstr, 0, drw_width) )) {
                      empty_line = FALSE;
                      return curvnp;
                   } 
                   else 
                      empty_line = TRUE;
            }
         }
         else if (curssp->itemtype==OBJ_SEQFEAT && itemsubtype==FEATDEF_PROT)
         {
            if (curssp->data !=NULL && !empty_line) {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                   
                   return curvnp;      
            }
         }
         else if (curssp->itemtype ==OBJ_SEQFEAT)
         {
            if ( !empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) { 
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<=stop && start<from_inseq+drw_width && stop>=from_inseq) 
                     return curvnp;      
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else if (itemsubtype == EDITDEF_CPL) 
         {
            if (curssp->data !=NULL && !empty_line) { 
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring (curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                        
                   return curvnp;      
            }
         }

         else if( itemsubtype==SEQFEAT_GENE || itemsubtype==SEQFEAT_RNA
              ||  itemsubtype==SEQFEAT_CDREGION)
         {
            if ( !empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<=stop && start<from_inseq+drw_width && stop>from_inseq) 
                     return curvnp;
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else {
/**
            ErrPostEx (SEV_ERROR, 0, 0, "fail in next_to_draw [46] subtype %ld type %ld", (long) itemsubtype, (long) curssp->itemtype); 
**/
            return NULL;
         }
         if (next)
            curvnp = curvnp->next;
         else
            curvnp = curvnp->prev;
         if (curvnp == NULL) 
         {
            if (next) {
               curvnp = adp->buffer;
               adp->hoffset += adp->visibleWidth;
               from_inbuf += drw_width;
               from_inseq += drw_width;
            } else {
               curvnp = adp->buffertail;
               adp->hoffset -= adp->visibleWidth;
               from_inbuf -= drw_width;
               from_inseq -= drw_width;
            } 
            if (from_inseq >= adp->length) break;
            drw_width = MIN ((Int4) adp->visibleWidth, 
                (Int4)(adp->bufferlength -(adp->hoffset -adp->bufferstart)));
         }
  }
  return NULL;
}

/********************************************************
***  SetupScrollBar   
***
***  
*********************************************************/
extern BaR SeqEdGetSlateScrollBar (PaneL pnl)
{
  EditAlignDataPtr   adp;

  if (pnl)
  {
     adp=GetAlignDataPanel (pnl);
     if(adp!=NULL)
     {
        if (adp->vscrollbar_mode) {
           return GetSlateVScrollBar ((SlatE)pnl);      
        }
        return GetSlateHScrollBar ((SlatE)pnl); 
     }
  }
  return NULL;
}

extern Int4 SeqEdGetValueScrollBar (PaneL pnl)
{
  return GetBarValue(SeqEdGetSlateScrollBar(pnl));
}

extern void SeqEdSetValueScrollBar (PaneL pnl, Int4 value)
{
  SetBarValue(SeqEdGetSlateScrollBar(pnl), value);
}

extern void SeqEdCorrectBarPage (PaneL pnl, Int4 page1, Int4 page2)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarPage (sb, page1, page2);
}

extern void SeqEdCorrectBarValue (PaneL pnl, Int4 value)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarValue (sb, value);
}

extern void SeqEdCorrectBarMax (PaneL pnl, Int4 value)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarMax (sb, value);
}

extern void SeqEdSetCorrectBarMax (PaneL pnl, EditAlignDataPtr adp, float hratio)
{
  BaR          sb;
  Int4         cbm = 0;

  if (adp->nlines < 11) {
     adp->voffset = 0;
     adp->hoffset = 0;
     adp->firstssp = get_firstline (NULL, adp->buffer);
  }
  else {
     adp->hoffset = (Int4)(hratio * (float)adp->length); 
     adp->voffset = hoffset2voffset (adp, adp->anp_list, adp->visibleWidth, 0, adp->length-1, adp->hoffset);
  }

  if (adp->nlines < 0) 
     cbm = 0;
  else 
     cbm = MAX ((Int4) 0, (Int4) (adp->nlines -1));
  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarMax (sb, cbm);
  SetBarValue (sb, (Int4)(adp->voffset));
}

static void count_feature_buf_line (ValNodePtr fnp_list, Int4 g_left, Int4 g_right, ValNodePtr PNTR feature_line)
{
        FeatNodePtr fnp;
        Int4 c_left, c_right;
        ValNodePtr vnp;
        Boolean found;

        if(fnp_list == NULL)
                return;
        
        while(fnp_list)
        {
           fnp = (FeatNodePtr)fnp_list->data.ptrvalue;
           c_left = fnp->extremes.left;
           c_right = fnp->extremes.right;
           if(!(c_left > g_right || c_right < g_left))
           {
                found = FALSE;
                for(vnp = *feature_line; vnp != NULL && !found; vnp = vnp->next)
                {
                        if(vnp->data.intvalue == (Int4)(fnp->itemID))
                                found = TRUE;
                }
                if(!found)
                        ValNodeAddInt(feature_line, 0, (Int4)(fnp->itemID));
           }
           fnp_list = fnp_list->next;
        }
}

static Int4 CountTextAlignNodeNum(AlignNodePtr anp, Int4 m_left, Int4 m_right)
{
        Int4 num_line = 0;
        Int4 g_left, g_right;

        AlignSegPtr asp;
        ValNodePtr feature_line, curr;  /*the number of lines for a feature*/

        g_left = anp->extremes.left;
        g_right = anp->extremes.right;
        if(m_left > g_right || m_right < g_left)
                return 0;

        num_line = 1;
        feature_line = NULL;

        /*process  the GAPs and the DIAGs segs*/
        for(asp = anp->segs; asp !=NULL; asp = asp->next)
        {
           g_left = asp->gr.left;
           g_right = asp->gr.right;
           if(!(g_left > m_right || g_right < m_left))
           {
              switch(asp->type)
              {  
                case GAP_SEG:
                   break;

                case REG_SEG:
                case DIAG_SEG:
                   g_left = MAX(m_left, g_left);
                   g_right = MIN(m_right, g_right);
                   count_feature_buf_line (asp->cnp, g_left, g_right, &feature_line);
                   break;
                default:
                   break;
              }
           }
           if(g_left > m_right)
                break;
        }
        if(feature_line != NULL)
        {
           for(curr = feature_line; curr != NULL; curr = curr->next)
                ++num_line;
           ValNodeFree(feature_line);
        }
 
        return num_line;
}
static Int4 addline_perblock (EditAlignDataPtr adp, Int4 diffs)
{
  Int4        line = 0;
  ValNodePtr  vnp;
  SeqParamPtr prm;
  Int1        j;

  if (adp->draw_scale) line += (Int4) diffs;
  if (adp->draw_bars)  line += (Int4) diffs;
  for (vnp = adp->params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if ( prm->complement ) line += (Int4) diffs;
     for (j=0; j<=6; j++) 
        if (prm->rf[j]) line += (Int4) diffs;
  }
  return line;
}  
 
static Int4 feat_linenum (Int4 slp_start, Int4 slp_stop, Int4 line_len, Int4 left,
Int4 right)
{
  Int4         modstart;
  Int4         modstop;
 
  slp_start = MAX (slp_start, left);
  modstart = slp_start % line_len;
  if ( modstart > 0) slp_start -= modstart;
 
  slp_stop = MIN (slp_stop, right);
  modstop = slp_stop % line_len;
  if ( modstop > 0) slp_stop += line_len;
 
  return (Int4)((slp_stop - slp_start) / line_len);
}
 
static Int4 CountFeatNum (ValNodePtr adpfeat, Int4 line_len, Int4 left, Int4 right){
  ValNodePtr   vnp;
  SelEdStructPtr sesp;
  SeqLocPtr    slp;
  Int4         line = 0;
 
  for (vnp = adpfeat; vnp != NULL; vnp = vnp->next)
  {
     sesp = (SelEdStructPtr) vnp->data.ptrvalue;
     for (; sesp != NULL; sesp = sesp->next) 
     {
        if (vnp->choice ==FEATDEF_CDS && sesp->regiontype ==OM_REGION_SEQLOC 
        && sesp->region !=NULL)
        {
           slp = (SeqLocPtr) sesp->region;;
           if (SeqLocStart(slp) > right || SeqLocStop(slp) < left) 
              line+=0;
           else {
              line+= feat_linenum (SeqLocStart(slp), SeqLocStop(slp), line_len, left, right);
           }
        }
     }
  }
  return line;
}
static Int4 count_nline (EditAlignDataPtr adp, ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4 voffset)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        Int4 h_block = 0;
        ValNodePtr curr;

        if(anp_list == NULL)
                return h_block;
        if(voffset == 0)
                return h_block;
        anp = (AlignNodePtr)anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return h_block;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return h_block;
        c_start = left;
        while(line_num < voffset)
        {
                c_stop = MIN(right, (c_start+line_len-1));
                for(curr = anp_list; curr != NULL; curr = curr->next)
                {
                        anp = (AlignNodePtr)curr->data.ptrvalue;
                        line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
                }
                line_num += (Int4) addline_perblock (adp, 1);
                line_num += (Int4) CountFeatNum (adp->feat, line_len, c_start, c_stop);
                line_num += (Int4) CountFeatNum (adp->seqfeat, line_len, c_start, c_stop);
                if (line_num > voffset) break;
                ++h_block;
                c_start = c_stop+1;
        }
  return c_start;
}

extern Int4 hoffset2voffset (EditAlignDataPtr adp, ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4 hoffset)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        Int4 preline;
        Int4 h_block = 0;
        ValNodePtr curr;
 
        if(anp_list == NULL)
                return h_block;
        if(hoffset == 0)
                return h_block;
        anp = (AlignNodePtr)anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return h_block;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return h_block;
        c_start = left;
        while(c_start<hoffset)
        {
                preline = line_num;
                c_stop = MIN(right, (c_start+line_len-1));
                for(curr = anp_list; curr != NULL; curr = curr->next)
                {
                        anp = (AlignNodePtr)curr->data.ptrvalue;
                        line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
                }
                line_num += (Int4) addline_perblock (adp, 1);
                line_num += (Int4) CountFeatNum (adp->feat, line_len, c_start, c_stop);
                line_num += (Int4) CountFeatNum (adp->seqfeat, line_len, c_start, c_stop);
                ++h_block;
                c_start = c_stop+1;
        }

preline = line_num;
  return preline;
}


/*******************************************
***   Scrolling functions  
***
*** function test whether in/out buffer: 
***          not used 
***
********************************************/
extern void VscrlProc (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  EditAlignDataPtr adp;
  RecT         r;
  Int4         pixels;
  WindoW       tempPort;
  Int4         temp;
  Int4         x;
  Int4         oldhoffset;

  if ( s == NULL ) { 
    return; 
  }
  if ( (adp = GetAlignDataPanel ((PaneL) s)) == NULL ) 
     return;
  if ( adp->seqnumber == 0 ) 
     return;
  tempPort = SavePort ((PaneL) s);
  Select ((PaneL) s);
  ObjectRect ((PaneL) s, &r);
  if ((newval > oldval && newval - oldval <= adp->vPage ) 
   || (newval < oldval && oldval - newval <= adp->vPage )) 
  {
     InsetRect (&r, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
     pixels = (oldval - newval) * adp->lineheight;
     r.bottom = r.top + adp->pnlLine * adp->lineheight +1;
     r.top = r.top + 1;
     ScrollRect (&r, 0, pixels);
  } 
  adp->voffset = GetBarValue (sb);
  if (abs(newval - oldval) == 1)
  {
     oldhoffset = adp->hoffset;
     temp = oldhoffset  + adp->visibleWidth + adp->visibleLength;
     temp = MIN (temp, adp->bufferstart + adp->bufferlength);
     if ((oldhoffset + adp->visibleWidth > adp->bufferstart 
     && temp < adp->bufferstart + adp->bufferlength) 
     || temp == adp->length)
     {
         adp->firstssp = next_to_draw (adp, (Boolean)((newval-oldval)>0));
     }
     else {
         data_collect_arrange (adp, TRUE);
         if (temp < adp->length)
            InvalRect (&r);
     }
     if (adp->hoffset == 0 && adp->firstssp->prev == NULL && adp->voffset != 0) 
     {
         adp->voffset = 0;
         CorrectBarValue (sb, (Int4) 0);
     }
     if (adp->seqnumber == 1
         && adp->hoffset+adp->visibleWidth>adp->length
         && adp->voffset<adp->nlines)
     {
        adp->voffset = adp->nlines;
        CorrectBarValue (sb, (Int4) adp->nlines);
     }
  } 
  else {
  	 Int4 tmp;
     x = adp->seqnumber;
     if (adp->draw_scale) x++;
     if (adp->draw_bars) x++;
     
     if (adp->showfeat)
     {
	     adp->hoffset = count_nline (adp, adp->anp_list, adp->visibleWidth, 0, adp->length-1, (Int4)((FloatLo)adp->voffset));
     }
     else 
     {
         if (adp->seqnumber > 1)
     	     adp->hoffset = adp->visibleWidth * (Int4)((FloatLo)adp->voffset / (adp->seqnumber + 2));
         else
             adp->hoffset = adp->visibleWidth * (Int4)((FloatLo)adp->voffset); 
     }
     
     /* Use adp->int4value2 to disable line number counting when scrolling */
     tmp = adp->int4value2;
     adp->int4value2 = -1;
     if (adp->hoffset > adp->bufferstart 
     && adp->hoffset+adp->visibleLength+ adp->visibleWidth < adp->bufferstart +adp->bufferlength) 
     {
         data_collect_arrange (adp, FALSE);
     } else {
         data_collect_arrange (adp, TRUE);
     }
     adp->int4value2 = tmp;
     
     if (x == 0) {
        adp->firstssp=go_to_next_to_draw(adp, TRUE, (Int4)0);
     } else {
        adp->firstssp=go_to_next_to_draw(adp, TRUE, (Int4)(adp->voffset%x));
     }
     InvalRect (&r);
  }
  RestorePort(tempPort);
  Update ();
}

extern void HscrlProc (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  EditAlignDataPtr adp;
  RecT         r;
  Int4         pixels;
  WindoW       tempPort;

  if ( s == NULL ) {
    return;
  }
  if ( (adp = GetAlignDataPanel ((PaneL) s)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  tempPort = SavePort ((PaneL) s);
  Select ((PaneL) s);
  ObjectRect ((PaneL) s, &r);
  if ((newval > oldval && newval - oldval <= adp->hPage )
   || (newval < oldval && oldval - newval <= adp->hPage ))
  {
     InsetRect (&r, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
     pixels = (oldval - newval) * adp->charw +1;
     ScrollRect (&r, pixels, 0);
  } else {
     InvalRect (&r);
  }
  adp->hoffset = GetBarValue (sb);
  RestorePort(tempPort);
  Update ();
}

/*******************************************************************
***  
***    do_resize
***  
*******************************************************************/
extern void do_resize_panel (PaneL pnl, EditAlignDataPtr adp, Int4 width, Int4 height, Boolean rearrange)
{
  Int4         old_voffset;
  Int4         new_buffer;
  Int4         j;
  Int4         lg;
  Int4         old_visibleWidth;
  Int4         x, y;
  float hratio;
  
  x = (width - adp->margin.right) / adp->charw;
  if (x < 0) 
     x = 0;
  hratio = (float)adp->hoffset / (float)adp->length;
  adp->pnlWidth = x;  
  if (adp->pnlWidth < adp->marginleft + 10) {
     adp->firstssp = NULL;
     return;
  }
  x = (height - adp->margin.bottom) / adp->lineheight;
  if (x < 0) 
     x = 0;
  adp->pnlLine = x; 
  if (adp->pnlLine < 3) {
     adp->firstssp = NULL;
     return;
  }
  y = 0; x = 0;
  if (adp->columnpcell > 0) {
     y = (Int4) (adp->pnlWidth -adp->marginleft) / (Int4) adp->columnpcell;
     x = (Int4) (adp->pnlWidth -adp->marginleft -y) % (Int4)(adp->columnpcell);
     if (x == 9) 
        x = -1;
  }
  old_visibleWidth = adp->visibleWidth;
  adp->visibleWidth = (Int4) (adp->pnlWidth -adp->marginleft -y -x);
  if (adp->visibleWidth < 10) {
     adp->firstssp = NULL;
     return;
  }
  if ( adp->seqnumber == 0 ) 
     return;
  if (old_visibleWidth != adp->visibleWidth) {
     old_voffset = adp->voffset;
     adp->voffset = (Int4)(((float)old_visibleWidth/(float)adp->visibleWidth) * (float)old_voffset);
  }
  new_buffer = adp->pnlLine * adp->visibleWidth;
  if (new_buffer * 3 > adp->minbufferlength)
  {
     adp->minbufferlength = new_buffer * 3;
     if ( adp->colonne != NULL ) 
        MemFree (adp->colonne);
     adp->colonne = NULL;
     lg = adp->minbufferlength + adp->editbuffer + 4;
     adp->colonne = (Int4Ptr) MemNew ((size_t) (lg * sizeof(Int4)));
     for (j=0; j<adp->minbufferlength +adp->editbuffer; j++) adp->colonne[j] = -1;
     rearrange = TRUE;
  }
  adp->vPage = /* adp->pnlLine - 1; */ 1;
  adp->hPage = adp->visibleWidth - 1;
  data_collect_arrange (adp, rearrange);
  SeqEdSetCorrectBarMax (pnl, adp, hratio);
  SeqEdCorrectBarPage (pnl, adp->vPage, adp->vPage);
  SeqEdCorrectBarValue (pnl, SeqEdGetValueScrollBar (pnl));
}

extern void do_resize_window (PaneL pnl, EditAlignDataPtr adp, Boolean rearrange)
{
  SeqEditViewFormPtr wdp;
  WindoW       w;
  RecT         rw;    /* window rect        */
  RecT         rp;    /* panel rect         */
  RecT         rb;    /* buttons rect       */
  RecT         rct;   /* new rect for panel */
  Int4         buttonwidth,
               buttonheight;
  Int4         x, y;

  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  if (wdp == NULL) 
     return;
  ObjectRect (w, &rw);
  get_client_rect (pnl, &rp);
  GetPosition (pnl, &rp);
  x = adp->xoff;
  y = rw.bottom - rw.top;
  LoadRect ( &rct, x, adp->yoff, (Int4)(rw.right - rw.left - adp->x), 
                                 (Int4)(y - adp->y));
  SetPosition (pnl, &rct );
  AdjustPrnt (pnl, &rct, FALSE);
  ObjectRect (pnl, &rp );

  GetPosition ((GrouP) wdp->btngp, &rb);
  buttonwidth  = rb.right  - rb.left;
  buttonheight = rb.bottom - rb.top;
  LoadRect (&rb, x, (Int4)(y -buttonheight - adp->ybutt - 1),
                     (Int4)(x+ buttonwidth), (Int4)(y-adp->ybutt-1));
  SetPosition(wdp->btngp, &rb);
  AdjustPrnt (wdp->btngp, &rb, FALSE);

  ResetClip ();
  Update ();
  InsetRect (&rp, 4, 4);
  do_resize_panel (pnl, adp, (Int4)(rp.right - rp.left), (Int4)(rp.bottom - rp.top), rearrange);
  return;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************
***  Draw Scale 
***         draw_scale
***         draw_bars
***
***********************************************************************/
static void draw_scale (EditAlignDataPtr adp, Int4 hoffset, Int4 scalelength, PoinT *ptlh)
{
  Char   str[128];
  Int4   scal;
  Int4   ptx, pty;
  Int4   j;
  Int4   marqueediff;

  if ( !adp->draw_scale ) return;
  SetColor (adp->colorRefs[COLOR_SCALE]);
  ptx = ptlh->x + adp->margin.left + (Int4)(adp->charw * 0.5) - 2;
  pty = ptlh->y + adp->ascent;
  marqueediff = (Int4)( 2.00 / 6.00 * adp->lineheight);
  for ( j = hoffset; j < hoffset + scalelength; j++) 
  {
         if ( adp->colonne[j] > -1) 
         {
                scal = (Int4)(adp->gr.left + 1 + adp->colonne[j]);
                if (scal % 10 == 0) 
                {
                       sprintf (str, "%d", (int)scal);
                       MoveTo ((Int4)(ptx - StringWidth(str) + (adp->charw/2) +1), 
                               (Int4)(pty + marqueediff));
                       PaintString (str);
                }
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int4) j % (Int4) adp->columnpcell == 0) ptx += adp->charw;
         ptx += adp->charw;
  }
  Black();
  return ;
}

static void draw_bars (EditAlignDataPtr adp, Int4 hoffset, Int4 scalelength, PoinT *ptlh)
{
  Int4   scal;
  Int4   ptx;
  Int4   y;
  Int4   j;
  Int4   marqueelong, marqueeshort, marqueediff;

  if ( !adp->draw_bars ) return;
  SetColor (adp->colorRefs[COLOR_SCALE]);
  ptx = ptlh->x + adp->margin.left + (Int4)(adp->charw * 0.5) - 1;
  y = ptlh->y + (Int4) (2.00 / 6.00 * adp->lineheight);
  marqueelong = (Int4) (4.00 / 6.00 * adp->lineheight);
  marqueeshort= (Int4) (2.00 / 6.00 * adp->lineheight);
  marqueediff = (Int4) (2.00 / 6.00 * adp->lineheight);
  for ( j = hoffset; j < hoffset + scalelength; j++) 
  {
         if ( adp->colonne[j] > -1) 
         {
                scal = (Int4)(adp->gr.left + 1 + adp->colonne[j]);
                if (scal % 10 == 0) 
                {
                       MoveTo (ptx , y );
                       LineTo (ptx , (Int4)(y + marqueelong));
                }
                else if (scal % 5 == 0) 
                {
                       MoveTo (ptx , (Int4)(y + marqueediff));
                       LineTo (ptx , (Int4)(y + marqueediff + marqueeshort));
                }
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int4) j % (Int4) adp->columnpcell == 0) ptx += adp->charw;
         ptx += adp->charw;
  }
  Black();
  return ;
}
         
/***********************************************************************
***  draw_id  
************************************************************************/
static void draw_id (EditAlignDataPtr adp, PoinT *pt, Int4 index, CharPtr strid, Uint1 strand, Int4 pos, Uint2 itemtype, Boolean idselected, Boolean is_master, Int4 group)
{
  Char     str[128], str1[128];
  RecT     rct;
  CharPtr  tmp;
  Int4     stringlens;
  Int4     j;
  Int4     total = 0;
  Int4     posspace;

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!*/ idselected = FALSE;

  if (adp->length < 1000) { posspace = 4;
  } else if (adp->length < 10000) { posspace = 5;
  } else if (adp->length < 100000) { posspace = 6;
  } else if (adp->length < 1000000) { posspace = 7;
  } else { posspace = 8;
  }

  if (adp->marginwithindex) total += 3;
  if (adp->marginwithgroup) total += 4;
  if (adp->marginwithpos) total += posspace;
  total += 2;

  str[0] = '\0';
  tmp = str;
  if (adp->marginwithindex && index > 0) {
     sprintf (str1, "%3d ",(int) index);
  } else {
     sprintf (str1, " ");
  }
  tmp = StringMove (tmp, str1);
  tmp[0] = '\0';
  
  adp->marginwithgroup = FALSE;
  if (adp->marginwithgroup) {
     if (group >= 0)
        sprintf (str1, "%4d ",(int) group);
     else sprintf (str1, "    ");
     tmp = StringMove (tmp, str1);
     tmp[0] = '\0';
  }

  if (strid != NULL && strid[0]!='\0') {
     stringlens = (Int4)(StringLen (strid));
/******!!!!!!!!!!!!!!!
     if (stringlens > 4) {
         if (strid[0] == 'l' && strid[1] == 'c' && strid[2] == 'l') {
            strid += 4;
            stringlens -= 4;
         }
     }
!!!!!!!!!!!!!!!!!*********/
     for (j=0; j<stringlens; j++) 
        str1 [j] = strid [j]; 
     str1 [j] = '\0';
     if (stringlens < adp->size_labels) {
        for (j=stringlens; j<adp->size_labels; j++)
              str1 [j] = ' ';
        str1 [j] = '\0';
        stringlens = (Int4)(StringLen (str1));
     }
     stringlens = MIN (adp->marginleft - total, stringlens);
     str1 [stringlens] = '\0';
     tmp = StringMove (tmp, str1);         
  } 
  else if (adp->marginleft > 20) {
     sprintf (str1, "          ");
     tmp = StringMove (tmp, str1);
  }
/*
  symbol[0] = ' ';
  if (strand == Seq_strand_minus) 
     symbol[1] = '<';
  else 
     symbol[1] = '>';
  symbol[2] = '\0';
  tmp = StringMove (tmp, symbol);
*/
  if (adp->marginwithpos) {
     if (posspace <= 4) {
        if (pos > 0)
           sprintf (str1, "%4ld ", (long) pos);
        else sprintf (str1, "     ");
     } else if (posspace == 5) {
        if (pos > 0) 
           sprintf (str1, "%5ld ", (long) pos); 
        else sprintf (str1, "      ");
     } else if (posspace == 6) {
        if (pos > 0)  
           sprintf (str1, "%6ld ", (long) pos);  
        else sprintf (str1, "       "); 
     } else {
        if (pos > 0)  
           sprintf (str1, "%7ld ", (long) pos);
        else sprintf (str1, "        ");
     } 
     tmp = StringMove (tmp, str1);
  }
  *tmp = '\0';

  stringlens = (Int4)(StringLen(str));
  if (stringlens < adp->marginleft -1)
         for (j = stringlens; j < adp->marginleft; j++, tmp++) 
            *tmp = ' '; 
  str [adp->marginleft -1] = ' ';
  str [adp->marginleft] = '\0';
  if ( itemtype != OBJ_BIOSEQ && adp->displaytype )
         SelectColor (162, 163, 82);
  else if (is_master && !adp->all_sequences)
         SetColor (adp->colorRefs[COLOR_ID_MASTER]);
  else 
         SetColor (adp->colorRefs[COLOR_ID]);
  if ( !idselected ) 
  {
Black();
         MoveTo (pt->x, (Int4)(pt->y + adp->ascent));
         PaintString (str);
  }
  else  {
         InvertColors ();
         LoadRect(&rct, pt->x, pt->y, (Int4)(pt->x +stringlens *adp->charw -2), 
                         (Int4)(pt->y + adp->lineheight));
         EraseRect (&rct);
         MoveTo (pt->x, (Int4)(pt->y +adp->ascent));
         PaintString (str);
         InvertColors();
  } 
  Black();
  return;
}

/******************************************************************
***      paint_caret,  draw caret : draws caret
***
*** Now caret is a T between basis of letters
*** Before it was a bar between letters:

  pt1.x = pt.x + ( column + marginleft ) * charw -1;
  pt1.y = pt.y;
  pt2.x = pt1.x;
  pt2.y = pt1.y + lineheight -2;
***
*******************************************************************/
static void paint_caret (PoinT pt, Int4 column, Int4 charw, Int4 lineheight, Int4 marginleft)
{
  PoinT pt1, pt2;
  Int4  lenv = 3;
  Int4  lenh = 2;
  Int4  row;

  Red ();
  WidePen (2);
  row = pt.x + ( column + marginleft ) * charw -1;
  pt1.x = row; 
  pt1.y = pt.y  + lineheight -1;
  pt2.x = row;
  pt2.y = pt1.y  -lenv;
  DrawLine (pt1, pt2);
  pt2.x = row +lenh;
  pt2.y = pt1.y;
  DrawLine (pt1, pt2);
  pt2.x = row -lenh;
  pt2.y = pt1.y;
  DrawLine (pt1, pt2);
  WidePen (1);
  Black ();
} 

static void draw_caret (EditAlignDataPtr adp, SelEdStructPtr sesp, PoinT pt, Int4 line)
{
  SeqLocPtr    slp;
  Int4         column;  
  
  if (!adp->display_panel) 
  {
     if ( is_samess_ses (&(adp->caret), sesp))
     {  
        slp = (SeqLocPtr) adp->caret.region;
        if(SeqPosInLineColumn (SeqLocStart(slp), adp->alignline[line], &column, adp->hoffset, adp) )
        {
           paint_caret (pt, column, adp->charw, adp->lineheight, adp->marginleft);
        }
     }
  }
}

/*********************************************************************
***      PaintSubseq : paint visible sequence  
**********************************************************************/
/*browse a structure to find out the colour number idx (zero-based)
 seq_color is a field in MediaInfo structure; see cn3dmsg.h*/

static ResidueColorCellPtr GiveClrFromIdx(Int4 idx,ValNodePtr seq_color)
{
Int4  i;

	i=0;
	while(seq_color){
		if (i==idx) {
			return((ResidueColorCellPtr)seq_color->data.ptrvalue);	
		}
		seq_color=seq_color->next;
		i++;
	}

   return NULL;
}

static void PaintSubseq (Int4 k, Int4 lg, Int4 from, PoinT *pt, CharPtr vstr, Boolean invert, EditAlignDataPtr adp, Boolean draw_diff, CharPtr *masterstr, Boolean cmplt, ValNodePtr seq_color,SeqIdPtr sip) 
{
  CharPtr      strPtr;
  Char         str[512];
  Uint4        curColor;
  Int4         to = k + lg;
  Int4         caret_color = 0;
  CharPtr      masterstrptr = NULL;
  RecT         rct;
  Uint4        blackColor = GetColorRGB(0,0,0), whiteColor = GetColorRGB(255, 255, 255),
               newColor;

  Boolean      pretty = (Boolean) (seq_color!=NULL);
  ResidueColorCellPtr rgb;
  Int4 from2;	
	
  dashedstring ((CharPtr) str, 512);
  strPtr = str;

  if (pretty) {
  	from2=AlignCoordToSeqCoord (from, sip, 
			(SeqAlignPtr)adp->sap_align->data, adp->sqloc_list, 0);	

	if (from2!=GAP_RESIDUE)
		rgb=GiveClrFromIdx(from2,seq_color);
	else rgb=NULL;

    
	if(rgb != NULL ){
	   curColor = GetColorRGB (rgb->rgb[0], rgb->rgb[1],rgb->rgb[2]);
    }
    else curColor =blackColor;
  }
  else {
     curColor = adp->colorRefs[(Uint1)(vstr[k] - '*')];
  } 
  if (curColor == (Uint4)172)
     curColor = blackColor;

  if (masterstr != NULL) {
     if (*masterstr != NULL) 
        masterstrptr = *masterstr;
  }

  while ( k < to ) 
  {
     if (pretty) {
		from2=AlignCoordToSeqCoord (k+from, sip, 
				(SeqAlignPtr)adp->sap_align->data, adp->sqloc_list, 0);	

		if (from2!=GAP_RESIDUE)
			rgb=GiveClrFromIdx(from2,seq_color);
		else rgb=NULL;
		
    	if(rgb != NULL ){
		   newColor= GetColorRGB (rgb->rgb[0], rgb->rgb[1],rgb->rgb[2]);
    	}
    	else newColor=blackColor ;
		if (vstr[k] == '-') newColor=blackColor ;
     }
     else {
        newColor = adp->colorRefs[(Uint1)(vstr[k] - '*')];
     }
/*   if (newColor == (Uint4)172)
        newColor = blackColor; */ 
                             /* yanli comment it out so that residues */ 
                            /* with magenta color can be seen in salsa */
     if (newColor == whiteColor) newColor = blackColor;
                            /* yanli added this so that residues colored with white */
       /* in Cn3D(with color by Domain) are drawn as black in salsa to be visible */

     if ( adp->colonne [from +k -adp->bufferstart] < 0 ) 
     {
        *strPtr = 0;
        strPtr = str;
        if (cmplt) 
           strPtr = complement_string (strPtr);
        MoveTo ((Int4)(pt->x +adp->margin.left), (Int4)(pt->y + adp->ascent));
	PaintString (strPtr);
        pt->x += caret_color * adp->charw;
        caret_color = 1;
        strPtr = str;
        Black(); 
        if ((from +k -adp->bufferstart) > adp->length) 
           return;
        pt->x += adp->intersalpwidth * adp->charw;
        k += adp->intersalpwidth;
        SetColor (curColor);
     }
     else if (curColor != newColor)
     {
        *strPtr = 0;
        strPtr = str;
        if (cmplt) 
           strPtr = complement_string (strPtr);
        MoveTo ((Int4)(pt->x +adp->margin.left), (Int4)(pt->y +adp->ascent));
        if (invert)
        {
           SetColor (adp->colorRefs[COLOR_SELECT]);
           InvertColors ();
           LoadRect (&rct, (Int4)(pt->x+adp->margin.left), pt->y, (Int4)(pt->x +adp->margin.left+ caret_color*adp->charw), (Int4)(pt->y +adp->lineheight -1));
           EraseRect (&rct);
           InvertColors ();
           if (adp->colorRefs[COLOR_SELECT]==blackColor)
              White();
           else  {
              Black ();
           }
        } 
        else {
           SetColor (curColor);
        }
        PaintString (strPtr);
        pt->x += caret_color *adp->charw;
        caret_color = 1;
        strPtr = str;
        if ( draw_diff && masterstrptr != NULL ) 
        {
           if (vstr[k] !='-') {
              if (ISA_na(adp->mol_type)) {
                 if ((*masterstrptr == vstr[k]) && (vstr[k]<65 || vstr[k]>91))
                    *strPtr = '.';
                 else 
                   *strPtr = vstr[k];
              } 
              else {
                 if (TO_LOWER(*masterstrptr) == TO_LOWER(vstr[k])) {
                    *strPtr = '.';
                 } else {
                    *strPtr = vstr[k];
                 }
              } 
           }
           else 
              *strPtr = vstr[k];
           masterstrptr++;
        } 
        else 
           *strPtr = vstr[k];
        strPtr++;
        if (adp->columnpcell > 0) 
        {
           if ((Int4) (k+1) % (Int4) adp->columnpcell == 0) {
              *strPtr = ' ';
              strPtr++;
              caret_color++;
           }
        }
        curColor = newColor;
        k++;
     } 
     else {
        if ( draw_diff && masterstrptr != NULL ) 
        {
           if (vstr[k] !='-') {
              if (ISA_na(adp->mol_type)) {
                 if ((*masterstrptr == vstr[k]) && (vstr[k]<65 || vstr[k]>91))
                    *strPtr = '.';
                 else 
                    *strPtr = vstr[k];
              } 
              else {
                 if (TO_LOWER(*masterstrptr) == TO_LOWER(vstr[k])) {
                    *strPtr = '.';
                 } else {
                    *strPtr = vstr[k];
                 }
              }
           }
           else 
              *strPtr = vstr[k];
           masterstrptr++;
        } 
        else 
           *strPtr = vstr[k];
        strPtr++;
        caret_color++;
        if (adp->columnpcell > 0) 
        {
           if ((Int4) (k+1) % (Int4) adp->columnpcell == 0) 
           {
              *strPtr = ' ';
              strPtr++;
              caret_color++;
           }
        }
        k++;
     }
  }
  *strPtr = 0;
  strPtr = str;
  if (cmplt) 
     strPtr = complement_string (strPtr);
  MoveTo ((Int4)(pt->x +adp->margin.left), (Int4)(pt->y +adp->ascent));
  if (invert)
  {
     SetColor (adp->colorRefs[COLOR_SELECT]); 
     InvertColors ();
     LoadRect (&rct, (Int4)(pt->x+adp->margin.left), pt->y, (Int4)(pt->x +adp->margin.left+ caret_color*adp->charw), (Int4)(pt->y +adp->lineheight -1));
     EraseRect (&rct);
     InvertColors ();
     if(adp->colorRefs[COLOR_SELECT]==blackColor)
        White();
     else {
        Black ();   
     }
  }
  else {
     SetColor (curColor);
  }
  PaintString (strPtr);
  Black ();
  pt->x += caret_color *adp->charw;
  *masterstr = masterstrptr;
  return;
} 

/*****************************************************************
***  draw_seq
******************************************************************/
static ValNodePtr deleteseqloc (ValNodePtr head,  ValNodePtr delp)
{
  ValNodePtr next = NULL,
             pre = NULL,
             tmp = NULL;
  
  tmp = head; 
  while (tmp!=NULL) 
  {
     next = tmp->next;
     if (tmp == delp) 
     {
        if (pre==NULL) {
           head = next;
        }
        else {
           pre->next = next;
        }
        tmp->next = NULL;
        SeqLocFree ((SeqLocPtr)tmp->data.ptrvalue);
        tmp->data.ptrvalue = NULL;
        ValNodeFree (tmp);
     }
     else 
        pre = tmp;
     tmp = next;
  }
  return head;
}


static ValNodePtr simplifySeqLocList (ValNodePtr valnode)
{
  ValNodePtr   vnpa = NULL,
               vnpb = NULL;
  SeqLocPtr    slpa,
               slpb;
  SeqIntPtr    sint;
  Boolean      check = TRUE,
               loopin = TRUE;

  while (check)
  {
     check = FALSE;
     loopin = TRUE;
     vnpa = valnode;
     while (vnpa != NULL && loopin)
     {
           for (vnpb = vnpa->next; vnpb != NULL; vnpb = vnpb->next)
           {
                 slpa = (SeqLocPtr)vnpa->data.ptrvalue;
                 slpb = (SeqLocPtr)vnpb->data.ptrvalue;
                 if (SeqLocCompare (slpa, slpb) == SLC_A_IN_B) {
                    valnode = deleteseqloc (valnode, vnpa);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_B_IN_A) {
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_EQ_B) {
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_OVERLAP_B) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    if (SeqLocStart(slpa) < SeqLocStart(slpb)) {
                       sint->to = SeqLocStop(slpb);
                    }
                    else {
                       sint = (SeqIntPtr) slpa->data.ptrvalue;
                       sint->from = SeqLocStart(slpb);
                    }
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocStop(slpa) == SeqLocStart(slpb)-1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->to = SeqLocStop(slpb);
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocStart(slpa) == SeqLocStop(slpb)+1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->from = SeqLocStart(slpb);
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
           }
           if (loopin && vnpa != NULL)
              vnpa = vnpa->next;
     }
  }
  return valnode;
}


static Int4 getnextpos (ValNodePtr vnp, Int4 *stop)
{
  ValNodePtr tmp;
  SeqLocPtr  slp;
  SeqIntPtr  sint;
  Int4       valmin = (Int4)(INT4_MAX-2), 
             val = (Int4)-1,
             val2 = (Int4)-1;

  for (tmp=vnp; tmp!=NULL; tmp=tmp->next) 
  {
     slp = (SeqLocPtr) tmp->data.ptrvalue;
     if (SeqLocStart(slp) < valmin) {
        valmin = SeqLocStart(slp);
     }
  }
  if (valmin < (Int4)(INT4_MAX-2)) 
  {
     for (tmp=vnp; tmp!=NULL; tmp=tmp->next) 
     {
        slp = (SeqLocPtr) tmp->data.ptrvalue;
        if (SeqLocStart(slp) == valmin) {
           val = SeqLocStart(slp);
           val2 = SeqLocStop (slp);
           sint = (SeqIntPtr) slp->data.ptrvalue;
           sint->from = INT4_MAX;
           sint->to = INT4_MAX;
           break;
        }
     }
  }
  *stop = val2;
  return val;
}

static void draw_seq (EditAlignDataPtr adp, SelEdStructPtr sesp, Int4 from, Int4 drw_width, PoinT ptlh, CharPtr drw_str, Boolean draw_diff, CharPtr masterstr, Boolean cplmt, Uint2 itemtype, Boolean check_selection)
{
  SelStructPtr  ssptmp;
  SeqIdPtr      sip;
  Int4          start, stop;
  Int4          left, right;
  Int4          seg_lg;
  Int4          k,
                k1;
  Int4          chklocp;
  Int4          iCount = 0;
  SeqLocPtr     slp;
  ValNodePtr    startp = NULL;

SeqIdPtr sesp_sip= NULL;
SeqIntPtr sinp = NULL;

	SeqEditViewProcsPtr svpp;
	ValNodePtr          vnp_seqinfo;
	ValNodePtr          vnp_color = NULL;
	SeqIdPtr            vnpcn3d_sip = NULL;
	MediaInfoPtr        MIPtr;


  adp->start_select = -1;
  start = stop = -1;
  ssptmp = NULL;
  if (check_selection) 
  {
    ssptmp = ObjMgrGetSelected (); 
    for (; ssptmp != NULL; ssptmp = ssptmp->next) 
    {
     if ( checkssp_for_editor (ssptmp) && is_samess_ses (ssptmp, sesp) ) 
     {
         sip = SeqLocId ((SeqLocPtr) ssptmp->region);
         start = SeqLocStart ((SeqLocPtr) ssptmp->region);
         adp->start_select = start;
         chklocp =chkloc(sip, start, adp->sqloc_list, &start);
         start = SeqCoordToAlignCoord (start, sip, (SeqAlignPtr) adp->sap_align->data, 0, chklocp);
         stop  = SeqLocStop  ((SeqLocPtr) ssptmp->region);
         chklocp =chkloc(sip, stop, adp->sqloc_list, &stop);
         stop = SeqCoordToAlignCoord (stop, sip, (SeqAlignPtr) adp->sap_align->data, 0, chklocp);
/******
         if (start<=stop && (start < from + drw_width && (stop + 1) > from)) 
******/
         if ((start < from + drw_width && (stop + 1) > from)) 
         {
            if(start<=stop)
               slp = SeqLocIntNew (start, stop, Seq_strand_plus, sip); 
            else
               slp = SeqLocIntNew (stop, start, Seq_strand_plus, sip); 
            ValNodeAddPointer (&startp, 0, (Pointer)slp);
            iCount ++;
         }
     }
    }
  }


  svpp = (SeqEditViewProcsPtr) GetAppProperty ("SeqEditDisplayForm");

  if(sesp->regiontype == 1){
     slp = (SeqLocPtr)sesp->region;
     sinp = (SeqIntPtr)slp->data.ptrvalue;
     sesp_sip = sinp->id;
  }               
  if (svpp) {
     if (svpp->Cn3D_On && sesp_sip) {
	vnp_seqinfo = svpp->seqinfo;

	while(vnp_seqinfo){
		MIPtr=(MediaInfoPtr)vnp_seqinfo->data.ptrvalue;
/*      if (MIPtr->entityID==sesp->entityID && MIPtr->itemID==sesp->itemID){
			vnp_color=MIPtr->seq_color;
                        vnpcn3d_sip = MIPtr->sip;
  	        break;  
		}         */    /* yanli comment it out, instead doing the following */

        if(SeqIdForSameBioseq(sesp_sip, MIPtr->sip)){
           vnp_color=MIPtr->seq_color;
           vnpcn3d_sip = MIPtr->sip;
                     /* do not break so as to get the latest one */
        }    
		vnp_seqinfo=vnp_seqinfo->next;
	}
     }
  }
  if (iCount > 0) {
     k = 0;
     k1 = from;
     startp = simplifySeqLocList (startp);
     start = getnextpos (startp, &stop);
     while ( k < drw_width && start >= 0) {
        if ( k + from < start ) 
        {
           left = (Int4) k1;
           right= start;
           seg_lg = MIN (right - left, drw_width);
           PaintSubseq (k, seg_lg, from, &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, vnp_color,vnpcn3d_sip);
           k += seg_lg;
        } 
        else if ( k +from >= start && k +from < stop + 1) 
        {
           left = (Int4) MAX (k1, start);
           right= (Int4) MIN (from + drw_width, stop + 1);
           seg_lg = MIN (right - left, drw_width);
           PaintSubseq (k, seg_lg, from,  &ptlh, drw_str, TRUE, adp, draw_diff, &masterstr, cplmt, vnp_color,vnpcn3d_sip);
           k += seg_lg;
           k1 = stop + 1;
           start = getnextpos (startp, &stop);
        }         
     }
     if ( k < drw_width) 
     {
        left = (Int4) k1;
        right= (Int4) from + drw_width;
        seg_lg = MIN (right - left, drw_width);
        PaintSubseq (k, seg_lg, from,  &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, vnp_color,vnpcn3d_sip);
     }
  }
  else {
     PaintSubseq (0, drw_width, from, &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, vnp_color,vnpcn3d_sip);
  }
  if (startp != NULL)
     ValNodeFree (startp);
  return;
}

/*****************************************************************
***  draw_line
******************************************************************/
static void draw_line (EditAlignDataPtr adp, Int4 start, Int4 stop, Uint1 strand, PoinT *pt, Int4 from, Int4 drw_width, Int4 line, Int4 alignline, Uint4 color, Uint2 wideline, Boolean selected, Boolean partialstart, Boolean partialstop, BoolPtr gapline)
{
  PoinT         pt1, pt2,
                pttmp;
  PoinT         oldpt1, oldpt2;
  FloatLo       hgt = (FloatLo)(2.5 / 4.0);
  RecT          leftrct, rightrct;
  Int4          left, right;
  Int4          yline;  
  Int4          column;  
  Int4          above, below;
  Int4          j, k;
  Int4          col1;
  Uint1         startin, stopin;
  BoolPtr       gap;
  
  left = (Int4) MAX ((Int4) from, (Int4) start);
  SeqPosInLineColumn (left, alignline, &column, adp->hoffset, adp);
  yline = pt->y + (Int4) (hgt * adp->ascent);
  pt1.x = pt->x + (column + adp->marginleft) * adp->charw;
  pt1.y = yline;

  if (adp->columnpcell > 0) {
     col1 = column - (column-1)/adp->columnpcell;
  }
  right = (Int4) MIN ((Int4)from+drw_width-1, (Int4)stop);
  SeqPosInLineColumn (right, alignline, &column, adp->hoffset, adp);
  pt2.x = pt->x + (column + adp->marginleft + 1) * adp->charw;
  pt2.y = yline;

  startin = stopin = 0;
  if ( start >= from ) {
     if (strand == Seq_strand_minus)  
        startin = 1;
     else
        startin = 2;
  } 
  if (stop < from+drw_width && stop != start) {
     if (strand == Seq_strand_minus) 
        stopin = 1;
     else 
        stopin = 2;
  }
  oldpt1.x = pt1.x;
  oldpt1.y = pt1.y;
  oldpt2.x = pt2.x;
  oldpt2.y = pt2.y;
  if (start<stop) {
   if (startin)
     pt1.x += 5;
  
   SetColor (color);
   if (gapline == NULL) {
     if (stopin)
        pt2.x -= 5;
     WidePen (wideline);
     DrawLine (pt1, pt2);
   } 
   else {
     pttmp.y = yline;
     j=k=left;
     for (j=0, gap=gapline; j<col1; j++) 
        gap++;
     while (j<right) {
        if (*gap) {
           while (*gap && k<right) {
              gap++;
              k++;
           }
           WidePen (1);
           Dashed ();
        }
        else {
           while (!(*gap) && k<right) {
              gap++;
              k++;
           }
           WidePen (wideline);
        }
        SeqPosInLineColumn (k, alignline, &column, adp->hoffset, adp);
        if (k==right && stopin == 0)  {
           column++;
        }
        pttmp.x = pt->x + (column + adp->marginleft) * adp->charw;
        DrawLine (pt1, pttmp);
        Solid ();
        pt1.x = pttmp.x;
        j = k;
     }
   }
  }
  if (selected) {
#ifdef WIN_MAC
         above = 2;
         below = 4;
#endif
#ifdef WIN_MSWIN
         above = 3;
         below = 3;
#endif
#ifdef WIN_MOTIF
         above = 3;
         below = 3;
#endif
         WidePen (1);
         Black ();
         pt1.y = pt2.y = yline + below;
         DrawLine (pt1, pt2);
         pt1.y = pt2.y = yline - above;
         DrawLine (pt1, pt2);
         White ();
         pt1.y = pt2.y = yline + below - 1;
         DrawLine (pt1, pt2);
         pt1.y = pt2.y = yline - above + 1;
         DrawLine (pt1, pt2);
  }
  WidePen (1);
  SetColor (color);
  pt1.x = oldpt1.x;
  pt1.y = oldpt1.y;
  pt2.x = oldpt2.x;
  pt2.y = oldpt2.y;
  if (startin == 1) {
#ifdef WIN_MAC
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int4)(pt1.y +4) );
     OffsetRect (&leftrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4) );
#endif
#ifdef WIN_MOTIF
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2) (pt1.y +4) );
#endif
     if (partialstart)
        CopyBits (&leftrct, leftTriSym);
     else 
        CopyBits (&leftrct, leftTriFillSym);
  } 
  else if (startin == 2) {
#ifdef WIN_MAC
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
     OffsetRect (&leftrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
#endif
     if (partialstart)
        CopyBits (&leftrct, rectSym);
     else
        CopyBits (&leftrct, rectFillSym);
  }
  if (stopin == 1) {
#ifdef WIN_MAC
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
     OffsetRect (&rightrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
     if (partialstop)
        CopyBits (&leftrct, rectSym);
     else 
        CopyBits (&rightrct, rectFillSym);
  }
  else if (stopin == 2) {
#ifdef WIN_MAC
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
     OffsetRect (&rightrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
     if (partialstop)
        CopyBits (&rightrct,rightTriSym);
     else 
        CopyBits (&rightrct,rightTriFillSym);
  }
  Black ();
  WidePen (1);
}

/*********************************************************************
***  draw_cds
*********************************************************************/
static void draw_trans (EditAlignDataPtr adp, SelEdStructPtr cds, CharPtr trans, PoinT *pt, Int4 from, Int4 drw_width, Uint2 offset)
{
  SeqLocPtr     slp;
  CharPtr       transtr;
  CharPtr       strPtr, transPtr;
  Int4          left, right;
  Int2          k;  

  if ( trans == NULL || cds == NULL ) return;
  slp = (SeqLocPtr) cds->region;
  if (slp == NULL) return;
  left  = (Int4) MAX ((Int4) from, (Int4) SeqLocStart(slp));
  right = (Int4) MIN ((Int4) from+drw_width, (Int4) SeqLocStop(slp)+1);
  if (right < left) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in draw_trans: left %ld right %ld", (long) left, (long) right);
  }
  transtr = (CharPtr)MemNew ((size_t)(512 * (sizeof (Char))));
  emptystring (transtr, 512);
  transtr[0] = ' ';
  transtr[511] = '\0';
  strPtr = transtr;
  strPtr += (Int2) (left - from);
  transPtr = trans;
  transPtr += cds->offset;   /*!!!!!!!!!!!!! += cds->codonstart ?? */
  if (SeqLocStart(slp) < from) 
     transPtr += (Int2)(from - SeqLocStart(slp));  
  if (adp->prot_mode == MPROTAA) {
     for (k =0; k < (right - left); k++, strPtr++, transPtr++) {
        if ( *transPtr=='M' || *transPtr=='*') *strPtr = *transPtr;
     }
  } 
  else {
     for (k =0; k < (right - left); k++, strPtr++, transPtr++) {
        *strPtr = *transPtr;
     }
  }
  strPtr++;
  *strPtr= '\0';
  draw_seq (adp, cds, from, drw_width, *pt, transtr, FALSE, NULL, TRUE, OBJ_SEQFEAT, FALSE);
  transtr = (CharPtr)MemFree (transtr);
  return;
}

/*********************************************************************
***  draw_pept
*********************************************************************/
static void draw_pept (EditAlignDataPtr adp, SelEdStructPtr cds, PoinT *pt, Int4 from, Int4 drw_width, Int2 line, Uint2 offset)
{
  ValNodePtr  pept;
  
  if (cds!=NULL) {
     pept = cds->data;
     if (pept != NULL)  { 
        if (pept->data.ptrvalue != NULL)  { 
           draw_trans (adp, cds, (CharPtr) pept->data.ptrvalue, pt, from, drw_width, offset);
        }
     }
  }
}

/*********************************************************************
***  draw_feat
*********************************************************************/
static void draw_feat (EditAlignDataPtr adp, Uint2 eid, Uint2 iid, Uint2 it, Uint2 choice, Int4 start, Int4 stop, Uint1 strand, CharPtr label, PoinT *pt, Int4 from, Int4 drw_width, Int2 line, Int4 alignline, Uint4 color, Boolean partialstart, Boolean partialstop, BoolPtr gapline)
{
  Char      str[128];
  Uint2     wideline;
  Boolean   featselect;

  if ( start > stop || start >= from + drw_width || stop < from ) 
     return;
  if (it == OBJ_VIRT) 
  {
     wideline = 5;
     if (choice == SEQFEAT_CDREGION ) 
         sprintf (str, "CDS %d", (int) iid);
     else if (choice == SEQFEAT_GENE) 
         sprintf (str, "gene %d", (int) iid);
     else if (choice == SEQFEAT_RNA)
         sprintf (str, "mRNA %d", (int) iid);
     str[15] = '\0';
     draw_id (adp, pt, -1, str, strand, 0, 3, FALSE, FALSE, -1);
  }
  else {
     wideline = 3;
     draw_id (adp, pt, -1, label, strand, 0, 3, FALSE, FALSE, -1);
  }
  featselect = (Boolean) (is_selectedbyID (eid, iid, it) != NULL);
  draw_line (adp, start, stop, strand, pt, from, drw_width, line, alignline, color, wideline, featselect, partialstart, partialstop, gapline);
}

static void what_inline (EditAlignDataPtr adp, Int4 line, Uint2 ei, Uint2 ii, Uint2 it, Uint2 ist, Uint2 al)
{ 
  adp->seqEntity_id [line]  = ei;
  adp->item_id [line]  = ii;
  adp->itemtype[line] = it;
  adp->itemsubtype[line] = ist;
  adp->alignline [line]= al;
}

static Boolean id_is_invalid (EditAlignDataPtr adp, PoinT ptlh)
{
  RecT  rect; 
  LoadRect(&rect, ptlh.x, ptlh.y, 
                  (Int2)(ptlh.x + adp->marginleft * adp->charw), 
                  (Int2)(ptlh.y + adp->lineheight));
  return RectInRgn (&rect, updateRgn);
}

static Boolean row_is_invalid (EditAlignDataPtr adp, PoinT ptlh)
{
  RecT  rect; 
  LoadRect (&rect, (Int2)(ptlh.x + adp->margin.left - adp->charw), ptlh.y,  
                   (Int2)(ptlh.x + adp->pnlWidth * adp->charw + 1), 
                   (Int2)(ptlh.y + adp->lineheight));
  return RectInRgn (&rect, updateRgn);
}

static Int2 rang (Uint2 format, SelEdStructPtr ssp)
{
  if (format == OBJ_BIOSEQ)
     return ssp->entityID;
  if (format == OBJ_SEQALIGN)
     return ssp->itemID;
  return 0;
}

static Int2 getparam (ValNodePtr vnprm, Uint2 eID, Uint2 iID, Uint1 choice)
{
  ValNodePtr   vnp;
  SeqParamPtr  prm;

  for (vnp = vnprm; vnp != NULL; vnp = vnp->next) {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == eID && prm->itemID == iID)
     {
        if (choice == 1)
           return prm->group;
     }
  }
  return -1;
}

static CharPtr seqid_tolabel (SeqIdPtr sip, Uint2 choice)
{
  BioseqPtr bsp;
  SeqIdPtr  tmp = NULL;
  Char      str[120];
  CharPtr   strp;
  Uint4     val;

  str[0] = '\0';
  if (choice && choice <= PRINTID_GIcc)
  {
     bsp = BioseqLockById(sip);
     if (bsp!=NULL) 
     {
        if (choice==PRINTID_GIcc)
           tmp = SeqIdFindBest (bsp->id, 0);
        else 
           tmp = bsp->id;
        if (tmp)
        {
           SeqIdWrite (tmp, str, PRINTID_TEXTID_ACCESSION, 120);
           val=WHICH_db_accession(str);
           if (val==0)
              choice = PRINTID_FASTA_SHORT;

           SeqIdWrite (tmp, str, choice, 120);
           BioseqUnlock(bsp);
           strp = StringSave (str);
           return strp;
        }
        BioseqUnlock(bsp);
     }
  }
  if (choice!=0 && sip) {
     SeqIdWrite (sip, str, choice, 120);
     strp = StringSave (str);
     return strp;
  } 
  strp = StringSave (str);
  return strp;
}

/*****************************************************************
***
***     Draw Alignment in Panel p 
***
******************************************************************/
static Char strcmpl[] = "complement";

extern void on_draw (PaneL p)
{
  PoinT            ptlh;        /* PoinT at the top of lineheight */
  RecT             rp;          /* Panel Rect */
  EditAlignDataPtr adp;
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  SeqLocPtr        curslp;
  SeqLocPtr        sip;
  SeqAlignPtr      salp;

  Int4          from_inbuf;  /* alignment coordinates in buffer */
  Int4          from_inseq;  /* alignment coordinates in buffer */
  Int4          stringlens;
  Int4          start,  stop;
  Int4          start2, stop2;

  Int4          line;        /* current line */
  Int4          curalgline;
  Int4          drw_width;   /* length of drw_str */
  Int4          group;
  Int4          chklocp,
                chklocp2;
  Uint2         itemsubtype;
  Uint4         offset;
 
  CharPtr       stridp;
  CharPtr       curstr = NULL;
  CharPtr       masterbuf =NULL, 
                masterstr =NULL;
  Int2          oldstyle;

  Boolean       invalid_id;
  Boolean       invalid_row;
  Boolean       is_master, draw_diff;
  Boolean       first = TRUE;
  Boolean       line_drawn;
  Boolean       empty_line;

  Int1          gapinline;
  BoolPtr       gapline,
                gaplinep;
  
  if ( (adp = GetAlignDataPanel (p)) == NULL ) 
     return;
  if ( adp->seqnumber == 0 ) 
     return;

  salp = (SeqAlignPtr) adp->sap_align->data;

  if ( adp->firstssp == NULL ) {
     adp->firstssp = get_firstline (NULL, adp->buffer);
  }
  if ( adp->firstssp == NULL ) 
     return;
  if ( adp->firstssp->region == NULL ) 
     return;
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int4) adp->visibleWidth, 
                   (Int4) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 || adp->bufferstart < 0 ) {
         return;
  }

  oldstyle = GetMuskCurrentSt ();
  SetMuskCurrentSt (GetMuskStyleName (adp->styleNum));

  line = (drw_width +2) * sizeof(Boolean);
  gapline = (BoolPtr) MemNew ((size_t)(line));
  MemSet ((Pointer)gapline, (int)(FALSE), (size_t)(line)); 
  gaplinep = NULL;
  gapinline = LINE_NOGAP;
  empty_line = FALSE;
 
  curssp = (SelEdStructPtr) adp->firstssp->region;
  if (curssp->itemtype ==OBJ_SEQFEAT) {
     for (curvnp=adp->firstssp; curvnp!=NULL; curvnp=curvnp->prev) {
        curssp = (SelEdStructPtr) curvnp->region;
        if (curssp->itemtype==OBJ_BIOSEQ)   /*&&curvnp->itemtype==FEATDEF_BAD)*/
        {
           if (curssp->data !=NULL) {
              vnp = (ValNodePtr) curssp->data;
              curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
              curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
              if ( curstr != NULL )
              {
                 gapinline=getgapsfromstring(curstr, 0, drw_width, &gapline);
                 empty_line =(Boolean)(gapinline == LINE_ONLYGAP);
                 if (gapinline == LINE_NOGAP)
                    gaplinep = NULL;
                 else
                    gaplinep = gapline;
              }
           }
           break; 
        }
     }
  }
  get_client_rect (p, &rp);
  SelectFont ((FonT)(adp->font));
  ptlh.x = rp.left;
  ptlh.y = rp.top;
  line = 0;
  curalgline = 0;
  adp->visibleLength = drw_width;
  masterbuf = get_master (adp->linebuff, adp->master.entityID, 
                          adp->master.itemID, adp->master.itemtype);

  curvnp = adp->firstssp;
  while ( line < adp->pnlLine )
  {

     curssp = (SelEdStructPtr) curvnp->region;
     if (curssp !=NULL) {
         adp->lastses = curssp;
         itemsubtype = curvnp->itemtype;
         invalid_row = row_is_invalid(adp, ptlh);
         invalid_id  = id_is_invalid (adp, ptlh);
         if (itemsubtype == EDITDEF_SCA) 
         {
                what_inline(adp, line, LINE0, LINE0, LINE0, LINE0, curalgline);
                if ( invalid_row ) {
                   draw_scale (adp, from_inbuf, adp->visibleWidth, &ptlh);
                }
                line++;
                ptlh.y += adp->lineheight;
         }
         else if (itemsubtype == EDITDEF_SCB) 
         {
                what_inline(adp, line, LINE0, LINE0, LINE0, LINE0, curalgline);
                if ( invalid_row ) {
                   draw_bars (adp, from_inbuf, adp->visibleWidth, &ptlh);
                }
                line++;
                ptlh.y += adp->lineheight;
         }
         else if (itemsubtype>=EDITDEF_RF1 && itemsubtype<=EDITDEF_RF6) 
         {
            if (curssp->data !=NULL  && !empty_line) {
               curslp = (SeqLocPtr) curssp->region;
               vnp = (ValNodePtr) curssp->data;
               if ( SeqLocStop (curslp) > from_inseq && vnp != NULL) 
               {
                  what_inline (adp, line, curssp->entityID, curssp->itemID,
                                curssp->itemtype, itemsubtype, curalgline);
                  curstr = (CharPtr) vnp->data.ptrvalue;
                  offset = (Uint2)(itemsubtype - EDITDEF_RF1) % (Uint2)3;
                  if (itemsubtype>=EDITDEF_RF1 && itemsubtype<EDITDEF_RF4)
                     draw_id (adp, &ptlh, rang(adp->input_format, curssp), NULL,  Seq_strand_plus, 0, curssp->itemtype, FALSE, FALSE, -1);
                  else
                     draw_id (adp, &ptlh, rang(adp->input_format, curssp), NULL, Seq_strand_minus, 0, curssp->itemtype, FALSE, FALSE, -1);

                  if ( invalid_row ) 
                  {
                    draw_trans (adp, curssp, curstr, &ptlh, from_inseq, drw_width, offset);
                  }
               }
               line++;
               ptlh.y += adp->lineheight;
            }
         }
         else if ( itemsubtype == FEATDEF_TRSL ) 
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while ( curssp != NULL && curssp->data != NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp);
                  if (start < from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) {
                        draw_pept (adp, curssp, &ptlh, from_inseq, drw_width, line, 0);
                     }
                     if (first) {
                        what_inline (adp, line, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if ( itemsubtype == EDITDEF_CPL ) 
         {
            if (curssp->data !=NULL && !empty_line) 
            {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring (curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                         /* !!!!!!! */
                {
                   what_inline (adp, line, curssp->entityID, curssp->itemID,  
                            curssp->itemtype, itemsubtype, curalgline);
                   if ( invalid_id && (( !adp->displaytype 
                   && curssp->itemtype != OBJ_BIOSEQ) || adp->displaytype)) 
                      draw_id (adp, &ptlh, -1, strcmpl, 255, 0, curssp->itemtype, FALSE, FALSE, -1);
                   if ( invalid_row && (( !adp->displaytype 
                   && curssp->itemtype != OBJ_BIOSEQ) || adp->displaytype) ) 
                   {
                      draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, TRUE, itemsubtype, FALSE);
                   }
                   line++;
                   ptlh.y += adp->lineheight;
                }
            }
         }
         else if (curssp->itemtype == OBJ_BIOSEQ)
         {
            if (curssp->data !=NULL) {
               vnp = (ValNodePtr) curssp->data;
               curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
               curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
               if ( curstr != NULL )         
               {
                 gapinline=getgapsfromstring(curstr, 0, drw_width, &gapline);
                 empty_line =(Boolean)(gapinline == LINE_ONLYGAP);
                 if (gapinline == LINE_NOGAP) 
                    gaplinep = NULL;
                 else
                    gaplinep = gapline;
                 if (adp->draw_emptyline || (!adp->draw_emptyline && !empty_line)) {
                   what_inline (adp, line, curssp->entityID, curssp->itemID,  
                                curssp->itemtype, itemsubtype, curalgline);
                   is_master = is_samess_ses (&(adp->master), curssp);
                   draw_diff=!(!adp->charmode || (adp->charmode && is_master));
                   masterstr = NULL;
                   if ( draw_diff ) {
                      stringlens = StringLen (masterbuf);
                      if ( from_inbuf < stringlens )
                            masterstr = masterbuf + from_inbuf;
                   }
                   if ( invalid_id ) {
                     if (adp->displaytype) 
                     {
                        if (adp->marginwithgroup) 
                           group=getparam(adp->params,curssp->entityID, curssp->itemID, 1);
                        else group = -1;
                        sip = SeqLocId((SeqLocPtr) curssp->region);
                        start=AlignCoordToSeqCoord(from_inseq, sip,salp, adp->sqloc_list, 0) +1;
                        stridp = seqid_tolabel (sip, adp->printid);
                        /* seqid_to_label calls BioseqLockById, which will
                         * refresh the Desktop if it is open and will change
                         * the font setting - need to put the font back.
                         */
                        SelectFont ((FonT)(adp->font));
                        draw_id (adp, &ptlh, rang(adp->input_format, curssp), stridp, curtdp->strand, start, curssp->itemtype, FALSE, is_master, group);
                     }
                   }
                   if ( invalid_row && adp->displaytype )
                   {
                      draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, FALSE, curssp->itemtype, TRUE);
                      draw_caret (adp, curssp, ptlh, line);
                   }
                   line++;
                   ptlh.y += adp->lineheight;
                 }
               }
            }
         }
         else if (curssp->itemtype==OBJ_SEQFEAT && itemsubtype==FEATDEF_PROT)
         {
            if (curssp->data !=NULL && !empty_line) 
            {
               vnp = (ValNodePtr) curssp->data;
               curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
               curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
               if ( curstr != NULL )                   
               {
                  what_inline (adp, line, curssp->entityID, curssp->itemID,  curssp->itemtype, itemsubtype, curalgline);
                  draw_id (adp, &ptlh, -1, curtdp->label, curtdp->strand, 0, curssp->itemtype, FALSE, FALSE, -1);
                  if ( invalid_row ) {
                     draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, FALSE, curssp->itemtype, FALSE);
                  }
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if ( curssp->itemtype ==OBJ_SEQFEAT )
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while (curssp!=NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp = chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp2 = chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp2);
                  if (start <= stop && start < from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) 
                     {
                        draw_feat (adp, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, start, stop, SeqLocStrand(curslp), curssp->label, &ptlh, from_inseq, drw_width, line, curalgline, (Uint4) getcolor_fromstyles(itemsubtype), (Boolean)(chklocp==-1), (Boolean)(chklocp2==APPEND_RESIDUE), gaplinep);
                     }
                     if (first) {
                        what_inline(adp, line, curssp->entityID, curssp->itemID,
                                     curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if( itemsubtype==SEQFEAT_GENE || itemsubtype==SEQFEAT_RNA
              ||  itemsubtype==SEQFEAT_CDREGION)
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while (curssp!=NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp2 =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp);
                  if (start <= from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) 
                     {
                        draw_feat (adp, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, start, stop, SeqLocStrand(curslp), curssp->label, &ptlh, from_inseq, drw_width, line, curalgline, (Uint4) getcolorforvirtfeat(curssp->itemID), (Boolean)(chklocp==-1), (Boolean)(chklocp2==APPEND_RESIDUE), gaplinep);
                     }
                     if (first) {
                        what_inline (adp, line, curssp->entityID, curssp->itemID,
                                     curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else {
/**
            ErrPostEx (SEV_ERROR, 0, 0, "fail in on_draw [46] subtype %ld type %ld", (long) itemsubtype, (long) curssp->itemtype); 
**/
            return;
         }
         if (line == adp->pnlLine) 
            break;
         curvnp = curvnp->next;
         if (curvnp == NULL) 
         {
            if (adp->rowpcell > 0 
            && ((Int4)(curalgline+1) % (Int4)adp->rowpcell) == 0) 
            {
                what_inline(adp,line,LINE0,LINE0, LINE0, LINE0, curalgline);
                ptlh.y += adp->lineheight;
                line++;
            }
            curvnp = adp->buffer;
            curalgline++;
/*****************
*******************/
if (curalgline>500)
break;
/*****************
*******************/

            from_inbuf += drw_width;
            from_inseq += drw_width;
            if (from_inseq >= adp->length) 
                break;
            drw_width = MIN ((Int4) adp->visibleWidth, (Int4)(adp->bufferlength 
                      - (adp->hoffset -adp->bufferstart) - adp->visibleLength));
            adp->visibleLength += drw_width;
         }
     }
  }
  MemFree (gapline);
  SetMuskCurrentSt (GetMuskStyleName (oldstyle));
  return;
}

/* This section of code is used for feature propagation */
typedef struct fprdata {
  FORM_MESSAGE_BLOCK
  BioseqPtr           bsp;
  SeqAlignPtr         salp;
  Uint2               selFeatItemID;
  Int4                aln_length;
  Int4                log10_aln_length;
  VieweR              details;
  SegmenT             dtpict;
  Int4                scaleX;
  GrouP               allOrSel;
  GrouP               gapSplit;
  DialoG              sequence_list_dlg;
  ButtoN              stopCDS;
  ButtoN              transPast;
  ButtoN              fixCDS;
  ButtoN              fuseJoints;
  ButtoN              accept;
} FprData, PNTR FprDataPtr;

typedef struct ivalinfo {
  Int4  start1;
  Int4  stop1;
  Int4  start2;
  Int4  stop2;
  struct ivalinfo PNTR next;
} IvalInfo, PNTR IvalInfoPtr;

static IvalInfoPtr IvalInfoFree (
  IvalInfoPtr ival
)

{
  IvalInfoPtr  next;

  while (ival != NULL) {
    next = ival->next;
    MemFree (ival);
    ival = next;
  }
  return NULL;
}

static IvalInfoPtr GetAlignmentIntervals (SeqAlignPtr sap, Int4 row1, Int4 row2, Int4 from, Int4 to)
{
   AlnMsg2Ptr    amp1;
   AlnMsg2Ptr    amp2;
   IvalInfoPtr  ival;
   IvalInfoPtr  ival_head = NULL;
   IvalInfoPtr  ival_prev = NULL;
   Int4         from_aln;
   Int4         start;
   Int4         stop;
   Int4         tmp;
   Int4         to_aln;
   Int4         seg_i, seg_n, seg_start, seg_stop;

   if (sap == NULL || sap->saip == NULL)
      return NULL;
   AlnMgr2GetNthSeqRangeInSA(sap, row1, &start, &stop);
   if (from < start)
      from = start;
   if (to > stop)
      to = stop;
   from_aln = AlnMgr2MapBioseqToSeqAlign(sap, from, row1);
   to_aln = AlnMgr2MapBioseqToSeqAlign(sap, to, row1);
   if (from_aln > to_aln)
   {
      tmp = from_aln;
      from_aln = to_aln;
      to_aln = tmp;
   }
   seg_n = AlnMgr2GetNumSegs(sap);
   for (seg_i=1; seg_i<=seg_n; seg_i++) {
     AlnMgr2GetNthSegmentRange(sap, seg_i, &seg_start, &seg_stop);
     if (seg_start > to_aln) continue;
     if (seg_stop < from_aln) continue;
     if (seg_start < from_aln) seg_start = from_aln;
     if (seg_stop > to_aln) seg_stop = to_aln;
     
     amp1 = AlnMsgNew2();
     amp1->from_aln = seg_start;
     amp1->to_aln = seg_stop;
     amp1->row_num = row1;
     amp2 = AlnMsgNew2();
     amp2->from_aln = seg_start;
     amp2->to_aln = seg_stop;
     amp2->row_num = row2;
     AlnMgr2GetNextAlnBit(sap, amp1);
     AlnMgr2GetNextAlnBit(sap, amp2);
     if (amp1->type == AM_SEQ && amp2->type == AM_SEQ) {
       ival = (IvalInfoPtr)MemNew(sizeof(IvalInfo));
       ival->start1 = amp1->from_row;
       ival->stop1 = amp1->to_row;
       ival->start2 = amp2->from_row;
       ival->stop2 = amp2->to_row;
       if (ival_head == NULL)
         ival_head = ival_prev = ival;
       else {
         ival_prev->next = ival;
         ival_prev = ival;
       }
     }
     AlnMsgFree2(amp1);
     AlnMsgFree2(amp2);
   }
   return ival_head;
}

static IvalInfoPtr MergeAdjacentIntervals (
  IvalInfoPtr list
)

{
  IvalInfoPtr  curr, last, next;

  if (list != NULL) {
    curr = list->next;
    last = list;
    while (curr != NULL) {
      next = curr->next;
      if (curr->start2 == last->stop2 + 1) {
        last->stop2 = MAX (curr->stop2, last->stop2);
        MemFree (curr);
        last->next = next;
      } else {
        last = curr;
      }
      curr = next;
    }
  }
  return list;
}

static SeqLocPtr MergeAdjacentSeqLocIntervals (SeqLocPtr location)

{
  SeqLocPtr    head = NULL;
  SeqIntPtr    sinp, last_sinp;
  SeqLocPtr    slp;
  SeqPntPtr    spp, last_spp;
  SeqLocPtr    last_slp = NULL;
  Int4         last_from = -1, last_to = -1;
  Uint1        last_strand = 0;
  SeqIdPtr     last_id = NULL;

  head = NULL;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    switch (slp->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
	    if (last_slp != NULL && last_strand == sinp->strand && last_id != NULL
			&& SeqIdComp (sinp->id, last_id) == SIC_YES
			&& ((last_from <= sinp->from + 1 && last_to >= sinp->from - 1)
			  || (last_from <= sinp->to + 1 && last_to >= sinp->to - 1)
			  || (last_from >= sinp->from - 1 && last_to <= sinp->to + 1)
			  || (last_from <= sinp->from + 1 && last_to >= sinp->to - 1)))
		{
		  /* intervals are adjacent, so expand previous interval */
  	      if (last_slp->choice == SEQLOC_INT) {
		    last_sinp = last_slp->data.ptrvalue;
		    if (last_sinp != NULL) {
		      last_from = MIN (sinp->from, last_sinp->from);
			  last_sinp->from = last_from;
			  last_to = MAX (sinp->to, last_sinp->to);
			  last_sinp->to = last_to;
			}
		  } else if (last_slp->choice == SEQLOC_PNT) {
		    /* change previous entry from point to interval that includes point */
		    spp = (SeqPntPtr)last_slp->data.ptrvalue;
		    if (spp != NULL) {
              last_sinp = SeqIntNew ();
			  last_from = MIN (spp->point, sinp->from);
			  last_sinp->from = last_from;
			  last_to = MAX (sinp->to, spp->point);
			  last_sinp->to = last_to;
			  last_sinp->strand = sinp->strand;
			  last_sinp->id = SeqIdDup (sinp->id);
			  last_slp->choice = SEQLOC_INT;
			  last_slp->data.ptrvalue = last_sinp;
			  SeqPntFree (spp);
			}
		  }		  	
		} else {
		  /* add new interval */
		  last_sinp = SeqIntNew ();
          last_sinp->from = sinp->from;
		  last_from = sinp->from;
          last_sinp->to = sinp->to;
		  last_to = sinp->to;
          last_sinp->strand = sinp->strand;
		  last_strand = sinp->strand;
		  last_sinp->id = SeqIdDup (sinp->id);
          last_id = last_sinp->id;
          last_slp = ValNodeAddPointer (&head, SEQLOC_INT, (Pointer) last_sinp);
        }
	  }
	  break;
	case SEQLOC_PNT:
	  spp = (SeqPntPtr)slp->data.ptrvalue;
	  if (spp != NULL) {
	    if (last_slp != NULL && last_strand == spp->strand && last_id != NULL
			&& SeqIdComp (spp->id, last_id) == SIC_YES
			&& last_to >= spp->point - 1)
		{
		  /* intervals are adjacent, so expand previous interval */
		  if (last_slp->choice == SEQLOC_INT) {
		    last_sinp = last_slp->data.ptrvalue;
			if (last_sinp != NULL) {
			  last_to = MAX (spp->point, last_sinp->to);
			  last_sinp->to = last_to;
			  last_from = MIN (spp->point, last_sinp->from);
			  last_sinp->from = last_from;
			}
		  } else if (last_slp->choice == SEQLOC_PNT) {
		    /* change previous entry from point to interval that includes point */
		    last_spp = (SeqPntPtr)last_slp->data.ptrvalue;
			if (last_spp != NULL) {
              last_sinp = SeqIntNew ();
			  last_from = MIN (spp->point, last_spp->point);
			  last_sinp->from = last_from;
			  last_to = MAX (last_spp->point, spp->point);
			  last_sinp->to = last_to;
			  last_sinp->strand = spp->strand;
			  last_sinp->id = SeqIdDup (spp->id);
			  last_slp->choice = SEQLOC_INT;
			  last_slp->data.ptrvalue = last_sinp;
			  SeqPntFree (spp);
			}
		  }
		} else {
		  /* add new point */
          last_spp = SeqPntNew ();
          last_spp->point = spp->point;
		  last_to = spp->point;
		  last_from = spp->point;
          last_spp->strand = spp->strand;
		  last_strand = spp->strand;
          last_spp->id = SeqIdDup (spp->id);
          last_slp = ValNodeAddPointer (&head, SEQLOC_PNT, (Pointer) last_spp);
		}
	  }
	  break;
	default:
	  break;
	}
    slp = SeqLocFindNext (location, slp);
  }

  if (head == NULL) return NULL;
  if (head->next == NULL) return head;

  slp = ValNodeNew (NULL);
  slp->choice = SEQLOC_MIX;
  slp->data.ptrvalue = (Pointer) head;
  slp->next = NULL;

  return slp;
}

/* We need to be able to propagate features with locations on multiple segments
 * whose IDs may be found in separate alignments.
 * To do this, we will take the location of the master feature and for each sublocation,
 * we will find the sequence ID of the sublocation.  We will find the alignment that
 * contains the master feature sublocation sequence ID and construct a propagated 
 * sublocation for each sequence in the same alignment as the master feature sublocation 
 * sequence ID and add this to a list.
 * Once we have created all of the propagated sublocations, we will reconstitute the
 * total locations for the propagated features by associating each propagated sublocation
 * with other propagated sublocations with the same sequence ID or having a sequence ID
 * that belongs to the same segmented set.
 * We can then use the list of total locations to construct the propagated features.
 */

static SeqLocPtr MapSubLoc 
(SeqLocPtr   master_subloc,
 SeqAlignPtr salp, 
 Int4        master_row,
 Int4        prop_row,
 Boolean     gapSplit) 
{
  SeqLocPtr   prop_loc = NULL;
  SeqIdPtr    prop_sip;
  SeqIntPtr   sinp;
  IvalInfoPtr ival_head, ival;
  SeqPntPtr   spp;
  Uint1       strand;
  BioseqPtr   prop_bsp;
   
  if (master_subloc == NULL || salp == NULL || master_row < 1 || prop_row < 1)
  {
    return NULL;
  }
  
  prop_sip = AlnMgr2GetNthSeqIdPtr (salp, prop_row);
  if (prop_sip == NULL) return NULL;
  prop_bsp = BioseqFind (prop_sip);
  if (prop_bsp != NULL)
  {
    prop_sip = SeqIdFindBest (prop_bsp->id, SEQID_GENBANK);
    prop_sip = SeqIdDup (prop_sip);
  }
  
  switch (master_subloc->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) master_subloc->data.ptrvalue;
      if (sinp != NULL) {
        strand = sinp->strand;
        ival_head = GetAlignmentIntervals (salp, master_row, prop_row, sinp->from, sinp->to);
        ival_head = MergeAdjacentIntervals (ival_head);
        if (ival_head != NULL) {

          /* what if one or the other interval maps into a gap? */

          if (gapSplit) {
            for (ival = ival_head; ival != NULL; ival = ival->next) {
              sinp = SeqIntNew ();
              sinp->from = ival->start2;
              sinp->to = ival->stop2;
              sinp->strand = strand;
              sinp->id = prop_sip;
              prop_sip = NULL;
              ValNodeAddPointer (&prop_loc, SEQLOC_INT, (Pointer) sinp);
            }
          } else {
            sinp = SeqIntNew ();
            sinp->from = ival_head->start2;
            for (ival = ival_head; ival->next != NULL; ival = ival->next) continue;
            sinp->to = ival->stop2;
            sinp->strand = strand;
            sinp->id = prop_sip;
            prop_sip = NULL;
            ValNodeAddPointer (&prop_loc, SEQLOC_INT, (Pointer) sinp);
          }

        }
        IvalInfoFree (ival_head);
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) master_subloc->data.ptrvalue;
      if (spp != NULL) {
        strand = spp->strand;
        ival_head = GetAlignmentIntervals (salp, master_row, prop_row, spp->point, spp->point);
        if (ival_head != NULL) {

          spp = SeqPntNew ();
          spp->point = ival_head->start2;
          spp->strand = strand;
          spp->id = prop_sip;
          prop_sip = NULL;
          ValNodeAddPointer (&prop_loc, SEQLOC_PNT, (Pointer) spp);

        }
        IvalInfoFree (ival_head);
      }
      break;
    case SEQLOC_PACKED_PNT :
      /* not yet implemented */
      break;
    default :
      break;
  }
  prop_sip = SeqIdFree (prop_sip);
  
  return prop_loc;
}

typedef struct loclist
{
  SeqIdPtr sip_list;
  SeqLocPtr slp;
} LocListData, PNTR LocListPtr;

static LocListPtr FindLocListForSeqId (ValNodePtr loc_list, SeqIdPtr sip)
{
  SeqIdPtr   search_sip, find_sip;
  LocListPtr llp;
  ValNodePtr vnp;
  
  for (vnp = loc_list; vnp != NULL; vnp = vnp->next)
  {
    llp = (LocListPtr) vnp->data.ptrvalue;
    if (llp != NULL && llp->slp != NULL)
    {
      for (search_sip = llp->sip_list; search_sip != NULL; search_sip = search_sip->next)
      {
        for (find_sip = sip; find_sip != NULL; find_sip = find_sip->next)
        {
          if (SeqIdComp (find_sip, search_sip) == SIC_YES)
          {
            return llp;
          } 
        }
      }
    }
  }
  return NULL;
}

/* Compare the Sequence ID of this location with the sequence IDs of other locations in the
 * list.  If the sequence ID matches that of another loclist, add to that loclist,
 * otherwise add a new loclist entry with the sequence ID of this sequence and the sequence IDs
 * of other sequences in the same segset if the segset flag is set.
 */ 
static ValNodePtr AssociatePropagatedSubloc (ValNodePtr subloc_list, SeqLocPtr prop_loc, Boolean segset)
{
  ValNodePtr   vnp;
  LocListPtr   llp = NULL;
  SeqIdPtr     prop_sip;
  BioseqPtr    prop_bsp;
  BioseqSetPtr parent_set;
  SeqEntryPtr  sep;
  
  if (prop_loc == NULL) return subloc_list;
  prop_sip = SeqLocId (prop_loc);
  if (prop_sip == NULL) 
  {
    SeqLocFree (prop_loc);
    return subloc_list;
  }
  
  llp = FindLocListForSeqId (subloc_list, prop_sip);
  if (llp != NULL)
  {
    ValNodeAddPointer (&(llp->slp), prop_loc->choice, prop_loc->data.ptrvalue);
    prop_loc = ValNodeFree (prop_loc);
  }
  else
  {
    llp = (LocListPtr) MemNew (sizeof (LocListData));
    if (llp != NULL)
    {
      llp->sip_list = SeqIdDup (prop_sip);
      llp->slp = prop_loc;
      if (segset)
      {
        prop_bsp = BioseqFind (prop_sip);
        if (prop_bsp != NULL && prop_bsp->idx.parenttype == OBJ_BIOSEQSET)
        {
          parent_set = prop_bsp->idx.parentptr;
          if (parent_set != NULL && parent_set->_class == BioseqseqSet_class_parts)
          {
            /* add other IDs from set to list */
            for (sep = parent_set->seq_set; sep != NULL; sep = sep->next)
            {
              if (IS_Bioseq (sep))
              {
                prop_bsp = (BioseqPtr) sep->data.ptrvalue;
                prop_sip = SeqIdDup (prop_bsp->id);
                ValNodeAddPointer (&(llp->sip_list), prop_sip->choice, prop_sip->data.ptrvalue);
                prop_sip = ValNodeFree (prop_sip);
              }
            }
          }
        }
      }
      vnp = ValNodeAddPointer (&subloc_list, 0, llp);
    }
  }
  return subloc_list;
}

static ValNodePtr FreeLocList (ValNodePtr loc_list)
{
  LocListPtr llp;
  
  if (loc_list == NULL) return NULL;
  loc_list->next = FreeLocList (loc_list->next);
  llp = loc_list->data.ptrvalue;
  if (llp != NULL)
  {
    SeqIdFree (llp->sip_list);
    SeqLocFree (llp->slp);
    MemFree (llp);
  }
  ValNodeFree (loc_list);
  return NULL;
}

static Int4 GetMasterRow (SeqAlignPtr salp, SeqIdPtr sip)
{
  Int4     master_row = -1;
  SeqIdPtr sip_next;
  
  if (salp == NULL || sip == NULL)
  {
    return -1;
  }
  
  while (sip != NULL && master_row == -1)
  {
    sip_next = sip->next;
    sip->next = NULL;
    master_row = AlnMgr2GetFirstNForSip (salp, sip);
    sip->next = sip_next;
    sip = sip_next;
  }
  return master_row;
}

static Boolean 
NthAlignmentSequenceInSeqPropList
(SeqAlignPtr salp,
 Int4        n,
 ValNodePtr  seq_for_prop)
{
  SeqIdPtr    sip;
  ValNodePtr  seq_vnp;
  BioseqPtr   bsp;
  SeqEntryPtr sep;
  BioseqSetPtr bssp;
  
  if (salp == NULL || seq_for_prop == NULL)
  {
    return FALSE;
  }
  
  sip = AlnMgr2GetNthSeqIdPtr(salp, n);
  for (seq_vnp = seq_for_prop; seq_vnp != NULL; seq_vnp = seq_vnp->next)
  {
    if (SeqIdIn (sip, (SeqIdPtr) seq_vnp->data.ptrvalue))
    {
      return TRUE;
    }
    /* check for segments */
    bsp = BioseqFind (seq_vnp->data.ptrvalue);
    if (bsp != NULL && bsp->repr == Seq_repr_seg)
    {
      sep = SeqMgrGetSeqEntryForData (bsp);
      sep = sep->next; /* parts should be next */
      if (sep != NULL && IS_Bioseq_set (sep))
      {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts)
        {
          for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
          {
            if (IS_Bioseq (sep))
            {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              if (bsp != NULL && SeqIdIn (sip, bsp->id))
              {
                return TRUE;
              }
            }
          }
        }
      }
    }
  }
  
  
  return FALSE;
}

static ValNodePtr GetPropagatedSublocations 
(SeqLocPtr master_location,
 Boolean gap_split,
 ValNodePtr prop_loc_list,
 ValNodePtr seq_for_prop,
 BoolPtr warned_about_master)
{
  SeqLocPtr   master_subloc;
  SeqIdPtr    master_subloc_id;
  BioseqPtr   master_subloc_bsp;
  SeqLocPtr   tmp_loc;
  SeqAlignPtr salp;
  Int4        master_row, num_rows, prop_row; 
  SeqLocPtr   prop_loc;

  if (master_location == NULL) return prop_loc_list;
  
  master_subloc = SeqLocFindNext (master_location, NULL);
  while (master_subloc != NULL) 
  {
    master_subloc_id = SeqLocId (master_subloc);
    master_subloc_bsp = BioseqFind (master_subloc_id);
    master_subloc_id = master_subloc_bsp->id;
    if (master_subloc_bsp != NULL)
    {
      if (master_subloc_bsp->repr == Seq_repr_seg)
      {
        if (warned_about_master != NULL && ! *warned_about_master)
        {
          Message (MSG_OK, "Warning - you are propagating a feature that "
                  "contains locations on the master sequence.  These locations"
                  " will be mapped to the segments in the propagated features.");
          *warned_about_master = TRUE;
        }
        /* this is a location on the master segment */
        tmp_loc = SegLocToParts (master_subloc_bsp, master_subloc);
        prop_loc_list = GetPropagatedSublocations (tmp_loc, gap_split, prop_loc_list, 
                                                   seq_for_prop, warned_about_master);
        tmp_loc = SeqLocFree (tmp_loc);
      }
      else
      {
        salp = FindAlignmentsForBioseq (master_subloc_bsp);
        while (salp != NULL)
        {
          master_row = GetMasterRow (salp, master_subloc_id);
          num_rows = AlnMgr2GetNumRows (salp);
          for (prop_row = 1; prop_row <= num_rows; prop_row++)
          {
            if (prop_row == master_row) continue;
            if (! NthAlignmentSequenceInSeqPropList (salp, prop_row, seq_for_prop))
            {
              continue;
            }

            prop_loc = MapSubLoc (master_subloc, salp, master_row, prop_row, gap_split);
            prop_loc_list = AssociatePropagatedSubloc (prop_loc_list, prop_loc, TRUE);           
          }
          salp = salp->next;
        }
      }
    }
    master_subloc = SeqLocFindNext (master_location, master_subloc);
  }
  return prop_loc_list;
}

static ValNodePtr MapLocForProp
(SeqLocPtr   master_location,
 Boolean     gapSplit,
 ValNodePtr  seq_for_prop,
 BoolPtr     warned_about_master)
{
  ValNodePtr  prop_loc_list = NULL, vnp;
  LocListPtr  llp;
  SeqLocPtr   slp;
  
  if (master_location == NULL) 
  {
    return NULL;
  }

  prop_loc_list = GetPropagatedSublocations (master_location, gapSplit, prop_loc_list, 
                                             seq_for_prop, warned_about_master);
  
  /* now fix locations in prop_loc_list (add SEQLOC_MIX header to the locations that
   * are mixed)
   */
  for (vnp = prop_loc_list; vnp != NULL; vnp = vnp->next)
  {
    llp = (LocListPtr) vnp->data.ptrvalue;
    if (llp != NULL)
    {
      if (llp->slp != NULL)
      {
        if (llp->slp->next != NULL)
        {
          slp = ValNodeNew (NULL);
          slp->choice = SEQLOC_MIX;
          slp->data.ptrvalue = llp->slp;
          llp->slp = slp;
        }
      }
    }
  }
  return prop_loc_list;
}


static SeqIdPtr GetSegIdList (BioseqPtr master_seg)
{
  SeqIdPtr  sip_list = NULL, sip_last = NULL, sip;
  SeqLocPtr slp;
  if (master_seg == NULL)
  {
    return NULL;
  }
  if (master_seg->repr != Seq_repr_seg)
  {
    sip_list = SeqIdDupList (master_seg->id);
  }
  else
  {    
    for (slp = master_seg->seq_ext; slp != NULL; slp = slp->next)
    {
      sip = SeqIdDup (SeqLocId (slp));
      if (sip_last == NULL)
      {
        sip_list = sip;
      }
      else
      {
        sip_last->next = sip;
      }
      sip_last = sip;
    }
  }
  return sip_list;
}

static void PropagateCodeBreaks 
(CdRegionPtr crp,
 SeqIdPtr sip,
 ValNodePtr codebreak_location_list,
 ValNodePtr codebreak_choice_list)
{
  CodeBreakPtr cbp, last_cbp = NULL;
  ValNodePtr   choice_vnp, cbp_vnp;
  LocListPtr   llp;
  BioseqPtr    cds_bsp;
  SeqIdPtr     sip_list;
  
  if (crp == NULL || codebreak_location_list == NULL || crp->code_break == NULL)
  {
    return;
  }
  
  crp->code_break = CodeBreakFree (crp->code_break);
  
  cds_bsp = BioseqFind (sip);
  if (cds_bsp == NULL)
  {
    return;
  }

  sip_list = GetSegIdList (cds_bsp);
  
  for (cbp_vnp = codebreak_location_list, choice_vnp = codebreak_choice_list;
       cbp_vnp != NULL && choice_vnp != NULL;
       cbp_vnp = cbp_vnp->next, choice_vnp = choice_vnp->next)
  {
    llp = FindLocListForSeqId (cbp_vnp->data.ptrvalue, sip_list);
    if (llp != NULL)
    {
      cbp = CodeBreakNew ();
      if (cbp != NULL)
      {
        cbp->loc = llp->slp;
        llp->slp = NULL;
        MemCpy (&(cbp->aa), choice_vnp->data.ptrvalue, sizeof (cbp->aa));
        if (last_cbp == NULL)
        {
          crp->code_break = cbp;
        }
        else
        {
          last_cbp->next = cbp;
        }
      }
    }
  }
  SeqIdFree (sip_list);
}

static void PropagateAnticodons 
(tRNAPtr trp,
 SeqIdPtr sip,
 ValNodePtr anticodon_location_list)
{
  LocListPtr llp;
  BioseqPtr  trna_bsp;
  SeqIdPtr   sip_list;
  
  if (trp == NULL || anticodon_location_list == NULL || trp->anticodon == NULL)
  {
    return;
  }
  
  trp->anticodon = SeqLocFree (trp->anticodon);
  
  trna_bsp = BioseqFind (sip);
  if (trna_bsp == NULL)
  {
    return;
  }
  sip_list = GetSegIdList (trna_bsp);
  
  llp = FindLocListForSeqId (anticodon_location_list, sip_list);
  if (llp != NULL)
  {
    trp->anticodon = llp->slp;
    llp->slp = NULL;
  }
  
  SeqIdFree (sip_list);
}


static void ExtendLocToEnd (
  SeqLocPtr location,
  BioseqPtr bsp,
  Uint1 strand
)

{
  SeqIntPtr  sinp;
  SeqLocPtr  slp, last = NULL;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    last = slp;
    slp = SeqLocFindNext (location, slp);
  }
  if (last == NULL) return;

  switch (last->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) last->data.ptrvalue;
      if (sinp != NULL) {
        if (strand == Seq_strand_minus) {
          sinp->from = 0;
        } else {
          sinp->to = bsp->length - 1;
        }
      }
    case SEQLOC_PNT :
      /* not yet implemented */
      break;
    case SEQLOC_PACKED_PNT :
      /* not yet implemented */
      break;
    default :
      break;
  }
}

static void TruncateCDS (
  SeqFeatPtr sfp,
  Uint1 frame,
  BioseqPtr pbsp
)

{
  Int4       len;
  SeqIntPtr  sinp;
  SeqLocPtr  slp;
  Int4       total = 0;

  if (frame > 0) {
    frame--;
  }
  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) {
    len = SeqLocLen (slp);

    if (len + total - frame <= (pbsp->length + 1) * 3) {
      total += len;
    } else {
      if (slp->choice == SEQLOC_INT) {
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          len = (pbsp->length + 1) * 3 - total;
          if (sinp->strand == Seq_strand_minus) {
            sinp->from = sinp->to - len + 1;
          } else {
            sinp->to = sinp->from + len - 1;
          }
        }
      }
      return;
    }

    slp = SeqLocFindNext (sfp->location, slp);
  }
}

/*------------------------------------------------------------------*/
/*                                                                  */
/* PropagateCDS () - Called from DoFeatProp() for CDS-specific      */
/*                   feature propagation.                           */
/*                                                                  */
/*------------------------------------------------------------------*/

static void PropagateCDS (SeqFeatPtr dup,
			  ProtRefPtr prp,
			  BioseqPtr  newbsp,
			  Boolean    stopCDS,
			  Boolean    transPast,
			  Boolean    cds3end,
			  Uint1      frame,
			  Uint1      strand)
{
  Uint2           entityID;
  MolInfoPtr      mip;
  Boolean         partial3;
  Boolean         partial5;
  BioseqPtr       pbsp;
  SeqDescrPtr     sdp;
  SeqIdPtr        sip;
  SeqFeatXrefPtr  xref;
  Boolean         xtend;

  /* Check parameters */

  if (dup == NULL || dup->data.choice != SEQFEAT_CDREGION || prp == NULL)
    return;

  /* Extend the location to the end if that was checked */

  if (transPast && cds3end) 
    ExtendLocToEnd (dup->location, newbsp, strand);

  /**/

  prp = AsnIoMemCopy ((Pointer) prp,
                       (AsnReadFunc) ProtRefAsnRead,
                       (AsnWriteFunc) ProtRefAsnWrite);

  xref = SeqFeatXrefNew ();
  if (xref == NULL)
    return;
  xref->data.choice = SEQFEAT_PROT;
  xref->data.value.ptrvalue = (Pointer) prp;
  xref->next = dup->xref;
  dup->xref = xref;

  entityID = ObjMgrGetEntityIDForPointer (newbsp);
  PromoteXrefsEx (dup, newbsp, entityID, (Boolean) (! stopCDS), FALSE, FALSE);

  /* Truncate new CDS based on new protein length */

  sip = SeqLocId (dup->product);
  if (sip == NULL)
    return;

  pbsp = BioseqFindCore (sip);
  if (pbsp == NULL)
    return;

  TruncateCDS (dup, frame, pbsp);

  /**/

  CheckSeqLocForPartial (dup->location, &partial5, &partial3);
  if (cds3end) {
    xtend = FALSE;
    if (strand == Seq_strand_minus) {
      if (SeqLocStop (dup->location) == 0) {
        xtend = TRUE;
      }
    } else {
      if (SeqLocStop (dup->location) == newbsp->length - 1) {
        xtend = TRUE;
      }
    }
    if (xtend) {
      partial3 = TRUE;
      SetSeqLocPartial (dup->location, partial5, partial3);
      for (sdp = pbsp->descr; sdp != NULL; sdp = sdp->next) {
        if (sdp->choice != Seq_descr_molinfo) continue;
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip == NULL) continue;
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        }
      }
    }
  }

  /* Set partial flag */

  dup->partial = (Boolean) (partial5 || partial3);
}


static SeqAlignPtr CheckForProteinAlignment (ByteStorePtr bs, BioseqPtr match_prot)
{
  BioseqPtr   newBsp;
  SeqAlignPtr salp;
  SeqEntryPtr nwsep;
  Boolean     revcomp;
  ErrSev      oldsev;
  
  if (bs == NULL || match_prot == NULL)
  {
    return NULL;
  }
  newBsp = BioseqNew ();
  if (newBsp == NULL) {
    return NULL;
  }

  newBsp->id = SeqIdParse ("lcl|ProtAlign");
  newBsp->repr = Seq_repr_raw;
  newBsp->mol = Seq_mol_aa;
  newBsp->seq_data_type = Seq_code_ncbieaa;
  newBsp->seq_data = bs;
  newBsp->length = BSLen (bs);

  /* create SeqEntry for temporary protein bioseq to live in */
  nwsep = SeqEntryNew ();
  nwsep->choice = 1;
  nwsep->data.ptrvalue = newBsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) newBsp, nwsep);

  oldsev = ErrSetMessageLevel (SEV_MAX);
  salp = Sequin_GlobalAlign2Seq (match_prot, newBsp, &revcomp);
  ErrSetMessageLevel (oldsev);

  SeqMgrDeleteFromBioseqIndex (newBsp);
  /* the calling function will free bs */
  newBsp->seq_data = NULL;
  nwsep = SeqEntryFree (nwsep);
  return salp;
}

/*------------------------------------------------------------------*/
/*                                                                  */
/* CalculateReadingFrame () -- Calculates a sequence's reading      */
/*                             frame by seeing which frame          */
/*                             generates the most amino acids       */
/*                             when converted to a protein.         */
/*                                                                  */
/*------------------------------------------------------------------*/

static Int2 CalculateReadingFrame (SeqFeatPtr sfp, Boolean partial3, BioseqPtr match_prot)
{
  ByteStorePtr  bs;
  CdRegionPtr   crp;
  Int4          len;
  Int4          max;
  Uint1         frame;
  Int2          i;
  CharPtr       protstr;
  SeqAlignPtr   salp = NULL;
  Boolean       best_is_aligned = FALSE;

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;

  max = 0;
  frame = 0;

  if (! partial3 && crp->frame != 0) 
  {
  	bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  	if (bs != NULL)
  	{
      protstr = BSMerge (bs, NULL);
      BSFree (bs);
      if (protstr != NULL) {
        len = StringLen (protstr);
        if (len > 0 && protstr [len - 1] == '*') {
          MemFree (protstr);
          return crp->frame; 	
        }
        MemFree (protstr);
      }
  	}
  }
  for (i = 1; i <= 3; i++) {
    crp->frame = (Uint1) i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    salp = CheckForProteinAlignment (bs, match_prot);
    len = BSLen (bs);
    BSFree (bs);
    if (salp == NULL) 
    {
      /* if no alignment, take longest protein */
      if (! best_is_aligned && len > max) {
        max = len;
        frame = (Uint1) i;
      }
    } else if (!best_is_aligned) {
      /* if this is the first alignment encountered, use this */
      max = len;
      frame = (Uint1) i;
      best_is_aligned = TRUE;
    } else if (len > max) {
      /* if this is not the first alignment, but the sequence is longer, use this */
      max = len;
      frame = (Uint1) i;
    }
    salp = SeqAlignFree (salp);
  }

  return frame;
}

/* we don't want to propagate Prot features if the target
 * sequence already has one.
 */
static Boolean OkToPropagate(SeqFeatPtr sfp, BioseqPtr bsp)
{
  SeqMgrFeatContext fcontext;
  
  if (sfp == NULL || bsp == NULL) return FALSE;
  /* never ok if a gap */
  if (sfp->idx.subtype == FEATDEF_gap) return FALSE;
  /* always ok if not a Prot feature */
  if (sfp->idx.subtype != FEATDEF_PROT) return TRUE;
  /* always ok if not a protein sequence */
  if (ISA_na (bsp->mol)) return TRUE;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_PROT, &fcontext);
  if (sfp == NULL)
    return TRUE;
  else
    return FALSE;
}


static SeqIdPtr GetSegSetId (SeqLocPtr slp)
{
  SeqLocPtr    sub_slp;
  SeqIdPtr     master_sip, part_sip;
  BioseqPtr    seg_bsp;
  BioseqSetPtr parent_set;
  BioseqSetPtr last_parent_set = NULL;
  SeqEntryPtr  sep;
  BioseqPtr    parent_bsp;
  
  if (slp == NULL) return NULL;
  master_sip = SeqLocId (slp);
  if (master_sip != NULL || slp->choice != SEQLOC_MIX)
  {
    return master_sip;
  }
  /* make sure all parts are from the same segmented set */
  for (sub_slp = slp->data.ptrvalue; sub_slp != NULL; sub_slp = sub_slp->next)
  {
    part_sip = SeqLocId (sub_slp);
    seg_bsp = BioseqFind (part_sip);
    if (seg_bsp == NULL || seg_bsp->idx.parenttype != OBJ_BIOSEQSET 
        || seg_bsp->idx.parentptr == NULL)
    {
      return NULL;
    }
    parent_set = (BioseqSetPtr) seg_bsp->idx.parentptr;
    if (parent_set->_class != BioseqseqSet_class_parts 
        || parent_set->idx.parenttype != OBJ_BIOSEQSET
        || parent_set->idx.parentptr == NULL)
    {
      return NULL;
    }
    if (last_parent_set == NULL)
    {
      last_parent_set = parent_set;
    }
    else if (last_parent_set != parent_set)
    {
      return NULL;
    }
  }
  if (last_parent_set == NULL)
  {
    return NULL;
  }
  parent_set = (BioseqSetPtr) last_parent_set->idx.parentptr;
  if (parent_set->_class != BioseqseqSet_class_segset)
  {
    return NULL;
  }
  for (sep = parent_set->seq_set; sep != NULL && ! IS_Bioseq (sep); sep = sep->next)
  {
  }
  if (sep == NULL) return NULL;
  parent_bsp = sep->data.ptrvalue;
  if (parent_bsp == NULL) return NULL;
  master_sip = parent_bsp->id;
  return master_sip;  
}

/* for each feature, find all propagated locations and create features using
 * those locations.
 */
static void 
PropagateOneFeat
(SeqFeatPtr sfp,
 Boolean    gapSplit,
 Boolean    fuse_joints,
 Boolean    stopCDS,
 Boolean    transPast,
 Boolean    cds3end,
 ValNodePtr seq_for_prop,
 BoolPtr    warned_about_master)
{
  ValNodePtr      feature_location_list, vnp;
  ValNodePtr      codebreak_location_list = NULL;
  ValNodePtr      codebreak_choice_list = NULL;
  ValNodePtr      anticodon_location_list = NULL;
  CodeBreakPtr    cbp;
  CdRegionPtr     crp;
  SeqFeatPtr      dup;
  Uint1           frame = 0;
  BioseqPtr       newbsp;
  SeqLocPtr       newloc, mergedloc;
  Boolean         partial5;
  Boolean         partial3;
  RnaRefPtr       rrp;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  tRNAPtr         trp;
  LocListPtr      llp;
  Uint2           strand;
  ProtRefPtr      prp = NULL;
  BioseqPtr       pbsp = NULL;
  SeqFeatPtr      prot;
  
  if (sfp == NULL || sfp->location == NULL) return;
  
  feature_location_list = MapLocForProp (sfp->location, gapSplit, 
                                         seq_for_prop, warned_about_master);
  if (feature_location_list == NULL) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  
  /* also need to propagate locations for CDS code breaks and tRNA anticodons */
  if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL)
    {
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) 
      {
        vnp = MapLocForProp (cbp->loc, gapSplit, seq_for_prop, warned_about_master);
        ValNodeAddPointer (&codebreak_location_list, 0, vnp);
        ValNodeAddPointer (&codebreak_choice_list, 0, &(cbp->aa));
      }
    }
    sip = SeqLocId (sfp->product);
    if (sip != NULL)
    {
      pbsp = BioseqFindCore (sip);
      if (pbsp != NULL)
      {
        prot = SeqMgrGetBestProteinFeature (pbsp, NULL);
        if (prot != NULL && prot->data.choice == SEQFEAT_PROT)
        {
          prp = (ProtRefPtr) prot->data.value.ptrvalue;          
        }
      }
    }
  }
  else if (sfp->data.choice == SEQFEAT_RNA)
  {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 2) 
    {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL && trp->anticodon != NULL) 
      {
        anticodon_location_list = MapLocForProp (trp->anticodon, gapSplit, 
                                                 seq_for_prop, warned_about_master);
      }
    }
  }
  
  for (vnp = feature_location_list; vnp != NULL; vnp = vnp->next)
  {
    llp = (LocListPtr) vnp->data.ptrvalue;
    if (llp == NULL)
    {
      continue;
    }
    newloc = llp->slp;
    llp->slp = NULL;
    if (newloc == NULL)
    {
      continue;
    }

    mergedloc = NULL;
    if (fuse_joints) 
    {
      mergedloc = MergeAdjacentSeqLocIntervals (newloc);
	    if (mergedloc != NULL) 
	    {
        SeqLocFree (newloc);
        newloc = mergedloc;
	    }
    }
    
    /* if we have a mixed location with multiple sequences, SeqLocId will
     * return NULL.
     * We should check to see if we are trying to create a feature for
     * a segset, in which case we'll want the Bioseq for the master segment.
     */
    sip = GetSegSetId (newloc);

    dup = AsnIoMemCopy ((Pointer) sfp,
                        (AsnReadFunc) SeqFeatAsnRead,
                        (AsnWriteFunc) SeqFeatAsnWrite);
    SeqLocFree (dup->location);
    dup->location = newloc;
    SetSeqLocPartial (dup->location, partial5, partial3);
    
    /* do not propagate feature IDs or links */
    ClearFeatIDs (dup);
    ClearFeatIDXrefs (dup);

    /* clean up product before we look for (and maybe don't find) the newbsp */
    if (dup->product != NULL)
      dup->product = SeqLocFree (dup->product);

    switch (dup->data.choice) {
      case SEQFEAT_CDREGION :
        crp = (CdRegionPtr) dup->data.value.ptrvalue;
        if (crp != NULL) {
          crp->frame = CalculateReadingFrame (dup, partial3, pbsp);
          frame = crp->frame;
 
          PropagateCodeBreaks (crp, sip,
                               codebreak_location_list,         
                               codebreak_choice_list);         
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) dup->data.value.ptrvalue;
        if (rrp != NULL && rrp->ext.choice == 2) {
          trp = (tRNAPtr) rrp->ext.value.ptrvalue;
          if (trp != NULL && trp->anticodon != NULL) 
          {
            PropagateAnticodons (trp, sip, anticodon_location_list);
          }
        }
        break;
      default :
        break;
    }

    newbsp = BioseqFindCore (sip);
    if (newbsp == NULL)
      return;
  
    /* need to call OkToPropagate with sfp instead of dup
     * because dup has not been indexed yet, so subtype isn't set.
     */
    if (OkToPropagate(sfp, newbsp))
    {
      sep = SeqMgrGetSeqEntryForData (newbsp);
      if (sep == NULL)
        return;
      CreateNewFeature (sep, NULL, dup->data.choice, dup);

      /* If we're doing a CDS propagation, then */
      /* do the extra stuff related to that.    */

      if (SEQFEAT_CDREGION == dup->data.choice)
      {
        strand = SeqLocStrand (dup->location);

        PropagateCDS (dup, prp, newbsp, stopCDS, transPast, cds3end, frame, strand);        
      }
    }
    else
    {
      SeqFeatFree (dup);
    }
  }  
  
  feature_location_list = FreeLocList (feature_location_list);
  anticodon_location_list = FreeLocList (anticodon_location_list);
  for (vnp = codebreak_location_list; vnp != NULL; vnp = vnp->next)
  {
    vnp->data.ptrvalue = FreeLocList (vnp->data.ptrvalue);
  }
  codebreak_location_list = ValNodeFree (codebreak_location_list);
  codebreak_choice_list = ValNodeFree (codebreak_choice_list);
  
}


static Boolean CDSgoesToEnd (
  BioseqPtr bsp,
  SeqMgrFeatContext PNTR fcontext
)

{
  if (fcontext->strand == Seq_strand_minus) {
    if (fcontext->left == 0 && fcontext->partialR) return TRUE;
  } else {
    if (fcontext->right == bsp->length - 1 && fcontext->partialR) return TRUE;
  }
  return FALSE;
}

extern void DoFixCDS (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BaseFormPtr        bfp;
  ByteStorePtr       bs;
  BioseqPtr          bsp;
  Boolean            change_partials = FALSE;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  size_t             len;
  Boolean            partial5;
  Boolean            partial3;
  SeqIntPtr          sintp;
  SeqLocPtr          slp;
  CharPtr            str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  bfp = (BaseFormPtr) userdata;
  if (SeqMgrGetDesiredFeature (bfp->input_entityID, NULL,
                               0, 0, sfp, &context) != sfp) return;
  bsp = context.bsp;
  if (bsp == NULL) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp->frame > 1) {
    if (context.strand == Seq_strand_minus) {
      if (context.right == bsp->length - 1) {
        partial5 = TRUE;
        change_partials = TRUE;
      }
    } else {
      if (context.left == 0) {
        partial5 = TRUE;
        change_partials = TRUE;
      }
    }
  }
    bs = ProteinFromCdRegion (sfp, TRUE);
    if (bs != NULL) {
      str = BSMerge (bs, NULL);
      BSFree (bs);
      if (str != NULL) {
        if (*str == '-') {
          if (! partial5) {
            partial5 = TRUE;
            change_partials = TRUE;
          }
        }
        len = StringLen (str);
        if (len > 0 && str [len - 1] != '*') {
          if (context.strand == Seq_strand_minus) {
          } else {
            if (bsp->length - context.right < 3) {
              slp = SeqLocFindNext (sfp->location, NULL);
              while (slp != NULL) {
                if (slp->choice == SEQLOC_INT) {
                  sintp = (SeqIntPtr) slp->data.ptrvalue;
                  if (sintp != NULL) {
                    if (sintp->to == context.right) {
                      sintp->to = bsp->length - 1;
                    }
                  }
                }
                slp = SeqLocFindNext (sfp->location, slp);
              }
            }
            partial3 = TRUE;
            change_partials = TRUE;
          }
        }
        MemFree (str);
      }
    }
  if (change_partials) {
    SetSeqLocPartial (sfp->location, partial5, partial3);
    ResynchCDSPartials (sfp, NULL);
  }
}

extern SeqFeatPtr 
GetNextFeatureOnSegOrMaster 
(BioseqPtr bsp, SeqFeatPtr sfp, Uint4 itemID, Uint4 index, SeqMgrFeatContextPtr fcontext)
{
  BioseqSetPtr       bssp;
  SeqEntryPtr        sep;
  BioseqPtr          master_bsp = NULL;
  SeqFeatPtr         next_sfp;
  SeqLocPtr          slp;
  SeqIdPtr           loc_id;
  Boolean            on_this_segment = FALSE;
  
  if (bsp == NULL)
  {
    return NULL;
  }
  if (bsp->idx.parenttype != OBJ_BIOSEQSET || bsp->idx.parentptr == NULL)
  {
    return SeqMgrGetNextFeature (bsp, sfp, itemID, index, fcontext);
  }

  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_parts
      || bssp->idx.parenttype != OBJ_BIOSEQSET || bsp->idx.parentptr == NULL)
  {
    return SeqMgrGetNextFeature (bsp, sfp, itemID, index, fcontext);
  }
  
  bssp = bssp->idx.parentptr;
  if (bssp->_class != BioseqseqSet_class_segset)
  {
    return SeqMgrGetNextFeature (bsp, sfp, itemID, index, fcontext);
  }

  for (sep = bssp->seq_set; sep != NULL && master_bsp == NULL; sep = sep->next)
  {
    if (IS_Bioseq (sep))
    {
      master_bsp = sep->data.ptrvalue;
      if (master_bsp != NULL && master_bsp->repr != Seq_repr_seg)
      {
        master_bsp = NULL;
      }
    }
  }
  
  if (master_bsp == NULL)
  {
    return SeqMgrGetNextFeature (bsp, sfp, itemID, index, fcontext);
  }
  
  next_sfp = SeqMgrGetNextFeature (master_bsp, sfp, itemID, index, fcontext);
  if (next_sfp == NULL) return NULL;

  while (next_sfp != NULL && !on_this_segment)
  {
    for (slp = SeqLocFindNext (next_sfp->location, NULL);
         slp != NULL && ! on_this_segment; slp = SeqLocFindNext (next_sfp->location, slp))
    {
      loc_id = SeqLocId (slp);
      if (SeqIdIn (loc_id, bsp->id))
      {
        on_this_segment = TRUE;
      }
    }
    if (!on_this_segment)
    {
      next_sfp = SeqMgrGetNextFeature (master_bsp, next_sfp, itemID, index, fcontext);
    }
  }
  
  return next_sfp;  
}


static void AcceptFeatProp (
  ButtoN b
)

{
  BioseqPtr          bsp;
  Boolean            cds3end;
  DenseSegPtr        dsp;
  Uint2              entityID;
  SeqMgrFeatContext  fcontext;
  FprDataPtr         fdp;
  Boolean            fixCDS;
  Boolean            gapSplit;
  SeqAlignPtr        salp;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  SeqIdPtr           sip_head;
  SeqIdPtr           sip_new;
  SeqIdPtr           sip_prev;
  SeqIdPtr           sip_tmp;
  Boolean            stopCDS;
  Boolean            transPast;
  Boolean            fuse_joints = FALSE;
  Boolean            warned_about_master = FALSE;
  ValNodePtr         seq_for_prop = NULL;

  fdp = (FprDataPtr) GetObjectExtra (b);
  if (fdp == NULL) return;
  SafeHide (fdp->form);

  bsp = fdp->bsp;
  salp = fdp->salp;
  if (bsp == NULL || salp == NULL) {
    Remove (fdp->form);
    return;
  }

  if (GetValue (fdp->allOrSel) == 1) {
    fdp->selFeatItemID = 0;
  }
  if (GetValue (fdp->gapSplit) == 1) {
    gapSplit = FALSE;
  } else {
    gapSplit= TRUE;
  }
  if (GetStatus (fdp->stopCDS)) {
    stopCDS = TRUE;
  } else {
    stopCDS = FALSE;
  }
  if (GetStatus (fdp->transPast)) {
    transPast = TRUE;
  } else {
    transPast = FALSE;
  }
  if (GetStatus (fdp->fixCDS)) {
    fixCDS = TRUE;
  } else {
    fixCDS = FALSE;
  }
  if (GetStatus (fdp->fuseJoints)) {
    fuse_joints = TRUE;
  } else {
	  fuse_joints = FALSE;
  }
  
  seq_for_prop = DialogToPointer (fdp->sequence_list_dlg);
  if (seq_for_prop != NULL && seq_for_prop->next == NULL)
  {
    for (sip = bsp->id; sip != NULL; sip = sip->next)
    {
      if (SeqIdIn (sip, (SeqIdPtr) seq_for_prop->data.ptrvalue))
      {
        Message (MSG_ERROR, "Can't propagate from a sequence to itself!");
        SafeShow (fdp->form);
        return;
      }
    }
  }

  SeqEntrySetScope (NULL);

  /* need to find alignment for each feature and row within that alignment for the feature */

  dsp = (DenseSegPtr)(salp->segs);
  sip = SeqIdFindBest (bsp->id, 0);
  sip_tmp = dsp->ids;
  sip_head = sip_prev = NULL;
  while (sip_tmp != NULL)
  {
    if (SeqIdComp(sip_tmp, bsp->id) == SIC_YES)
       sip_new = SeqIdDup(sip);
    else
       sip_new = SeqIdDup(sip_tmp);
    if (sip_head != NULL)
    {
       sip_prev->next = sip_new;
       sip_prev = sip_new;
    } else
       sip_head = sip_prev = sip_new;
    sip_tmp = sip_tmp->next;
  }
  dsp->ids = sip_head;

  if (fdp->selFeatItemID != 0) {

    /* propagate single selected feature */

    sfp = SeqMgrGetDesiredFeature (0, bsp, fdp->selFeatItemID, 0, NULL, &fcontext);
    if (sfp != NULL) {
      cds3end = CDSgoesToEnd (bsp, &fcontext);
      PropagateOneFeat (sfp, gapSplit, fuse_joints, stopCDS, transPast, 
                        cds3end, seq_for_prop, &warned_about_master);
    }
  } else {

    /* propagate all features on bioseq */

    sfp = GetNextFeatureOnSegOrMaster (bsp, NULL, 0, 0, &fcontext);
    while (sfp != NULL) {
      cds3end = CDSgoesToEnd (bsp, &fcontext);
      PropagateOneFeat (sfp, gapSplit, fuse_joints, stopCDS, transPast,
                        cds3end, seq_for_prop, &warned_about_master);

      sfp = GetNextFeatureOnSegOrMaster (bsp, sfp, 0, 0, &fcontext);
    }
  }

  seq_for_prop = ValNodeFree (seq_for_prop);
  if (fixCDS) {
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    sep = GetTopSeqEntryForEntityID (entityID);
    fdp->input_entityID = entityID;
    /* reindex before calling DoFixCDS */
    SeqMgrIndexFeatures (entityID, NULL);
    VisitFeaturesInSep (sep, fdp, DoFixCDS);
  }

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);

  Remove (fdp->form);
}

static void SetFeaturePropagateAccept (Pointer userdata)
{
  FprDataPtr  fdp;
  ValNodePtr  err_list;
  
  fdp = (FprDataPtr) userdata;
  if (fdp == NULL)
  {
    return;
  }
  err_list = TestDialog (fdp->sequence_list_dlg);
  if (err_list == NULL)
  {
    Enable (fdp->accept);
  }
  else
  {
    Disable (fdp->accept);
  }
  ValNodeFree (err_list);
}

static void FeaturePropagateFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr        bfp;
  StdEditorProcsPtr  sepp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp == NULL) return;
  switch (mssg) {
    case VIB_MSG_CLOSE :
      Remove (f);
      break;
    default :
      sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
      if (sepp != NULL && sepp->handleMessages != NULL) {
        sepp->handleMessages (f, mssg);
      }
      break;
  }
}

#ifdef WIN_MAC
extern void UpdateSequenceFormActivated (WindoW w)

{
  HANDLE openItem;
  HANDLE closeItem;
  BoolPtr pInitialFormsActive;
  
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  
  pInitialFormsActive = GetAppProperty (SEQFORM_INIT_ACTIVE);
  if (pInitialFormsActive != NULL)
  {
    *pInitialFormsActive = FALSE;
  }
  
  openItem = (HANDLE) GetAppProperty (SEQFORM_OPEN_ITEM);
  closeItem = (HANDLE) GetAppProperty (SEQFORM_CLOSE_ITEM);

  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   NULL);
}

extern void UpdateSequenceFormActivate (WindoW w)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    if (bfp->activate != NULL) {
      bfp->activate (w);
    }
  }
}
#endif


extern ForM FeaturePropagateForm (
  BioseqPtr bsp,
  SeqAlignPtr salp,
  Uint2 selFeatItemID
)

{
  ButtoN      b;
  GrouP       c;
  GrouP       seq_choice_grp;
  FprDataPtr  fdp;
  GrouP       g;
  PrompT      ppt;
  SeqIdPtr    sip;
  Char        strid [41];
  Char        txt [128];
  WindoW      w;
  Uint2       entityID;

  if (bsp == NULL) return NULL;
  fdp = (FprDataPtr) MemNew (sizeof (FprData));
  if (fdp == NULL) return NULL;
  w = FixedWindow (-50, -33, -10, -10, "Feature Propagate", NULL);
  if (w == NULL) return NULL;

  SetObjectExtra (w, (Pointer) fdp, StdCleanupFormProc);
  fdp->form = (ForM) w;
  fdp->formmessage = FeaturePropagateFormMessage;

#ifdef WIN_MAC
  fdp->activate = UpdateSequenceFormActivated;
  SetActivate (w, UpdateSequenceFormActivate);
#endif

  fdp->bsp = bsp;
  fdp->salp = salp;
  fdp->selFeatItemID = selFeatItemID;

  sip = SeqIdFindWorst (bsp->id);
  SeqIdWrite (sip, strid, PRINTID_REPORT, sizeof (strid) - 1);
  if (ISA_na (bsp->mol)) {
    sprintf (txt, "Propagate from %s to", strid);
  } else {
    sprintf (txt, "Propagate from %s to", strid);
  }

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 5, 5);

  seq_choice_grp = HiddenGroup (g, 0, 2, NULL);
  ppt = StaticPrompt (seq_choice_grp, txt, 0, 0, programFont, 'c');
  entityID = ObjMgrGetEntityIDForPointer (bsp);
  fdp->sequence_list_dlg = SequenceSelectionDialog (seq_choice_grp, 
                                                    SetFeaturePropagateAccept,
                                                    fdp, 
                                                    TRUE, 
                                                    ISA_na (bsp->mol), 
                                                    ISA_aa (bsp->mol), 
                                                    entityID);

  fdp->allOrSel = HiddenGroup (g, 2, 0, NULL);
  RadioButton (fdp->allOrSel, "All Features");
  b = RadioButton (fdp->allOrSel, "Selected Feature");
  if (selFeatItemID > 0) {
    SetValue (fdp->allOrSel, 2);
  } else {
    Disable (b);
    SetValue (fdp->allOrSel, 1);
  }

  fdp->gapSplit = HiddenGroup (g, 2, 0, NULL);
  RadioButton (fdp->gapSplit, "Extend over gaps");
  RadioButton (fdp->gapSplit, "Split at gaps");
  SetValue (fdp->gapSplit, 1);

  fdp->stopCDS = CheckBox (g, "Stop CDS translation at internal stop codon", NULL);
  SetStatus (fdp->stopCDS, FALSE);

  fdp->transPast = CheckBox (g, "Translate CDS after partial 3' boundary", NULL);

  fdp->fixCDS = CheckBox (g, "Cleanup CDS partials after propagation", NULL);
  SetStatus (fdp->fixCDS, TRUE);

  fdp->fuseJoints = CheckBox (g, "Fuse adjacent propagated intervals", NULL);
  SetStatus (fdp->fuseJoints, FALSE);

  c = HiddenGroup (w, 4, 0, NULL);
  fdp->accept = DefaultButton (c, "Accept", AcceptFeatProp);
  SetObjectExtra (fdp->accept, (Pointer) fdp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) seq_choice_grp, (HANDLE) fdp->allOrSel,
                (HANDLE) fdp->gapSplit, (HANDLE) fdp->stopCDS,
                (HANDLE) fdp->transPast, (HANDLE) fdp->fixCDS,
				(HANDLE) fdp->fuseJoints,
                (HANDLE) c, NULL);
  RealizeWindow (w);
  SendMessageToDialog (fdp->sequence_list_dlg, NUM_VIB_MSG + 1);


  return (ForM) w;
}

typedef struct batchapplyfeaturedetailsdlg
{
  DIALOG_MESSAGE_BLOCK
  
  TexT           defline;
  TexT           geneName;
  TexT           geneDesc;
  TexT           protName;
  TexT           protDesc;
  TexT           rnaName;
  TexT           featcomment;
  DialoG         featdef_choice_dlg;
  PopuP          rnaSubType;
  PopuP          reading_frame;
  Int4           feattype;
  
} BatchApplyFeatureDetailsDlgData, PNTR BatchApplyFeatureDetailsDlgPtr;

extern BatchApplyFeatureDetailsPtr 
BatchApplyFeatureDetailsFree 
(BatchApplyFeatureDetailsPtr bafdp)
{
  if (bafdp != NULL)
  {
    bafdp->defline = MemFree (bafdp->defline);
    bafdp->geneName = MemFree (bafdp->geneName);
	bafdp->geneDesc = MemFree (bafdp->geneDesc);
    bafdp->protName = MemFree (bafdp->protName);
    bafdp->protDesc = MemFree (bafdp->protDesc);
    bafdp->rnaName = MemFree (bafdp->rnaName);
    bafdp->featcomment = MemFree (bafdp->featcomment);
    bafdp->featdef_name = MemFree (bafdp->featdef_name);
    bafdp = MemFree (bafdp);
  }
  return bafdp;
}

extern BatchApplyFeatureDetailsPtr BatchApplyFeatureDetailsNew (void)
{
  BatchApplyFeatureDetailsPtr    bafdp;

  bafdp = (BatchApplyFeatureDetailsPtr) MemNew (sizeof (BatchApplyFeatureDetailsData));
  if (bafdp != NULL)
  {
    bafdp->defline = NULL;
    bafdp->geneName = NULL;
	bafdp->geneDesc = NULL;
    bafdp->protName = NULL;
    bafdp->protDesc = NULL;
    bafdp->rnaName = NULL;
    bafdp->featcomment = NULL;
    bafdp->featdef_name = NULL;
    bafdp->featdef_choice = FEATDEF_GENE;
    bafdp->reading_frame = 4;
    bafdp->rnaSubType = 4;
  }
  return bafdp;
}

static ENUM_ALIST(rnax_subtype_alist)
  {" ",           99},
  {"unknown",      0},
  {"preRna",       1},
  {"mRNA",         2},
  {"tRNA",         3},
  {"rRNA",         4},
  {"snRNA",        5},
  {"scRNA",        6},
  {"snoRNA",       7},
  {"misc_RNA",   255},
END_ENUM_ALIST

static void BatchApplyFeatureDetailsToDialog (DialoG d, Pointer data)
{
  BatchApplyFeatureDetailsDlgPtr dlg;
  BatchApplyFeatureDetailsPtr    bafdp;
  ValNode                        vn;
  
  dlg = (BatchApplyFeatureDetailsDlgPtr) GetObjectExtra (d);
  bafdp = (BatchApplyFeatureDetailsPtr) data;
  
  if (dlg == NULL)
  {
    return;
  }
  
  vn.next = NULL;

  if (bafdp == NULL)
  {
    SafeSetTitle (dlg->defline, "");
    SafeSetTitle (dlg->geneName, "");
	SafeSetTitle (dlg->geneDesc, "");
    SafeSetTitle (dlg->protName, "");
    SafeSetTitle (dlg->protDesc, "");
    SafeSetTitle (dlg->rnaName, "");
    SafeSetTitle (dlg->featcomment, "");
    SafeSetValue (dlg->reading_frame, 4);
    if (dlg->featdef_choice_dlg != NULL)
    {
      vn.choice = FEATDEF_GENE;
      vn.data.ptrvalue = NULL;
      PointerToDialog (dlg->featdef_choice_dlg, &vn);
    }
    SafeSetValue (dlg->rnaSubType, 4);
  }
  else
  {
    if (StringHasNoText (bafdp->defline))
    {
      SafeSetTitle (dlg->defline, "");
    }
    else
    {
      SafeSetTitle (dlg->defline, bafdp->defline);
    }

    if (StringHasNoText (bafdp->geneName))
    {
      SafeSetTitle (dlg->geneName, "");
    }
    else
    {
      SafeSetTitle (dlg->geneName, bafdp->geneName);
    }
    
    if (StringHasNoText (bafdp->protName))
    {
      SafeSetTitle (dlg->protName, "");
    }
    else
    {
      SafeSetTitle (dlg->protName, bafdp->protName);
    }
    
    if (StringHasNoText (bafdp->protDesc))
    {
      SafeSetTitle (dlg->protDesc, "");
    }
    else
    {
      SafeSetTitle (dlg->protDesc, bafdp->protDesc);
    }

    if (StringHasNoText (bafdp->rnaName))
    {
      SafeSetTitle (dlg->rnaName, "");
    }
    else
    {
      SafeSetTitle (dlg->rnaName, bafdp->rnaName);
    }

    if (StringHasNoText (bafdp->featcomment))
    {
      SafeSetTitle (dlg->featcomment, "");
    }
    else
    {
      SafeSetTitle (dlg->featcomment, bafdp->featcomment);
    }

    SafeSetValue (dlg->reading_frame, bafdp->reading_frame);
    if (dlg->featdef_choice_dlg != NULL)
    {
      vn.choice = bafdp->featdef_choice;
      vn.data.ptrvalue = NULL;
      PointerToDialog (dlg->featdef_choice_dlg, &vn);
    }
    SafeSetValue (dlg->rnaSubType, bafdp->rnaSubType);
  }  
}

static Pointer BatchApplyFeatureDetailsDialogToData (DialoG d)
{
  BatchApplyFeatureDetailsDlgPtr dlg;
  BatchApplyFeatureDetailsPtr    bafdp;
  ValNodePtr                     vnp;
  UIEnum                         val;
  
  dlg = (BatchApplyFeatureDetailsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  bafdp = BatchApplyFeatureDetailsNew ();
  if (bafdp == NULL)
  {
    return NULL;
  }
  
  if (dlg->defline == NULL || TextHasNoText (dlg->defline))
  {
    bafdp->defline = NULL;
  }
  else
  {
    bafdp->defline = SaveStringFromTextAndStripNewlines (dlg->defline);
  }
  
  if (dlg->geneName == NULL || TextHasNoText (dlg->geneName))
  {
    bafdp->geneName = NULL;
  }
  else
  {
    bafdp->geneName = SaveStringFromText (dlg->geneName);
  }
  
  if (dlg->geneDesc == NULL || TextHasNoText (dlg->geneDesc))
  {
    bafdp->geneDesc = NULL;
  }
  else
  {
    bafdp->geneDesc = SaveStringFromText (dlg->geneDesc);
  }
  
 if (dlg->protName == NULL || TextHasNoText (dlg->protName))
  {
    bafdp->protName = NULL;
  }
  else
  {
    bafdp->protName = SaveStringFromText (dlg->protName);
  }
  
  if (dlg->protDesc == NULL || TextHasNoText (dlg->protDesc))
  {
    bafdp->protDesc = NULL;
  }
  else
  {
    bafdp->protDesc = SaveStringFromText (dlg->protDesc);
  }
  
  if (dlg->rnaName == NULL || TextHasNoText (dlg->rnaName))
  {
    bafdp->rnaName = NULL;
  }
  else
  {
    bafdp->rnaName = SaveStringFromText (dlg->rnaName);
  }
  
  if (dlg->featcomment == NULL || TextHasNoText (dlg->featcomment))
  {
    bafdp->featcomment = NULL;
  }
  else
  {
    bafdp->featcomment = SaveStringFromTextAndStripNewlines (dlg->featcomment);
  }
  
  if (dlg->reading_frame == NULL)
  {
    bafdp->reading_frame = 4;
  }
  else
  {
    bafdp->reading_frame = GetValue (dlg->reading_frame);
  }
  
  if (dlg->featdef_choice_dlg == NULL)
  {
    bafdp->featdef_choice = FEATDEF_GENE;
  }
  else
  {
    vnp = DialogToPointer (dlg->featdef_choice_dlg);
    if (vnp != NULL)
    {
      bafdp->featdef_choice = vnp->choice;
      bafdp->featdef_name = StringSave (vnp->data.ptrvalue);
      vnp = ValNodeFreeData (vnp);
    }
    else
    {
      bafdp->featdef_choice = FEATDEF_GENE;
      bafdp->featdef_name = StringSave ("Gene");
    }
  }

  if (dlg->rnaSubType == NULL)
  {
    bafdp->rnaSubType = 4;
  }
  else
  {
    if (GetEnumPopup (dlg->rnaSubType, rnax_subtype_alist, &val))
    {
      bafdp->rnaSubType = val;
    }
    else
    {
      bafdp->rnaSubType = 4;
    }
  }
  
  return bafdp;  
}

static void BatchApplyFeatureDetailsMessage (DialoG d, Int2 mssg)

{
  BatchApplyFeatureDetailsDlgPtr  dlg;

  dlg = (BatchApplyFeatureDetailsDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        PointerToDialog (d, NULL);
        break;
      case VIB_MSG_ENTER :
        if (dlg->feattype == ADD_TITLE) {
          Select (dlg->defline);
        } else if (dlg->feattype == ADD_RRNA) {
          Select (dlg->rnaName);
        } else {
          Select (dlg->geneName);
        }
        break;
      default :
        break;
    }
  }
}

/* This section of code is used for managing lists of features.
 * Sometimes the features will be displayed alphabetically, sometimes
 * they will be displayed alphabetically with a list of the most used features
 * also appearing at the top of the list.
 */

/* This is used to compare feature names with the special alphabetical order */
static int CompareFeatureNames (CharPtr cp1, CharPtr cp2)
{
  /* NULL name goes at the end */
  if (cp1 == NULL && cp2 == NULL) return 0;
  if (cp1 == NULL) return 1;
  if (cp2 == NULL) return -1;

  /* starts with a space goes at the beginning */
  if (cp1 [0] == ' ' && cp2 [0] == ' ') return 0;
  if (cp1 [0] == ' ') return -1;
  if (cp2 [0] == ' ') return 1;

  /* Is "All" or [ALL FEATURES] goes at the beginning */
  if ((StringCmp (cp1, "All") == 0
    || StringCmp (cp1, "[ALL FEATURES]") == 0)
    && (StringCmp (cp2, "All") == 0
    || StringCmp (cp2, "[ALL FEATURES]") == 0))
  {
    return 0;
  }
  if (StringCmp (cp1, "All") == 0
    || StringCmp (cp1, "[ALL FEATURES]") == 0)
  {
    return -1;
  }
  if (StringCmp (cp2, "All") == 0
    || StringCmp (cp2, "[ALL FEATURES]") == 0)
  {
    return 1;
  }

  /* starts with a number -> goes at the end */
  if (cp1 [0] >= '0' && cp1 [0] <= '9'
   && cp2 [0] >= '0' && cp2 [0] <= '9')
  {
    return StringICmp (cp1, cp2);
  }
  if (cp1 [0] >= '0' && cp1 [0] <= '9')
  {
    return 1;
  }
  if (cp2 [0] >= '0' && cp2 [0] <= '9')
  {
    return -1;
  }

  /* starts with a tilde or dash - sort with other tildes, put before numbers after alphas */
  if (cp1 [0] == '~' && cp2 [0] == '~') 
  {
    return StringICmp (cp1 + 1, cp2 + 1);
  }
  if (cp1 [0] == '~') return 1;
  if (cp2 [0] == '~') return -1;

  if (cp1 [0] == '-' && cp2 [0] == '-') 
  {
    return StringICmp (cp1 + 1, cp2 + 1);
  }
  if (cp1 [0] == '-') return 1;
  if (cp2 [0] == '-') return -1;

  return StringICmp (cp1, cp2);
}

extern int LIBCALLBACK CompareFeatureValNodeStrings (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);

  if (vnp1 == NULL || vnp2 == NULL) return 0;

  return CompareFeatureNames (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
}

extern int LIBCALLBACK CompareImpFeatEnumFieldAssoc (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr        vnp1, vnp2;
  EnumFieldAssocPtr ap1, ap2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  ap1 = (EnumFieldAssocPtr) vnp1->data.ptrvalue;
  ap2 = (EnumFieldAssocPtr) vnp2->data.ptrvalue;
  if (ap1 == NULL || ap2 == NULL) return 0;

  return CompareFeatureNames (ap1->name, ap2->name);
}

extern void SortEnumFieldAssocPtrArray (EnumFieldAssocPtr alist, CompareFunc compar)
{
  ValNodePtr        head, vnp;
  EnumFieldAssocPtr ap;
  Int4              index;

  /* first, create ValNode list so we can sort the data */
  head = NULL;
  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    vnp = ValNodeNew (head);
    if (vnp == NULL) return;
    vnp->data.ptrvalue = MemNew (sizeof (EnumFieldAssoc));
    if (vnp->data.ptrvalue == NULL) return;
    MemCpy (vnp->data.ptrvalue, ap, sizeof (EnumFieldAssoc));
    if (head == NULL) head = vnp;
  }

  /* Now sort the ValNode list */
  head = SortValNode (head, compar);

  /* Now repopulate the EnumFieldAssoc list */
  index = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next)
  {
    MemCpy (alist + index++, vnp->data.ptrvalue, sizeof (EnumFieldAssoc));
  }

  /* And free the ValNode list */
  ValNodeFreeData (head);
}

static void AddTextToComment (ButtoN b, CharPtr text)
{
  BatchApplyFeatureDetailsDlgPtr dlg;
  CharPtr                        orig_comment;
  CharPtr                        new_comment;
  
  dlg = (BatchApplyFeatureDetailsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL || StringHasNoText (text))
  {
    return;
  }
  
  orig_comment = SaveStringFromText (dlg->featcomment);
  if (StringHasNoText (orig_comment))
  {
    SetTitle (dlg->featcomment, text);
  }
  else
  {
    new_comment = (CharPtr) MemNew ((StringLen (orig_comment) + StringLen (text) + 3) * sizeof (Char));
    if (new_comment != NULL)
    {
      StringCpy (new_comment, orig_comment);
      StringCat (new_comment, "; ");
      StringCat (new_comment, text);
      SetTitle (dlg->featcomment, new_comment);
      new_comment = MemFree (new_comment);
    }
  } 
  orig_comment = MemFree (orig_comment);
}

static void Add18SITS28SToComment (ButtoN b)
{
  AddTextToComment (b, "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA");
}

static void Add16SIGS23SToComment (ButtoN b)
{
  AddTextToComment (b, "contains 16S ribosomal RNA, 16S-23S ribosomal RNA intergenic spacer, and 23S ribosomal RNA");
}

extern DialoG 
BatchApplyFeatureDetailsDialog (GrouP parent, Int4 feattype)
{
  BatchApplyFeatureDetailsDlgPtr dlg;
  GrouP                          p, r = NULL, text_group = NULL, comment_btns_grp = NULL;
  ButtoN                         comment_btn;
  Nlm_EnumFieldAssocPtr          ap;
  Int4                           j;
  Boolean                        is_indexer = FALSE;
  
  dlg = (BatchApplyFeatureDetailsDlgPtr) MemNew (sizeof (BatchApplyFeatureDetailsDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  if (GetAppProperty ("InternalNcbiSequin") != NULL)
  {
    is_indexer = TRUE;
  }

  p = HiddenGroup (parent, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = BatchApplyFeatureDetailsToDialog;
  dlg->fromdialog = BatchApplyFeatureDetailsDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  dlg->feattype = feattype;
  
  /* codon start controls */
  if (dlg->feattype == ADD_CDS)
  {
    r = HiddenGroup (p, 2, 0, NULL);
    StaticPrompt (r, "Codon Start", 0, dialogTextHeight, programFont, 'l');
    dlg->reading_frame = PopupList (r, TRUE, NULL);
    PopupItem (dlg->reading_frame, "1");
    PopupItem (dlg->reading_frame, "2");
    PopupItem (dlg->reading_frame, "3");
    PopupItem (dlg->reading_frame, "Best");
    SetValue (dlg->reading_frame, 4);
  }
  else if (dlg->feattype == ADD_RRNA) 
  {
    r = HiddenGroup (p, 2, 0, NULL);
    StaticPrompt (r, "RNA subtype", 0, dialogTextHeight, programFont, 'l');
    dlg->rnaSubType = PopupList (r, TRUE, NULL);
    InitEnumPopup (dlg->rnaSubType, rnax_subtype_alist, NULL);
    SetEnumPopup (dlg->rnaSubType, rnax_subtype_alist, (UIEnum) 4);
  }

  text_group = HiddenGroup (p, 0, 2, NULL);
  if (dlg->feattype == ADD_TITLE) {
    StaticPrompt (text_group, "Title", 0, 0, programFont, 'c');
    dlg->defline = ScrollText (text_group, 20, 4, programFont, TRUE, NULL);
  } else {
    text_group = HiddenGroup (p, 2, 0, NULL);
    if (dlg->feattype == ADD_CDS) {
      StaticPrompt (text_group, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      dlg->geneName = DialogText (text_group, "", 20, NULL);
	  if (is_indexer)
	  {
		StaticPrompt (text_group, "Gene Description", 0, dialogTextHeight, programFont, 'l');
        dlg->geneDesc = DialogText (text_group, "", 20, NULL);
	  }

      StaticPrompt (text_group, "Protein Name", 0, dialogTextHeight, programFont, 'l');
      dlg->protName = DialogText (text_group, "", 20, NULL);
      StaticPrompt (text_group, "Protein Description", 0, dialogTextHeight, programFont, 'l');
      dlg->protDesc = DialogText (text_group, "", 20, NULL);
    } else if (dlg->feattype == ADD_RRNA) {
      StaticPrompt (text_group, "RNA Name", 0, dialogTextHeight, programFont, 'l');
      dlg->rnaName = DialogText (text_group, "", 20, NULL);
      StaticPrompt (text_group, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      dlg->geneName = DialogText (text_group, "", 20, NULL);
	  if (is_indexer)
	  {
		StaticPrompt (text_group, "Gene Description", 0, dialogTextHeight, programFont, 'l');
        dlg->geneDesc = DialogText (text_group, "", 20, NULL);
	  }
    } else if (dlg->feattype == ADD_IMP) {
      StaticPrompt (text_group, "Type", 0, 6 * Nlm_stdLineHeight, programFont, 'l');
      ap = import_featdef_alist (FALSE, FALSE, FALSE);
      SortEnumFieldAssocPtrArray (ap, CompareImpFeatEnumFieldAssoc);
      /* replace first item with Gene */
      ap [0].name = MemFree (ap [0].name);
      ap [0].name = StringSave ("Gene");
      ap [0].value = FEATDEF_GENE;
      
      dlg->featdef_choice_dlg = EnumAssocSelectionDialog (text_group, ap, 
                                                  "feat_detail",
                                                   FALSE, NULL, NULL);
      /* clean up enumassoc list - not needed any more */
      for (j = 0; ap [j].name != NULL; j++) {
        MemFree (ap [j].name);
      }
      MemFree (ap);

      StaticPrompt (text_group, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      dlg->geneName = DialogText (text_group, "", 20, NULL);
	  if (is_indexer)
	  {
		StaticPrompt (text_group, "Gene Description", 0, dialogTextHeight, programFont, 'l');
        dlg->geneDesc = DialogText (text_group, "", 20, NULL);
	  }
    }
    StaticPrompt (text_group, "Comment", 0, 4 * Nlm_stdLineHeight, programFont, 'l');
    dlg->featcomment = ScrollText (text_group, 20, 4, programFont, TRUE, NULL);
    
  }
  
  if (dlg->feattype == ADD_RRNA && is_indexer)
  {
    comment_btns_grp = HiddenGroup (p, 2, 0, NULL);
    comment_btn = PushButton (comment_btns_grp, "Add '18S-ITS-5.8S-ITS-28S' to comment", Add18SITS28SToComment);
    SetObjectExtra (comment_btn, dlg, NULL);
    comment_btn = PushButton (comment_btns_grp, "Add '16S-IGS-23S' to comment", Add16SIGS23SToComment);
    SetObjectExtra (comment_btn, dlg, NULL);
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) text_group, (HANDLE) r, (HANDLE) comment_btns_grp, NULL);
  
  return (DialoG) p;
}

typedef struct alreadyhas {
  Boolean        rsult;
  Uint1          featchoice;
  Uint1          descchoice;
  Uint1          rnatype;
} AlreadyHas, PNTR AlreadyHasPtr;

static Boolean SeeIfAlreadyHasGatherFunc (GatherContextPtr gcp)

{
  AlreadyHasPtr  ahp;
  RnaRefPtr      rrp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;

  ahp = (AlreadyHasPtr) gcp->userdata;
  if (ahp == NULL ) return TRUE;

  if (gcp->thistype == OBJ_SEQFEAT && ahp->featchoice != 0) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == ahp->featchoice && sfp->data.value.ptrvalue != NULL) {
      if (sfp->data.choice == SEQFEAT_RNA) {
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp->type != ahp->rnatype) return TRUE;
      }
      ahp->rsult = TRUE;
      return FALSE;
    }
  } else if (gcp->thistype == OBJ_SEQDESC && ahp->descchoice != 0) {
    sdp = (ValNodePtr) gcp->thisitem;
    if (sdp != NULL && sdp->choice == ahp->descchoice && sdp->data.ptrvalue != NULL) {
      ahp->rsult = TRUE;
      return FALSE;
    }
  }
  return TRUE;
}

static Boolean AlreadyHasFeatOrDesc (SeqEntryPtr sep, Uint1 featchoice, Uint1 descchoice, Uint1 rnatype)

{
  AlreadyHas   ah;
  BioseqPtr    bsp;
  GatherScope  gs;
  SeqEntryPtr  nsep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  ah.rsult = FALSE;
  ah.featchoice = featchoice;
  ah.descchoice = descchoice;
  ah.rnatype = rnatype;
  if (sep == NULL) return FALSE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  if (descchoice != 0) {
    nsep = FindNucSeqEntry (sep);
    if (nsep != NULL && IS_Bioseq (nsep)) {
      bsp = (BioseqPtr) nsep->data.ptrvalue;
      if (bsp != NULL) {
        slp = ValNodeNew (NULL);
        slp->choice = SEQLOC_WHOLE;
        sip = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
        slp->data.ptrvalue = sip;
        gs.target = slp;
      }
    }
  }
  GatherSeqEntry (sep, (Pointer) (&ah), SeeIfAlreadyHasGatherFunc, &gs);
  gs.target = SeqLocFree (gs.target);
  return ah.rsult;
}

static void AddGeneXrefToFeat (SeqFeatPtr sfp, CharPtr str)
{
  SeqFeatXrefPtr    xref;
  GeneRefPtr        grp;
  
  if (sfp == NULL || StringHasNoText (str)) return;
  
  /* add gene xref to feature */
  xref = SeqFeatXrefNew ();
  if (xref != NULL)
  {
    grp = CreateNewGeneRef (str, NULL, NULL, FALSE);
    if (grp != NULL) 
    {
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = grp;
      xref->next = sfp->xref;
      sfp->xref = xref;
    }
  }
}

static SeqFeatPtr ApplyGene 
(CharPtr gene_name,
 BatchApplyFeatureDetailsPtr feature_details_data,
 SeqEntryPtr gene_sep, 
 SeqFeatPtr sfp,
 SeqLocPtr  slp)
{
  GeneRefPtr        grp;
  SeqFeatPtr        gene_sfp;
  SeqFeatXrefPtr    xref;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp = NULL;
  SeqFeatPtr        other_feat;
  SeqFeatPtr        overlap_gene;
  Boolean           added_xrefs = FALSE;
  SeqFeatPtr        misc_feat = NULL;
  SeqLocPtr         overlap_loc;
  Boolean           partial5, partial3;

  if (feature_details_data == NULL || gene_sep == NULL || slp == NULL 
	  || (StringHasNoText (gene_name) && StringHasNoText (feature_details_data->geneDesc))) return NULL;

  CheckSeqLocForPartial (slp, &partial5, &partial3);

  /* we need a location to use when we're checking for feature-stealing genes */
  if (sfp != NULL)
  {
    overlap_loc = sfp->location;
  }
  else
  {
    misc_feat = CreateNewFeature (gene_sep, NULL, SEQFEAT_COMMENT, NULL);
    if (NULL == misc_feat)
    return NULL;
    misc_feat->location = SeqLocCopy (slp);

    misc_feat->partial = (partial5 || partial3);
    overlap_loc = misc_feat->location;
  }
  
  /* first, add gene xrefs to all features on bioseq that are contained in the location */
  /* maintain list of features that had xrefs before, should not remove them later */
  if (IS_Bioseq (gene_sep))
  {
    bsp = (BioseqPtr) gene_sep->data.ptrvalue;
  }
  else if (sfp != NULL)
  {
    bsp = BioseqFindFromSeqLoc (sfp->location);
  }
  if (bsp != NULL)
  {
    other_feat = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (other_feat != NULL)
    {
      if (other_feat != sfp && other_feat->data.choice != SEQFEAT_GENE)
      {
        for (xref = other_feat->xref;
             xref != NULL && xref->data.choice != SEQFEAT_GENE;
             xref = xref->next)
        {}
        if (xref == NULL
            && SeqLocCompare (other_feat->location, overlap_loc) == SLC_A_EQ_B)
        {
          overlap_gene = SeqMgrGetOverlappingGene (other_feat->location, NULL);
          if (overlap_gene != NULL)
          {
            AddGeneXrefToFeat (other_feat, fcontext.label);
            added_xrefs = TRUE;
          }
        }
      }
      other_feat = SeqMgrGetNextFeature (bsp, other_feat, 0, 0, &fcontext);
    }   
  }
  
  if (misc_feat != NULL)
  {
    misc_feat->idx.deleteme = TRUE;
    DeleteMarkedObjects (0, OBJ_SEQENTRY, gene_sep);
  }
  
  grp = CreateNewGeneRef (gene_name, NULL, feature_details_data->geneDesc, FALSE);
  if (NULL == grp)
    return NULL;

  gene_sfp = CreateNewFeature (gene_sep, NULL, SEQFEAT_GENE, NULL);
  if (NULL == gene_sfp)
    return NULL;

  gene_sfp->data.value.ptrvalue = (Pointer) grp;
  gene_sfp->location = slp;
  gene_sfp->partial = partial5 | partial3;

  if (added_xrefs && sfp != NULL)
  {
    /* add gene xref to feature */
    AddGeneXrefToFeat (sfp, gene_name);
  }
  
  return gene_sfp;
}

static void 
SetFrameForCodingRegion 
(SeqFeatPtr                  sfp,
 BioseqPtr                   bsp,
 BatchApplyFeatureDetailsPtr feature_details_data,
 Int4Ptr                     errcount,
 ValNodePtr PNTR             ambigList)
{
  CdRegionPtr        crp;
  ByteStorePtr       bs;
  Uint1              frame;
  Int2               i;
  Int4               len;
  Int4               lens [4];
  Int4               max;
  Char               str [128];
  SeqIdPtr           sip;

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION 
      || sfp->data.value.ptrvalue == NULL
      || bsp == NULL
      || feature_details_data == NULL
      || errcount == NULL || ambigList == NULL)
  {
    return;
  }

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;

  if (feature_details_data->reading_frame < 1 
      || feature_details_data->reading_frame > 3)
  {
    max = 0;
    frame = 0;
    for (i = 1; i <= 3; i++) {
      crp->frame = (Uint1) i;
      bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
      len = BSLen (bs);
      BSFree (bs);
      lens [i] = len;
      if (len > max) {
        max = len;
        frame = (Uint1) i;
      }
    }
    str [0] = '\0';
    sip = SeqIdFindBest (bsp->id, 0);
    SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str) - 1);
    for (i = 1; i <= 3; i++) {
      if (lens [i] == max && i != frame) {
        (*errcount)++;
        ValNodeCopyStr (ambigList, 0, str);
      }
    }
    crp->frame = frame;
  }
  else
  {
    crp->frame = feature_details_data->reading_frame;
  }
}

static void 
AddProductForCDS 
(SeqFeatPtr                  sfp,
 BatchApplyFeatureDetailsPtr feature_details_data,
 SeqEntryPtr                 sep,
 SeqEntryPtr                 nsep,
 Uint2                       entityID)
{
  ByteStorePtr       bs;
  Char               ch;
  ValNodePtr         descr;
  Int2               i;
  MolInfoPtr         mip;
  SeqEntryPtr        old;
  CharPtr            prot;
  ProtRefPtr         prp;
  SeqEntryPtr        psep;
  CharPtr            ptr;
  ValNodePtr         vnp;
  SeqEntryPtr        parent_sep;
  SeqFeatPtr         prot_sfp;  
  BioseqPtr          bsp;
  Boolean            partial5, partial3;
  
  if (sfp == NULL
      || sfp->data.choice != SEQFEAT_CDREGION
      || feature_details_data == NULL
      || sep == NULL || nsep == NULL)
  {
    return;
  }
  
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  
  /* determine the parent of this sequence (for use when segmented) */
  parent_sep = NULL;
  if (IS_Bioseq (sep))
  {
    parent_sep = GetBestTopParentForData (entityID, sep->data.ptrvalue);
  }
  if (parent_sep == NULL)
  {
    parent_sep = sep;
  }
  
  /* Create corresponding protein sequence data for the CDS */

  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (NULL == bs)
    return;

  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (NULL == prot)
    return;

  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = (Int2) StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }

  /* Create the product protein Bioseq */
  
  bsp = BioseqNew ();
  if (NULL == bsp)
    return;
  
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (NULL);
  bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);
  
  /* Create a new SeqEntry for the Prot Bioseq */
  
  psep = SeqEntryNew ();
  if (NULL == psep)
    return;
  
  psep->choice = 1;
  psep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
  
  /* Add a descriptor to the protein Bioseq */
  
  mip = MolInfoNew ();
  if (NULL == mip)
    return;
  
  mip->biomol = 8;
  mip->tech = 8;
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
  if (NULL == vnp)
    return;
  
  vnp->data.ptrvalue = (Pointer) mip;
  
  /**/
  
  descr = ExtractBioSourceAndPubs (parent_sep);

  AddSeqEntryToSeqEntry (parent_sep, psep, TRUE);
  nsep = FindNucSeqEntry (parent_sep);
  ReplaceBioSourceAndPubs (parent_sep, descr);
  SetSeqFeatProduct (sfp, bsp);
  
  /* create a full-length protein feature for the new protein sequence */
  if (! StringHasNoText (feature_details_data->protName) 
      && ! StringHasNoText (feature_details_data->protDesc))
  {
    prp = CreateNewProtRef (feature_details_data->protName, 
                            feature_details_data->protDesc, 
                            NULL, NULL);
  }
  else if (!StringHasNoText (feature_details_data->protName))
  {
    prp = CreateNewProtRef (feature_details_data->protName, NULL, NULL, NULL);
  }
  else if (!StringHasNoText (feature_details_data->protDesc))
  {
    prp = CreateNewProtRef (NULL, feature_details_data->protDesc, NULL, NULL);
  }
  else
  { 
    prp = ProtRefNew ();
  }
  
  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    }
  }
}


static SeqFeatPtr 
ApplyOneCodingRegion 
(BatchApplyFeatureDetailsPtr feature_details_data,
 BioseqPtr                   bsp,
 SeqEntryPtr                 sep,
 SeqEntryPtr                 nsep,
 Uint2                       entityID,
 Boolean                     suppressDups)
{
  CdRegionPtr        crp;
  Int2               genCode;
  SeqFeatPtr         sfp;
  SeqEntryPtr        parent_sep;

  if (feature_details_data == NULL || sep == NULL || nsep == NULL)
  {
    return NULL;
  }

  /* If necessary then check for duplication before adding */
  if (suppressDups &&
      entityID > 0 &&
      AlreadyHasFeatOrDesc (sep, SEQFEAT_CDREGION, 0, 0))
    return NULL;
  
  /* determine the parent of this sequence (for use when segmented) */
  parent_sep = NULL;
  if (IS_Bioseq (sep))
  {
    parent_sep = GetBestTopParentForData (entityID, sep->data.ptrvalue);
  }
  if (parent_sep == NULL)
  {
    parent_sep = sep;
  }
  
  /*Create a new CDS feature */

  genCode = SeqEntryToGeneticCode (parent_sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (NULL == crp)
    return NULL;
  
  sfp = CreateNewFeature (nsep, NULL, SEQFEAT_CDREGION, NULL); 	

  if (NULL == sfp)
    return NULL;
  
  sfp->data.value.ptrvalue = (Pointer) crp;

  return sfp;
}

static SeqFeatPtr ApplyOneRNA 
(BatchApplyFeatureDetailsPtr feature_details_data,
 SeqEntryPtr                 sep, 
 SeqEntryPtr                 nsep, 
 Uint2                       entityID,
 Boolean                     suppressDups)
{
  RnaRefPtr  rrp;
  SeqFeatPtr sfp;
  
  if (feature_details_data == NULL || sep == NULL || nsep == NULL)
  {
    return NULL;
  }
  
  if (suppressDups && entityID > 0 &&
	     AlreadyHasFeatOrDesc (sep, SEQFEAT_RNA, 0, feature_details_data->rnaSubType))
  {
    return NULL;
  }
  
  rrp = RnaRefNew ();
  if (rrp == NULL)
  {
    return NULL;
  }

  rrp->type = feature_details_data->rnaSubType;
  if (! StringHasNoText (feature_details_data->rnaName)) {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave (feature_details_data->rnaName);
  }
  sfp = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, NULL);
  if (sfp != NULL) {
    sfp->data.value.ptrvalue = (Pointer) rrp;
  }
  return sfp;
}

static GBQualPtr GBQualListCopy (GBQualPtr gbqual_list)
{
  GBQualPtr gbqual_orig, gbqual_newlist = NULL, gbqual_new = NULL;
  
  gbqual_orig = gbqual_list;
  while (gbqual_orig != NULL)
  {
    if (gbqual_newlist == NULL)
    {
      gbqual_newlist = GBQualNew ();
      gbqual_new = gbqual_newlist;
    }
    else
    {
      gbqual_new->next = GBQualNew ();
      gbqual_new = gbqual_new->next;
    }
    
    if (gbqual_new == NULL)
    {
      gbqual_orig = NULL;
    }
    else
    {
      gbqual_new->qual = StringSave (gbqual_orig->qual);
      gbqual_new->val = StringSave (gbqual_orig->val);
      gbqual_orig = gbqual_orig->next;
    }
  }
    
  return gbqual_newlist;
}

static SeqFeatPtr ApplyOtherFeature 
(BatchApplyFeatureDetailsPtr feature_details_data,
 SeqEntryPtr                 nsep,
 GBQualPtr                   gbqual_list)
{
  ImpFeatPtr ifp;
  SeqFeatPtr sfp = NULL;
  
  if (feature_details_data == NULL 
      || nsep == NULL 
      || feature_details_data->featdef_choice == FEATDEF_GENE)
  {
    return NULL;
  }
    
  ifp = ImpFeatNew ();
  if (ifp == NULL) 
  {
    return NULL;
  }
  
  ifp->key = StringSave (feature_details_data->featdef_name);
  sfp = CreateNewFeature (nsep, NULL, SEQFEAT_IMP, NULL);
  if (sfp != NULL) {
    sfp->data.value.ptrvalue = (Pointer) ifp;
    sfp->qual = GBQualListCopy (gbqual_list);
  }
  return sfp;  
}

static GBQualPtr GetImpFeatureGBQuals (BatchApplyFeatureDetailsPtr feature_details_data)
{
  WindoW                w;
  DialoG                dlg;
  GrouP                 h, c;
  ButtoN                b;
  ModalAcceptCancelData acd;
  GBQualPtr             gbqual_list = NULL;
  
  if (feature_details_data == NULL 
      || feature_details_data->featdef_choice == FEATDEF_GENE)
  {
    return NULL;
  }
  
  w = MovableModalWindow(-20, -13, -10, -10, "Qualifiers", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  dlg = CreateImportFields (h, feature_details_data->featdef_name, NULL, FALSE);
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton(c, "OK", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton(c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  
  if (acd.accepted)
  {
    gbqual_list = DialogToPointer (dlg);
  }
  
  Remove (w);
  
  return gbqual_list;
}

typedef struct applyalignmentfeature
{
  FORM_MESSAGE_BLOCK
  DialoG      location_dlg;
  DialoG      feature_details;
  SeqEntryPtr sep;
  SeqAlignPtr salp;
  Int4        feattype;
  Uint2       entityID;
} ApplyAlignmentFeatureDlgData, PNTR ApplyAlignmentFeatureDlgPtr;

static Boolean 
OkToContinueWithFeatureApply 
(SeqLocPtr   seqloc_list,
 SeqAlignPtr salp)
{
  Int4       num_missing = 0, seq_num, start_num;
  SeqLocPtr  slp;
  SeqIdPtr   sip;
  Int4Ptr    found_list;
  Boolean    found_this;
  Int4       msg_len = 0;
  CharPtr    msg_txt = NULL;
  CharPtr    msg_start = "The following sequence(s) are all gaps at this position: ";
  CharPtr    msg_end = ". Do you want to continue?";
  ValNodePtr id_list = NULL, id_vnp;
  Char       id_txt [255];
  Boolean    rval = TRUE;
  BioseqPtr  bsp;
  
  if (salp == NULL || seqloc_list == NULL)
  {
    return FALSE;
  }
  
  found_list = (Int4Ptr) MemNew (salp->dim * sizeof (Int4));
  if (found_list == NULL)
  {
    return FALSE;
  }
  
  for (seq_num = 1; seq_num <= salp->dim; seq_num++)
  {
    found_list [seq_num - 1] = 0;
  }
  
  for (slp = seqloc_list, start_num = 1; slp != NULL; slp = slp->next, start_num++)
  {
    if (slp->choice == SEQLOC_NULL)
    {
      continue;
    }
    found_this = FALSE;
    if (start_num <= salp->dim)
    {
      sip = AlnMgr2GetNthSeqIdPtr (salp, start_num);
      if (SeqIdComp (sip, SeqLocId (slp)) == SIC_YES)
      {
        found_list [start_num - 1] = 1;
        found_this = TRUE;
      }
    }
    
    for (seq_num = 1; seq_num <= salp->dim && ! found_this; seq_num++)
    {
      sip = AlnMgr2GetNthSeqIdPtr (salp, seq_num);
      if (SeqIdComp (sip, SeqLocId (slp)) == SIC_YES)
      {
        found_list [seq_num - 1] = 1;
        found_this = TRUE;
      }
    }
  }
  
  for (seq_num = 1; seq_num <= salp->dim; seq_num++)
  {
    if (found_list [seq_num - 1] == 0)
    {
      sip = AlnMgr2GetNthSeqIdPtr (salp, seq_num);
      bsp = BioseqFind (sip);
      if (bsp != NULL && bsp->id != NULL)
      {
        sip = SeqIdFindBest (bsp->id, 0);
      }
      SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt));
      msg_len += StringLen (id_txt) + 3;
      ValNodeAddPointer (&id_list, 0, StringSave (id_txt));
    }
  }
  
  if (msg_len > 0)
  {
    msg_len += StringLen (msg_start) + StringLen (msg_end);
    msg_txt = (CharPtr) MemNew (msg_len * sizeof (Char));
    if (msg_txt != NULL)
    {
      StringCpy (msg_txt, msg_start);
      for (id_vnp = id_list; id_vnp != NULL; id_vnp = id_vnp->next)
      {
        StringCat (msg_txt, id_vnp->data.ptrvalue);
        if (id_vnp->next != NULL)
        {
          StringCat (msg_txt, ", ");
        }
      }
      StringCat (msg_txt, msg_end);
      if (ANS_NO == Message (MSG_YN, msg_txt))
      {
        rval = FALSE;
      }
      msg_txt = MemFree (msg_txt);
    }
  }
  id_list = ValNodeFreeData (id_list);
  found_list = MemFree (found_list);
  
  return rval;
}

static void DoApplyFeatureToAlignment (ButtoN b)
{
  ApplyAlignmentFeatureDlgPtr aafdp;
  BatchApplyFeatureDetailsPtr feature_details_data;
  SeqLocPtr                   seqloc_list, slp_next, slp_this, gene_slp;
  Int4                        errcount = 0;
  ValNodePtr                  ambigList = NULL;
  Boolean                     suppressDups = FALSE;
  Boolean                     partial5, partial3;
  SeqEntryPtr                 sep, nsep, gene_sep;
  BioseqPtr                   bsp;
  SeqFeatPtr                  sfp, gene_sfp;
  GBQualPtr                   gbqual_list = NULL;
  Boolean                     found_null = FALSE;
  
  aafdp = (ApplyAlignmentFeatureDlgPtr) GetObjectExtra (b);
  if (aafdp == NULL)
  {
    return;
  }
  
  feature_details_data = (BatchApplyFeatureDetailsPtr) DialogToPointer (aafdp->feature_details);
  seqloc_list = (SeqLocPtr) DialogToPointer (aafdp->location_dlg);

  if (! OkToContinueWithFeatureApply (seqloc_list, aafdp->salp))
  {
    /* Free data and loclist */
    feature_details_data = BatchApplyFeatureDetailsFree (feature_details_data);
    slp_this = seqloc_list;
    while (slp_this != NULL)
    {
      slp_next = slp_this->next;
      slp_this->next = NULL;
      slp_this = SeqLocFree (slp_this);
      slp_this = slp_next;
    }
    return;
  }
  
  if (feature_details_data == NULL || seqloc_list == NULL)
  {
    /* Free data and loclist */
    feature_details_data = BatchApplyFeatureDetailsFree (feature_details_data);
    slp_this = seqloc_list;
    while (slp_this != NULL)
    {
      slp_next = slp_this->next;
      slp_this->next = NULL;
      slp_this = SeqLocFree (slp_this);
      slp_this = slp_next;
    }
    return;
  }
  
  if (aafdp->feattype == ADD_IMP)
  {
    gbqual_list = GetImpFeatureGBQuals (feature_details_data);
  }

  slp_this = seqloc_list;
  while (slp_this != NULL)
  {
    slp_next = slp_this->next;
    slp_this->next = NULL;
    sfp = NULL;
    
    bsp = BioseqFindFromSeqLoc (slp_this);
    sep = SeqMgrGetSeqEntryForData (bsp);
  
    if (sep != NULL)
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL)
      { 
        if (aafdp->feattype == ADD_CDS)
        {    
          sfp = ApplyOneCodingRegion (feature_details_data, bsp, sep, nsep,
                                      aafdp->entityID, suppressDups);
        }
        else if (aafdp->feattype == ADD_RRNA)
        {
          sfp = ApplyOneRNA (feature_details_data, sep, nsep,
                             aafdp->entityID, suppressDups);
        }
        else if (aafdp->feattype == ADD_IMP)
        {
          sfp = ApplyOtherFeature (feature_details_data, nsep, gbqual_list);
        }
      }
    }
    if (sfp != NULL)
    {
      /* set location */
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp_this;
      CheckSeqLocForPartial (slp_this, &partial5, &partial3);
      sfp->partial = partial5 || partial3;
      
      if (aafdp->feattype == ADD_CDS)
      {
        SetFrameForCodingRegion (sfp, bsp, feature_details_data,
                                 &errcount, &ambigList);
        AddProductForCDS (sfp, feature_details_data, sep, nsep, aafdp->entityID);
      }
      
    }
    
    if (sfp != NULL
        || (aafdp->feattype == ADD_IMP
            && feature_details_data->featdef_choice == FEATDEF_GENE))
    {
      if (! StringHasNoText (feature_details_data->geneName) || ! StringHasNoText (feature_details_data->geneDesc)) 
      {
        /* Create a Gene ref feature on the nuc seq or segment */
        /* we can only create a feature where the sep->choice is 1 */
        if (sep->choice == 1)
        {
          gene_sep = sep;
        }
        else
        {
          /* if we have added a product, nsep may no longer point to the nucseq entry */
          gene_sep = FindNucSeqEntry (sep);
        }
      
        if (aafdp->entityID > 0 
            && suppressDups
            && AlreadyHasFeatOrDesc (gene_sep, SEQFEAT_GENE, 0, 0))
        {
          /* do not create */
        }
        else
        {
          CheckSeqLocForPartial (slp_this, &partial5, &partial3);
        
          if (sfp == NULL)
          {
            gene_slp = slp_this;
          }
          else
          {
            gene_slp = SeqLocMerge (bsp, slp_this, NULL, TRUE, FALSE, FALSE);
            SetSeqLocPartial (gene_slp, partial5, partial3);
          }
          gene_sfp = ApplyGene (feature_details_data->geneName,
                                feature_details_data,
                                gene_sep, sfp, gene_slp);
          gene_sfp->partial = partial5 || partial3;
                                
          if (sfp == NULL)
          {
            sfp = gene_sfp;
          }
        } 
      }
    }
    
    if (sfp == NULL)
    {
      /* remove location that will not be used */
      slp_this = SeqLocFree (slp_this);
    }
    else
    {
      /* add comment */
      if (!StringHasNoText (feature_details_data->featcomment))
      {
        sfp->comment = StringSave (feature_details_data->featcomment);
      }
    }
    
    slp_this = slp_next;
  }
  
  feature_details_data = BatchApplyFeatureDetailsFree (feature_details_data);
  
  /* free gbqual_list */
  gbqual_list = GBQualFree (gbqual_list);

  ObjMgrSetDirtyFlag (aafdp->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, aafdp->entityID, 0, 0);
    
  Remove (aafdp->form);
}

extern void ApplyFeatureToAlignment (Uint2 entityID, SeqAlignPtr salp, SeqLocPtr slp, Int4 feattype)
{
  ApplyAlignmentFeatureDlgPtr aafdp;
  WindoW                      w;
  GrouP                       h, c;
  Boolean                     nucsOK = TRUE;  /* LATER, figure out whether alignment */
  Boolean                     protsOK = TRUE; /* contains nucs or prots */
  ButtoN                      b;
  
  aafdp = (ApplyAlignmentFeatureDlgPtr) MemNew (sizeof (ApplyAlignmentFeatureDlgData));
  if (aafdp == NULL) return;
  
  aafdp->sep = GetTopSeqEntryForEntityID (entityID);
  if (aafdp->sep == NULL)
  {
    aafdp = MemFree (aafdp);
    return;
  }
  aafdp->entityID = entityID;
  aafdp->feattype = feattype;
  aafdp->salp = salp;
    
  w = FixedWindow (-50, -33, -10, -10, "Apply Features to Alignment", StdCloseWindowProc); 
  SetObjectExtra (w, aafdp, StdCleanupExtraProc); 
  SetObjectExtra (w, aafdp, NULL);
  aafdp->form = (ForM) w;
  aafdp->input_entityID = entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  aafdp->location_dlg = CreateIntervalEditorDialogExEx (h, "Location",
                                                        4, 2, aafdp->sep,
                                                        nucsOK, protsOK,
                                                        TRUE, TRUE, TRUE, NULL, NULL,
                                                        TRUE, FALSE, NULL, NULL,
                                                        TRUE);
  PointerToDialog (aafdp->location_dlg, slp);                                                      
  
  aafdp->feature_details = BatchApplyFeatureDetailsDialog (h, feattype);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoApplyFeatureToAlignment);
  SetObjectExtra (b, aafdp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) aafdp->location_dlg,
                              (HANDLE) aafdp->feature_details,
                              (HANDLE) c,
                               NULL);

  Show (w);
}

extern SeqAnnotPtr GetSeqAnnotForAlignment (SeqAlignPtr sap)
{
  SeqAnnotPtr  sanp = NULL;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  Boolean      found = FALSE;
  
  if (sap == NULL)
  {
    return NULL;
  }
  
  if (sap->idx.parenttype == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) sap->idx.parentptr;
    if (bsp != NULL)
    {
      sanp = bsp->annot;
    }
  }
  else if (sap->idx.parenttype == OBJ_BIOSEQSET)
  {
    bssp = (BioseqSetPtr) sap->idx.parentptr;
    if (bssp != NULL)
    {
      sanp = bssp->annot;
    }
  }
  else if (sap->idx.parenttype == OBJ_SEQANNOT)
  {
    sanp = (SeqAnnotPtr) sap->idx.parentptr;
  }
  while (sanp != NULL && !found)
  {
    if (sanp->type == 2 && sanp->data == sap)
    {
      found = TRUE;
    }
    else
    {
      sanp = sanp->next;
    }
  }
  return sanp;
}

extern void ConvertPairwiseToMultipleAlignment (SeqAlignPtr sap)
{
  SeqAlignPtr salp, salp_next;
  
  salp = sap;
  while (salp != NULL)
  {
     if (salp->saip != NULL)
     {
        SeqAlignIndexFree(salp->saip);
        salp->saip = NULL;
     }
     salp = salp->next;
  }
  AlnMgr2IndexSeqAlign (sap);

  salp = AlnMgr2GetSubAlign (sap, 0, -1, 0, FALSE);
  if (salp == NULL) return;
  AlnMgr2IndexSeqAlign (salp);

  if (sap->idx.prevlink != NULL) {
    *(sap->idx.prevlink) = (Pointer) salp;
  }

  while (sap != NULL) {
    salp_next = sap->next;
    sap->next = NULL;
    SeqAlignFree (sap);
    sap = salp_next;
  }
}
