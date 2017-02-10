#ifndef CONNECT___NCBI_BUFFER__H
#define CONNECT___NCBI_BUFFER__H

/*  $Id: ncbi_buffer.h,v 6.11 2004/10/27 18:43:45 lavr Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Denis Vakatov
 *
 * File Description:
 *   Memory-resident FIFO storage area (to be used e.g. in I/O buffering)
 *
 * Handle:  BUF
 *
 * Functions:
 *   BUF_SetChunkSize
 *   BUF_Size
 *   BUF_Prepend
 *   BUF_Append
 *   BUF_Write
 *   BUF_PushBack
 *   BUF_Peek
 *   BUF_PeekAt
 *   BUF_Read
 *   BUF_Destroy
 *
 */

#if defined(NCBIBUF__H)
#  error "<ncbibuf.h> and <ncbi_buffer.h> must never be #include'd together"
#endif

#include <connect/connect_export.h>
#include <stddef.h>     /* ...to define "size_t"... */


/** @addtogroup BuffServices
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


struct BUF_tag;
typedef struct BUF_tag* BUF;  /* handle of a buffer */


/*!
 * Set minimal size of a buffer memory chunk.
 * Return the actually set chunk size on success;  zero on error
 * NOTE:  if "*pBuf" == NULL then create it
 *        if "chunk_size" is passed 0 then set it to BUF_DEF_CHUNK_SIZE
 */
#define BUF_DEF_CHUNK_SIZE 1024
extern NCBI_XCONNECT_EXPORT size_t BUF_SetChunkSize
(BUF*        pBuf,
 size_t      chunk_size
 );


/*!
 * Return the number of bytes stored in "buf".
 * NOTE: return 0 if "buf" == NULL
 */
extern NCBI_XCONNECT_EXPORT size_t BUF_Size(BUF buf);


/*!
 * Prepend a block of data (of the specified size) to the
 * beginning of the buffer (to be read first).  Note that unlike
 * BUF_Pushback(), in this call the data is not copied into the buffer
 * but instead is just linked in from the original location.
 * Return non-zero (true) if succeeded, zero (false) if failed.
 */
extern NCBI_XCONNECT_EXPORT int/*bool*/ BUF_Prepend
(BUF*        pBuf,
 const void* data,
 size_t      size
);


/*!
 * Append a block of data (of the specified size) past the end
 * of the buffer (to be read last).  Note that unlike
 * BUF_Write(), in this call the data is not copied to the buffer
 * but instead is just linked in from the original location.
 * Return non-zero (true) if succeeded, zero (false) if failed.
 */
extern NCBI_XCONNECT_EXPORT int/*bool*/ BUF_Append
(BUF*        pBuf,
 const void* data,
 size_t      size
 );


/*!
 * Add new data to the end of "*pBuf" (to be read last).
 * On error (failed memory allocation), return zero value.
 * NOTE:  if "*pBuf" == NULL then create it.
 */
extern NCBI_XCONNECT_EXPORT /*bool*/int BUF_Write
(BUF*        pBuf,
 const void* data,
 size_t      size
 );


/*!
 * Write the data to the very beginning of "*pBuf" (to be read first).
 * On error (failed memory allocation), return zero value.
 * NOTE:  if "*pBuf" == NULL then create it.
 */
extern NCBI_XCONNECT_EXPORT /*bool*/int BUF_PushBack
(BUF*        pBuf,
 const void* data,
 size_t      size
 );


/*!
 * Equivalent to "BUF_PeekAt(buf, 0, data, size)", see description below.
 */
extern NCBI_XCONNECT_EXPORT size_t BUF_Peek
(BUF         buf,
 void*       data,
 size_t      size
 );


/*!
 * Copy up to "size" bytes stored in "buf" (starting at position "pos")
 * to "data".
 * Return the # of copied bytes (can be less than "size").
 * Return zero and do nothing if "buf" is NULL or "pos" >= BUF_Size(buf).
 * Do nothing and return min(BUF_Size(buf)-pos, size) if "data" is NULL.
 */
extern NCBI_XCONNECT_EXPORT size_t BUF_PeekAt
(BUF         buf,
 size_t      pos,
 void*       data,
 size_t      size
 );


/*!
 * Copy up to "size" bytes stored in "buf" to "data" and remove
 * copied data from the "buf".
 * Return the # of copied-and/or-removed bytes(can be less than "size")
 * NOTE: if "buf"  == NULL then do nothing and return 0
 *       if "data" == NULL then do not copy data anywhere(still, remove it)
 */
extern NCBI_XCONNECT_EXPORT size_t BUF_Read
(BUF         buf,
 void*       data,
 size_t      size
 );


/*!
 * Make the buffer empty.
 */
extern NCBI_XCONNECT_EXPORT void BUF_Erase(BUF buf);


/*!
 * Destroy all internal data.
 * NOTE: do nothing if "buf" == NULL
 */
extern NCBI_XCONNECT_EXPORT void BUF_Destroy(BUF buf);


#ifdef __cplusplus
}
#endif


/* @} */


/*
 * ---------------------------------------------------------------------------
 * $Log: ncbi_buffer.h,v $
 * Revision 6.11  2004/10/27 18:43:45  lavr
 * +BUF_Erase()
 *
 * Revision 6.10  2004/10/27 18:09:57  lavr
 * +BUF_Prepend(), +BUF_Append()
 *
 * Revision 6.9  2003/04/10 12:52:15  siyan
 * Changed group name for doxygen
 *
 * Revision 6.8  2003/04/09 17:58:40  siyan
 * Added doxygen support
 *
 * Revision 6.7  2003/01/08 01:59:32  lavr
 * DLL-ize CONNECT library for MSVC (add NCBI_XCONNECT_EXPORT)
 *
 * Revision 6.6  2002/09/19 17:59:53  lavr
 * Header file guard macro changed; log moved to the end
 *
 * Revision 6.5  2001/04/23 22:20:26  vakatov
 * BUF_PeekAt() -- special case for "data" == NULL
 *
 * Revision 6.4  2001/04/23 18:07:19  vakatov
 * + BUF_PeekAt()
 *
 * Revision 6.3  2000/02/23 22:33:37  vakatov
 * Can work both "standalone" and as a part of NCBI C++ or C toolkits
 *
 * Revision 6.2  1999/10/12 16:30:10  vakatov
 * include <string.h> to define "size_t"
 *
 * Revision 6.1  1999/08/17 19:45:22  vakatov
 * Moved all real code from NCBIBUF to NCBI_BUFFER;  the code has been cleaned
 * from the NCBI C toolkit specific types and API calls.
 * NCBIBUF module still exists for the backward compatibility -- it
 * provides old NCBI-wise interface.
 *
 * ===========================================================================
 */

#endif /* CONNECT___NCBI_BUFFER__H */
