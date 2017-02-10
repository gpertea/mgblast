#ifndef CONNECT___NCBI_SOCKET_UNIX__H
#define CONNECT___NCBI_SOCKET_UNIX__H

/*  $Id: ncbi_socket_unix.h,v 1.1 2004/10/26 14:44:44 lavr Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   TCP/IP socket API extension for UNIX
 *
 */

#include <connect/ncbi_socket.h>


/** @addtogroup Sockets
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


extern NCBI_XCONNECT_EXPORT EIO_Status LSOCK_CreateUNIX
(const char*    path,    /* [in]  filename of the named socket to create */
 unsigned short backlog, /* [in]  maximal # of pending connections       */
 LSOCK*         lsock,   /* [out] handle of the created listening socket */
 ESwitch        log      /* [in]  whether to do logging on this socket   */
 );


extern NCBI_XCONNECT_EXPORT EIO_Status SOCK_CreateUNIX
(const char*     file,     /* [in]  filename of the UNIX socket to connect to*/
 const STimeout* timeout,  /* [in]  connection timeout (infinite if NULL)    */
 SOCK*           sock,     /* [out] handle of the created socket             */
 const void*     init_data,/* [in]  initial output data segment (may be NULL)*/
 size_t          init_size,/* [in]  size of initial data segment (may be 0)  */
 ESwitch         log       /* [in]  whether to do logging on this socket     */
 );


#ifdef __cplusplus
} /* extern "C" */
#endif


/* @} */


/*
 * ---------------------------------------------------------------------------
 * $Log: ncbi_socket_unix.h,v $
 * Revision 1.1  2004/10/26 14:44:44  lavr
 * Initial revision
 *
 * ===========================================================================
 */

#endif /* CONNECT___NCBI_SOCKET_UNIX__H */
