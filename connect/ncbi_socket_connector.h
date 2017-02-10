#ifndef CONNECT___NCBI_SOCKET_CONNECTOR__H
#define CONNECT___NCBI_SOCKET_CONNECTOR__H

/*  $Id: ncbi_socket_connector.h,v 6.9 2003/05/14 03:48:51 lavr Exp $
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
 *   Implement CONNECTOR for a network socket(based on the NCBI "SOCK").
 *
 *   See in "connectr.h" for the detailed specification of the underlying
 *   connector("CONNECTOR", "SConnectorTag") methods and structures.
 *
 */

#include <connect/ncbi_connector.h>
#include <connect/ncbi_socket.h>


/** @addtogroup Connectors
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/* This is equivalent to SOCK_CreateConnectorEx(host, port, max_try, 0,0,0),
 * see below.
 */
extern NCBI_XCONNECT_EXPORT CONNECTOR SOCK_CreateConnector
(const char*    host,   /* server:  host                            */
 unsigned short port,   /* server:  service port                    */
 unsigned int   max_try /* max.number of attempts to establish conn */
 );


typedef enum {
    eSCC_DebugPrintout  = 1,
    eSCC_SetReadOnWrite = 2
} ESCC_Flags;
typedef unsigned int TSCC_Flags;  /* bitwise OR of "ESCC_Flags */


/* Create new CONNECTOR structure to handle connection to a socket.
 * Make up to "max_try" attempts to connect to the "host:port" before
 * giving up.
 * On successful connect, send the first "init_size" bytes from buffer
 * "init_data"(can be NULL -- then send nothing) to the opened connection.
 * NOTE:  the connector makes (and then uses) its own copy of the "init_data".
 * Return NULL on error.
 */
extern NCBI_XCONNECT_EXPORT CONNECTOR SOCK_CreateConnectorEx
(const char*    host,      /* server:  host                                  */
 unsigned short port,      /* server:  service port                          */
 unsigned int   max_try,   /* max.number of attempts to establish connection */
 const void*    init_data, /* data to send to server on connect              */
 size_t         init_size, /* size of the "init_data" buffer                 */
 TSCC_Flags     flags      /* bitwise OR of additional flags: see above      */
 );


/* Equivalent to SOCK_CreateConnectorOnTopEx(sock, max_try, 0,0,0),
 * see below.
 */
extern NCBI_XCONNECT_EXPORT CONNECTOR SOCK_CreateConnectorOnTop
(SOCK         sock,   /* socket object                                       */
 unsigned int max_try /* max.number of tries to re-establish the link if lost*/
 );


/* Create new CONNECTOR structure on top of existing socket object (SOCK),
 * acquiring the ownership of the socket, and overriding all timeouts
 * that might have been set already in it. Timeout values will be taken from
 * connection (CONN), after the connector is used in CONN_Create() call.
 * Please note that this function revokes all ownership of the socket, and
 * further assumes the passed socket is for the sole use of the connector.
 * A socket obtained as a result of accepting connection on a listening socket
 * (aka server-side socket) is not allowed to have reconnects (max_try = 0).
 */
extern NCBI_XCONNECT_EXPORT CONNECTOR SOCK_CreateConnectorOnTopEx
(SOCK         sock,      /* socket object                                    */
 unsigned int max_try,   /* max.# of tries to reconnect if disconnected      */
 TSCC_Flags   flags      /* bitwise OR of additional flags: see above        */
 );


#ifdef __cplusplus
}  /* extern "C" */
#endif


/* @} */


/*
 * --------------------------------------------------------------------------
 * $Log: ncbi_socket_connector.h,v $
 * Revision 6.9  2003/05/14 03:48:51  lavr
 * SOCKET_CreateConnectorOnTopEx(): init_data and init_size parameters removed
 *
 * Revision 6.8  2003/04/09 19:05:57  siyan
 * Added doxygen support
 *
 * Revision 6.7  2003/01/08 01:59:33  lavr
 * DLL-ize CONNECT library for MSVC (add NCBI_XCONNECT_EXPORT)
 *
 * Revision 6.6  2002/09/19 18:05:34  lavr
 * Header file guard macro changed; log moved to end
 *
 * Revision 6.5  2002/08/07 16:31:40  lavr
 * Added eSCC_SetReadOnWrite to constructor flags; log moved to end
 *
 * Revision 6.4  2002/01/12 22:15:28  lavr
 * Fixed function description to use not former but current parameter names
 *
 * Revision 6.3  2001/12/04 15:54:29  lavr
 * +SOCK_CreateConnectorOnTop(), +SOCK_CreateConnectorOnTopEx()
 *
 * Revision 6.2  2001/01/23 23:09:18  lavr
 * Flags added to 'Ex' constructor
 *
 * Revision 6.1  2000/04/07 20:05:37  vakatov
 * Initial revision
 *
 * ==========================================================================
 */

#endif /* CONNECT___NCBI_SOCKET_CONNECTOR__H */
