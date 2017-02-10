/*  $Id: asn_connection.c,v 1.8 2004/11/23 16:15:29 lavr Exp $
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
 *     Build C Toolkit ASN streams on top of CONN (connection).
 *
 */

#include <ctools/asn_connection.h>
#include "../connect/ncbi_priv.h"


#ifdef __cplusplus
extern "C" {
    static Int2 LIBCALLBACK s_AsnRead(Pointer p, CharPtr buff, Uint2 len);
    static Int2 LIBCALLBACK s_AsnWrite(Pointer p, CharPtr buff, Uint2 len);
    static void s_CloseAsnConn(CONN conn, ECONN_Callback type, void* data);
}
#endif


static Int2 LIBCALLBACK s_AsnRead(Pointer p, CharPtr buff, Uint2 len)
{
    size_t n_read;
    CONN_Read((CONN) p, buff, len, &n_read, eIO_ReadPlain);
    return (Int2) n_read;
}


static Int2 LIBCALLBACK s_AsnWrite(Pointer p, CharPtr buff, Uint2 len)
{
    size_t n_written;
    CONN_Write((CONN) p, buff, len, &n_written, eIO_WritePersist);
    return (Int2) n_written;
}


struct SAsnConn_Cbdata {
    SCONN_Callback cb;
    AsnIoPtr       ptr;
};


static void s_CloseAsnConn(CONN conn, ECONN_Callback type, void* data)
{
    struct SAsnConn_Cbdata* cbdata = (struct SAsnConn_Cbdata*) data;

    assert(type == eCONN_OnClose && cbdata && cbdata->ptr);
    CONN_SetCallback(conn, type, &cbdata->cb, 0); 
    AsnIoFree(cbdata->ptr, 0/*not a file - don't close*/);
    if ( cbdata->cb.func )
        (*cbdata->cb.func)(conn, type, cbdata->cb.data);
    free(cbdata);
}


static void s_SetAsnConn_CloseCb(CONN conn, AsnIoPtr ptr)
{
    struct SAsnConn_Cbdata* cbdata = (struct SAsnConn_Cbdata*)
        malloc(sizeof(*cbdata));

    assert( ptr );
    if ( cbdata ) {
        SCONN_Callback cb;
        cbdata->ptr = ptr;
        cb.func = s_CloseAsnConn;
        cb.data = cbdata;
        CONN_SetCallback(conn, eCONN_OnClose, &cb, &cbdata->cb);
    } else
        CORE_LOG(eLOG_Error,
                 "Cannot create cleanup callback for ASN conn-based stream");
}


AsnIoPtr CreateAsnConn(CONN               conn,
                       EAsnConn_Direction direction,
                       EAsnConn_Format    fmt)
{
    AsnIoPtr ptr;
    int flags;

    switch (fmt) {
    case eAsnConn_Binary:
        flags = ASNIO_BIN;
        break;
    case eAsnConn_Text:
        flags = ASNIO_TEXT;
        break;
    default:
        return 0;
    }

    switch (direction) {
    case eAsnConn_Input:
        ptr = AsnIoNew(flags | ASNIO_IN, 0, (void*) conn, s_AsnRead, 0);
        break;
    case eAsnConn_Output:
        ptr = AsnIoNew(flags | ASNIO_OUT, 0, (void*) conn, 0, s_AsnWrite);
        break;
    default:
        return 0;
    }

    if (ptr)
        s_SetAsnConn_CloseCb(conn, ptr);
    return ptr;
}


CONN CreateAsnConn_ServiceEx(const char*           service,
                             EAsnConn_Format       input_fmt,
                             AsnIoPtr*             input,
                             EAsnConn_Format       output_fmt,
                             AsnIoPtr*             output,
                             TSERV_Type            type,
                             const SConnNetInfo*   net_info,
                             const SSERVICE_Extra* params)
{
    CONN conn;
    CONNECTOR c = SERVICE_CreateConnectorEx(service, type, net_info, params);
    if (!c || CONN_Create(c, &conn) != eIO_Success)
        return 0/*failed*/;
    assert(conn);

    if (input)
        *input = CreateAsnConn(conn, eAsnConn_Input, input_fmt);

    if (output)
        *output = CreateAsnConn(conn, eAsnConn_Output, output_fmt);

    return conn;
}


CONN CreateAsnConn_Service(const char*     service,
                           EAsnConn_Format input_fmt,
                           AsnIoPtr*       input,
                           EAsnConn_Format output_fmt,
                           AsnIoPtr*       output)
{
    return CreateAsnConn_ServiceEx(service, input_fmt, input,
                                   output_fmt, output, fSERV_Any, 0, 0);
}


/*
 * --------------------------------------------------------------------------
 * $Log: asn_connection.c,v $
 * Revision 1.8  2004/11/23 16:15:29  lavr
 * -<connect/ncbi_service_connector.h> now included from <ctools/asn_connectionh>
 *
 * Revision 1.7  2004/03/23 02:29:19  lavr
 * Remove redundant inits in I/O callbacks
 *
 * Revision 1.6  2004/02/23 15:23:50  lavr
 * New (last) parameter "how" added in CONN_Write() API call
 *
 * Revision 1.5  2003/11/13 16:01:31  lavr
 * Included headers revised
 *
 * Revision 1.4  2002/08/07 16:38:18  lavr
 * EIO_ReadMethod enums changed accordingly; log moved to end
 *
 * Revision 1.3  2001/09/24 20:32:44  lavr
 * +SSERVICE_Extra* parameter in CreateAsnConn_ServiceEx()
 *
 * Revision 1.2  2001/06/29 14:04:17  lavr
 * Close callback removes itself when called; fixes not to add callbacks
 * in CreateAsnConn_ServiceEx, because they are already there!
 *
 * Revision 1.1  2001/06/28 21:58:23  lavr
 * Initial revision
 *
 * ==========================================================================
 */
