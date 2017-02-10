/*  $Id: test_ncbi_http_get.c,v 6.14 2005/04/20 15:51:19 lavr Exp $
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
 *   Retrieve a Web document via HTTP
 *
 */

#include "../ncbi_ansi_ext.h"
#include "../ncbi_priv.h"
#include <connect/ncbi_http_connector.h>
#include <connect/ncbi_util.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#ifdef NCBI_OS_UNIX
#  include <unistd.h>
#endif
/* This header must go last */
#include "test_assert.h"


int main(int argc, char* argv[])
{
    static const STimeout s_ZeroTmo = {0, 0};
    CONNECTOR     connector;
    SConnNetInfo* net_info;
    char          blk[512];
    EIO_Status    status;
    THCC_Flags    flags;
    CONN          conn;
    FILE*         fp;
    time_t        t;
    size_t        n;
    char*         s;

    CORE_SetLOGFormatFlags(fLOG_None          | fLOG_Level   |
                           fLOG_OmitNoteLevel | fLOG_DateTime);
    CORE_SetLOGFILE(stderr, 0/*false*/);

    if (argc < 2 || !*argv[1])
        CORE_LOG(eLOG_Fatal, "URL has to be supplied on the command line");
    if (argc > 3)
        CORE_LOG(eLOG_Fatal, "Command cannot take more than 2 arguments");
    if (argc == 3) {
        fp = strcmp(argv[2], "-") == 0 ? stdin : fopen(argv[2], "rb");
        if (!fp) {
            CORE_LOGF_ERRNO(eLOG_Error, errno, ("Cannot open \"%s\"",
                                                argv[2] ? argv[2] : ""));
        }
    } else
        fp = 0;

    CORE_LOG(eLOG_Note, "Creating network info structure");
    if (!(net_info = ConnNetInfo_Create(0)))
        CORE_LOG(eLOG_Fatal, "Cannot create network info structure");
    if ((s = getenv("CONN_TIMEOUT"))  &&  strcmp(s, "0") == 0) {
        memcpy(&net_info->tmo, &s_ZeroTmo, sizeof(s_ZeroTmo));
        net_info->timeout = &net_info->tmo;
    }

    CORE_LOGF(eLOG_Note,
              ("Parsing URL \"%s\" into network info structure", argv[1]));
    if (!ConnNetInfo_ParseURL(net_info, argv[1]))
        CORE_LOG(eLOG_Fatal, "URL parsing failed");

    if ((s = getenv("CONN_RECONNECT")) != 0 &&
        (strcasecmp(s, "TRUE") == 0 || strcasecmp(s, "1") == 0)) {
        flags = fHCC_AutoReconnect;
        CORE_LOG(eLOG_Note, "Reconnect mode will be used");
    } else
        flags = 0;

    CORE_LOG(eLOG_Note, "Creating HTTP connector");
    if (!(connector = HTTP_CreateConnector(net_info, 0, flags)))
        CORE_LOG(eLOG_Fatal, "Cannot create HTTP connector");

    CORE_LOG(eLOG_Note, "Creating connection");
    if (CONN_Create(connector, &conn) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Cannot create connection");
    CONN_SetTimeout(conn, eIO_Open,      net_info->timeout);
    CONN_SetTimeout(conn, eIO_ReadWrite, net_info->timeout);
    while (fp  &&  !feof(fp)) {
        n = fread(blk, 1, sizeof(blk), fp);
        status = CONN_Write(conn, blk, n, &n, eIO_WritePersist);
        if (status != eIO_Success) {
            CORE_LOGF(eLOG_Fatal, ("Write error: %s", IO_StatusStr(status)));
        }
    }
    if (fp)
        fclose(fp);

    t = time(0);
    for (;;) {
        status = CONN_Wait(conn, eIO_Read, net_info->timeout);
        if (status != eIO_Success) {
            if (status == eIO_Closed)
                break;
            if ((unsigned long)(time(0) - t) > 30)
                CORE_LOG(eLOG_Fatal, "Timed out");
#ifdef NCBI_OS_UNIX
            usleep(500);
#endif
            continue;
        }

        status = CONN_Read(conn, blk, sizeof(blk), &n, eIO_ReadPlain);
        if (status != eIO_Success && status != eIO_Closed)
            CORE_LOGF(eLOG_Fatal, ("Read error: %s", IO_StatusStr(status)));
        if (n) {
            fwrite(blk, 1, n, stdout);
            fflush(stdout);
        } else if (status == eIO_Closed) {
            break;
        } else
            CORE_LOG(eLOG_Fatal, "Empty read");
    }

    ConnNetInfo_Destroy(net_info);
    CORE_LOG(eLOG_Note, "Closing connection");
    CONN_Close(conn);

    CORE_LOG(eLOG_Note, "Completed");
    CORE_SetLOG(0);
    return 0;
}


/*
 * --------------------------------------------------------------------------
 * $Log: test_ncbi_http_get.c,v $
 * Revision 6.14  2005/04/20 15:51:19  lavr
 * Allow addtl file argument to attach as a body ('-' for stdin)
 *
 * Revision 6.13  2005/04/19 16:37:52  lavr
 * Allow HTTP method to be overridden from the environment
 *
 * Revision 6.12  2003/09/30 20:59:39  lavr
 * Fix typo in previous log message
 *
 * Revision 6.11  2003/09/30 20:57:15  lavr
 * Allow to set zero timeout via environment
 *
 * Revision 6.10  2003/05/20 23:52:01  lavr
 * Explicit cast "time_t"->"unsigned" to avoid GCC warning
 *
 * Revision 6.9  2003/05/19 16:58:37  lavr
 * Modified to use {0,0} timeouts in CONN_Wait() and use app timeout handling
 *
 * Revision 6.8  2003/05/14 03:58:43  lavr
 * Match changes in respective APIs of the tests
 *
 * Revision 6.7  2002/10/28 15:47:06  lavr
 * Use "ncbi_ansi_ext.h" privately
 *
 * Revision 6.6  2002/08/14 14:42:46  lavr
 * Use fwrite() instead of printf() when printing out the data fetched
 *
 * Revision 6.5  2002/08/07 16:38:08  lavr
 * EIO_ReadMethod enums changed accordingly; log moved to end
 *
 * Revision 6.4  2002/03/22 19:47:23  lavr
 * Test_assert.h made last among the include files
 *
 * Revision 6.3  2002/02/05 21:45:55  lavr
 * Included header files rearranged
 *
 * Revision 6.2  2002/01/16 21:23:15  vakatov
 * Utilize header "test_assert.h" to switch on ASSERTs in the Release mode too
 *
 * Revision 6.1  2001/12/30 19:42:06  lavr
 * Initial revision
 *
 * ==========================================================================
 */
