/*  $Id: test_ncbi_connutil_hit.c,v 6.14 2006/04/20 14:01:43 lavr Exp $
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
 *   Test "ncbi_connutil.c"::URL_Connect()
 *
 */

#include <connect/ncbi_connutil.h>
#include <connect/ncbi_util.h>
#include <stdlib.h>
/* This header must go last */
#include "test_assert.h"


int main(int argc, char** argv)
{
    /* Prepare to connect:  parse and check cmd.-line args, etc. */
    const char*    host        = (argc > 1) ? argv[1] : "";
    unsigned short port        = (argc > 2) ? (unsigned short) atoi(argv[2]):0;
    const char*    path        = (argc > 3) ? argv[3] : "";
    const char*    args        = (argc > 4) ? argv[4] : "";
    const char*    inp_file    = (argc > 5) ? argv[5] : "";
    const char*    user_header = (argc > 6) ? argv[6] : "";

    size_t   content_length;
    STimeout timeout;

    SOCK sock;
    EIO_Status status;
    char buffer[10000];

    CORE_SetLOGFormatFlags(fLOG_None          | fLOG_Level   |
                           fLOG_OmitNoteLevel | fLOG_DateTime);
    CORE_SetLOGFILE(stderr, 0/*false*/);

    fprintf(stderr, "Running...\n"
            "  Executable:      '%s'\n"
            "  URL host:        '%s'\n"
            "  URL port:         %hu\n"
            "  URL path:        '%s'\n"
            "  URL args:        '%s'\n"
            "  Input data file: '%s'\n"
            "  User header:     '%s'\n"
            " Reply(if any) from the hit URL goes to the standard output.\n\n",
            argv[0],
            host, (unsigned short) port, path, args, inp_file, user_header);

    if ( argc < 4 ) {
        fprintf(stderr,
                "Usage:   %s host port path args inp_file [user_header]\n"
                "Example: %s ............\n"
                "\nTwo few arguments.\n",
                argv[0], argv[0]);
        return 1;
    }

    {{
        FILE *fp = fopen(inp_file, "rb");
        long offset;

        if ( !fp ) {
            fprintf(stderr, "Non-existent file '%s'\n", inp_file);
            return 2;
        }
        if ( fseek(fp, 0, SEEK_END) != 0  ||  (offset = ftell(fp)) < 0 ) {
            fprintf(stderr, "Cannot obtain size of file '%s'\n", inp_file);
            return 2;
        }
        fclose(fp);
        content_length = (size_t) offset;
    }}

    timeout.sec  = 10;
    timeout.usec = 0;
    
    /* Connect */
    sock = URL_Connect(host, port, path, args,
                       eReqMethod_Any, content_length,
                       &timeout, &timeout, user_header, 1/*true*/, eDefault);
    if ( !sock )
        return 3;
    
    {{ /* Pump data from the input file to socket */
        FILE* fp = fopen(inp_file, "rb");
        if ( !fp ) {
            fprintf(stderr, "Cannot open file '%s' for reading\n", inp_file);
            return 4;
        }

        for (;;) {
            size_t n_written;
            size_t n_read = fread(buffer, 1, sizeof(buffer), fp);
            if ( n_read <= 0 ) {
                if ( content_length ) {
                    fprintf(stderr,
                            "Cannot read last %lu bytes from file '%s'\n",
                            (unsigned long) content_length, inp_file);
                    return 5;
                }
                break;
            }

            assert(content_length >= n_read);
            content_length -= n_read;
            status = SOCK_Write(sock, buffer, n_read,
                                &n_written, eIO_WritePersist);
            if ( status != eIO_Success ) {
                fprintf(stderr, "Error writing to socket(%s)\n",
                        IO_StatusStr(status));
                return 6;
            }
        }

        fclose(fp);
    }}

    /* Read reply from socket, write it to STDOUT */
    {{
        size_t n_read;
        for (;;) {
            status = SOCK_Read(sock, buffer, sizeof(buffer), &n_read,
                               eIO_ReadPlain);
            if (status != eIO_Success)
                break;

            fwrite(buffer, 1, n_read, stdout);
        }

        if ( status != eIO_Closed ) {
            fprintf(stderr,
                    "Error occurred after reading %ld bytes from socket(%s)\n",
                    (long) content_length, IO_StatusStr(status));
        }
        fprintf(stdout, "\n");
    }}

    /* Success:  close the socket, cleanup, and exit */
    SOCK_Close(sock);
    CORE_SetLOG(0);

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * $Log: test_ncbi_connutil_hit.c,v $
 * Revision 6.14  2006/04/20 14:01:43  lavr
 * Cleanup to demonstrate no leaks
 *
 * Revision 6.13  2005/07/22 16:09:43  lavr
 * Add severity and timestamp to log message format
 *
 * Revision 6.12  2002/08/12 15:10:43  lavr
 * Use persistent SOCK_Write()
 *
 * Revision 6.11  2002/08/07 16:38:08  lavr
 * EIO_ReadMethod enums changed accordingly; log moved to end
 *
 * Revision 6.10  2002/04/15 19:21:44  lavr
 * +#include "../test/test_assert.h"
 *
 * Revision 6.9  2002/03/22 19:48:56  lavr
 * Removed <stdio.h>: included from ncbi_util.h or ncbi_priv.h
 *
 * Revision 6.8  2002/02/05 21:45:55  lavr
 * Included header files rearranged
 *
 * Revision 6.7  2001/03/02 20:03:57  lavr
 * Typos fixed
 *
 * Revision 6.6  2001/01/23 23:21:21  lavr
 * Added new argument to URL_Connect
 *
 * Revision 6.5  2001/01/08 22:42:14  lavr
 * eRequestMethodAny -> eRequestMethod_Any
 *
 * Revision 6.4  2000/12/29 18:24:40  lavr
 * File size is now discovered without use of stat call.
 *
 * Revision 6.3  2000/09/27 13:49:29  lavr
 * URL_Connect args adjusted
 *
 * Revision 6.2  2000/03/29 17:21:47  vakatov
 * + CORE_SetLOG(0) at the program end.
 *
 * Revision 6.1  2000/03/24 22:53:38  vakatov
 * Initial revision
 *
 * ===========================================================================
 */
