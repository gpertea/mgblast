/*  $Id: test_ncbi_buffer.c,v 6.10 2003/02/25 16:26:58 lavr Exp $
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
 *   Test suite for "ncbi_buffer.[ch]", the memory-resident FIFO storage area
 *
 */

#include <connect/ncbi_buffer.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* This header must go last */
#include "test_assert.h"


static unsigned s_Rand(void)
{   /* a uniform random number generator */
    static unsigned s_Random = 1;
    s_Random = s_Random * 1103515245 + 12345;
    return (s_Random / 65536) % 32768;
}


extern int main(void)
{
#  define X_MAX_N_IO  (unsigned) 4
#  define X_MAX_READ  (size_t)   (3 * BUF_DEF_CHUNK_SIZE)
#  define X_TIMES     (unsigned) (s_Rand() % X_MAX_N_IO)
#  define X_BYTES     (size_t)   (s_Rand() % X_MAX_READ)

    BUF buf  = 0;
    BUF buf1 = 0;
    int/*bool*/ do_loop = 1 /* true */;

    /* a simple test */
    {{
        char charbuf[128];
        assert(BUF_PushBack(&buf, (const char*) "0", 1));
        assert(BUF_Write(&buf, (const char*) "1", 1));
        assert(BUF_Peek(buf, charbuf, sizeof(charbuf)) == 2);
        assert(BUF_PushBack(&buf, (const char*) "BB", 2));
        assert(BUF_PushBack(&buf, (const char*) "aa", 2));
        assert(BUF_Write(&buf, (const char*) "23", 3));
        assert(BUF_Read(buf, charbuf, sizeof(charbuf)) == 9);
        assert(strcmp(charbuf, (const char*) "aaBB0123") == 0);
        BUF_Destroy(buf);
        buf = 0;
    }}

    /* usage */
    fprintf(stderr, "Waiting for the data in STDIN...\n");

    /* read up to the very end of input stream */
    while ( do_loop ) {
        char charbuf[X_MAX_READ];
        unsigned i, n_times;

        /* read from the input stream, write to the NCBI IO-buf */
        n_times = X_TIMES;
        for (i = 0;  i < n_times;  i++) {
            size_t n_bytes = X_BYTES;
            if ( !n_bytes )
                continue;
            n_bytes = fread(charbuf, 1, n_bytes, stdin);
            fprintf(stderr, "STDIN     %5lu\n", (unsigned long) n_bytes);
            if ( !n_bytes ) {
                do_loop = 0 /* false */; /* end of the input stream */
                break;
            }
            assert(BUF_Write(&buf,  charbuf, n_bytes));
            assert(BUF_Write(&buf1, charbuf, n_bytes));
            fprintf(stderr, "BUF_Write %5lu\n", (unsigned long) n_bytes);
        }

        /* peek & read from the NCBI IO-buf, write to the output stream */
        n_times = X_TIMES;
        for (i = 0;  i < n_times  &&  BUF_Size(buf);  i++) {
            int/*bool*/ do_peek = (s_Rand() % 2 == 0);
            size_t n_peek = 0;
            size_t n_bytes = X_BYTES;
            if ( !n_bytes )
                continue;

            /* peek from the NCBI IO-buf */
            if ( do_peek ) {
                unsigned j, n_peek_times = s_Rand() % 3 + 1;
                for (j = 0;  j < n_peek_times;  j++) {
                    n_peek = BUF_Peek(buf, charbuf, n_bytes);
                    fprintf(stderr, "\tBUF_Peek %5lu\n",(unsigned long)n_peek);
                }
            }

            /* read (or just discard) the data */
            if (do_peek  &&  s_Rand() % 2 == 0)
                n_bytes = BUF_Read(buf, 0, n_bytes);
            else
                n_bytes = BUF_Read(buf, charbuf, n_bytes);

            fprintf(stderr, "\t\tBUF_Read %5lu\n", (unsigned long) n_bytes);
            assert(!do_peek  ||  n_bytes == n_peek);

            /* push back & re-read */
            if (s_Rand() % 3 == 0) {
                size_t n_pushback = s_Rand() % n_bytes;
                if (n_pushback == 0)
                    n_pushback = 1;
                assert(BUF_PushBack
                       (&buf, charbuf + n_bytes - n_pushback, n_pushback));
                assert(BUF_Read
                       (buf, charbuf + n_bytes - n_pushback, n_pushback));
            }

            /* write the read data to the output stream */
            assert(n_bytes == fwrite(charbuf, 1, n_bytes, stdout));
            fprintf(stderr, "\t\tSTDOUT   %5lu\n", (unsigned long) n_bytes);
        }
    }

    /* flush the IO-buf to the output stream */
    while ( BUF_Size(buf) ) {
        char charbuf[256];
        size_t n_bytes = BUF_Read(buf, charbuf, sizeof(charbuf));
        {{
            char   tmp[sizeof(charbuf)];
            size_t n_pushback = s_Rand() % 64;
            if (n_pushback > n_bytes)
                n_pushback = n_bytes;

            assert(BUF_PushBack
                   (&buf, charbuf + n_bytes - n_pushback, n_pushback));
            assert(BUF_Read
                   (buf, tmp, n_pushback) == n_pushback);
            memcpy(charbuf + n_bytes - n_pushback, tmp, n_pushback);
        }}
        fprintf(stderr, "\t\tBUF_Read/flush %5lu\n", (unsigned long) n_bytes);
        assert(n_bytes);
        assert(n_bytes == fwrite(charbuf, 1, n_bytes, stdout));
        fprintf(stderr, "\t\tSTDOUT  /flush %5lu\n", (unsigned long) n_bytes);
    }
    fflush(stdout);

    /* Test for "BUF_PeekAt()" */
    {{
        size_t buf1_size = BUF_Size(buf1);
        int n;
        assert(buf1_size > 0);

        for (n = 0;  n < 20;  n++) {
            size_t pos;
            size_t size;

            /* Erratically copy "buf1" to "buf" */
            for (pos = 0;  pos < buf1_size;  pos += size) {
                char temp_buf[BUF_DEF_CHUNK_SIZE * 2];
                size_t n_peeked;

                size = s_Rand() % (BUF_DEF_CHUNK_SIZE * 2);
                n_peeked = BUF_PeekAt(buf1, pos, temp_buf, size);
                if (pos + size <= buf1_size) {
                    assert(n_peeked == size);
                } else {
                    assert(n_peeked == buf1_size - pos);
                }
                assert(BUF_PeekAt(buf1, pos, temp_buf, size) == n_peeked);
                assert(BUF_Write(&buf, temp_buf, n_peeked));
            }

            /* Compare "buf" and "buf1"; empty up "buf" in process */
            assert(BUF_Size(buf1) == BUF_Size(buf));
            for (pos = 0;  pos < buf1_size;  pos += size) {
                char bb[1024];
                char b1[1024];
                assert(sizeof(bb) == sizeof(b1));

                size = BUF_Read(buf, bb, sizeof(bb));
                assert(BUF_PeekAt(buf1, pos, b1, size) == size);

                assert(size <= sizeof(b1));
                assert(memcmp(bb, b1, size) == 0);
            }

            assert(pos == buf1_size);
            assert(BUF_Size(buf1) == buf1_size);
            assert(BUF_Size(buf)  == 0);
        }
    }}


    /* cleanup & exit */
    BUF_Destroy(buf1);
    BUF_Destroy(buf);
    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * $Log: test_ncbi_buffer.c,v $
 * Revision 6.10  2003/02/25 16:26:58  lavr
 * Log moved to end; code re-indented
 *
 * Revision 6.9  2003/01/08 02:01:27  lavr
 * BUF_Destroy() made not returning a value
 *
 * Revision 6.8  2002/03/22 19:46:30  lavr
 * Test_assert.h made last among the include files
 *
 * Revision 6.7  2002/01/16 21:23:14  vakatov
 * Utilize header "test_assert.h" to switch on ASSERTs in the Release mode too
 *
 * Revision 6.6  2001/04/23 18:07:22  vakatov
 * + BUF_PeekAt()
 *
 * Revision 6.5  2000/03/24 23:12:12  vakatov
 * Starting the development quasi-branch to implement CONN API.
 * All development is performed in the NCBI C++ tree only, while
 * the NCBI C tree still contains "frozen" (see the last revision) code.
 *
 * Revision 6.4  2000/02/23 22:34:37  vakatov
 * Can work both "standalone" and as a part of NCBI C++ or C toolkits
 *
 * Revision 6.3  1999/11/22 16:12:54  vakatov
 * Allow to #include ncbi_buffer.h as a local
 *
 * Revision 6.2  1999/11/12 17:31:51  vakatov
 * Cosmetics -- get rid of an extra warning in Release build
 *
 * Revision 6.1  1999/10/12 16:33:50  vakatov
 * Moved all TEST suite code from "ncbi_buffer.c" to "test/test_ncbi_buffer.c"
 *
 * ===========================================================================
 */
