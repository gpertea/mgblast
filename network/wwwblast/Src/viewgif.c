/* $Id: viewgif.c,v 1.7 2005/12/28 20:56:14 merezhuk Exp $
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
* File Name:  $RCSfile: viewgif.c,v $
*
* Author:  Ilya Dondoshansky
*
* Initial Creation Date: 2002/12/02
*
* $Revision: 1.7 $
*
* File Description:
*        CGI program, part of standalone WWW Blast package, pipes graphic overview image back to user.
*
* $Log: viewgif.c,v $
* Revision 1.7  2005/12/28 20:56:14  merezhuk
* fix for broken gif;
*
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <fcntl.h>


static void SigAlrmHandler(int);
static void SigTermHandler(int);

static char FileName[128];

int main(void) 
{
    FILE *pp;
    char tmp_buff[1024];
    char *PidFile;
    char *server;
    int  bytes;
    char *ContentGif  =  "Content-type: image/gif\r\n\r\n";
    struct sigaction sa;
    sigset_t sigset;

    sigfillset(&sigset);
    sa.sa_mask = sigset;

    sa.sa_flags = SA_RESETHAND | SA_RESTART;
    sa.sa_handler = SigAlrmHandler;
    sigaction(SIGALRM, &sa, NULL);

    sa.sa_handler = SigTermHandler;
    sigaction(SIGTERM, &sa, NULL);
    sigaction(SIGPIPE, &sa, NULL);
    
    PidFile = (char *) getenv("QUERY_STRING"); 
    sprintf(FileName, "TmpGifs/%s", PidFile);
    
    server = (char *) getenv("SERVER_SOFTWARE");

    if((pp = fopen(FileName, "r")) == NULL) {
        /* Just do nothing */
        sprintf(tmp_buff, "HTTP/1.0 204 Not Modified\n");
        write(1, tmp_buff, strlen(tmp_buff));
        if (server)
           sprintf(tmp_buff, "Server: %s\n", server);
        write(1, tmp_buff, strlen(tmp_buff));
        sprintf(tmp_buff, "MIME-Version: 1.0\n");
        write(1, tmp_buff, strlen(tmp_buff));
        write(1, ContentGif, strlen(ContentGif)); 
    } else {
        sprintf(tmp_buff, "HTTP/1.0 200 OK\r\n");
        write(1, tmp_buff, strlen(tmp_buff));
        if (server)
           sprintf(tmp_buff, "Server: %s\n", server);
        write(1, tmp_buff, strlen(tmp_buff));
        sprintf(tmp_buff, "MIME-Version: 1.0\r\n");
        write(1, tmp_buff, strlen(tmp_buff));
        write(1, ContentGif, strlen(ContentGif)); 
        
        while ((bytes =fread(tmp_buff, 1, 256, pp)) >0)
            write(1, tmp_buff, bytes); 
    }
    remove(FileName); 
    return 0;
}
static void SigAlrmHandler(int id)
{
    
    char tmp_buff[1024];
    char *ContentGif  =  "Content-type: image/gif\r\n\r\n";
    
    sprintf(tmp_buff, "HTTP/1.0 204 Not Modified\n");
    write(1, tmp_buff, strlen(tmp_buff));
    sprintf(tmp_buff, "Server: %s\n", (char *) getenv("SERVER_SOFTWARE"));
    write(1, tmp_buff, strlen(tmp_buff));
    sprintf(tmp_buff, "MIME-Version: 1.0\n");
    write(1, tmp_buff, strlen(tmp_buff));
    write(1, ContentGif, strlen(ContentGif)); 
    remove(FileName); 
    exit(1);   
}

static void SigTermHandler(int id)
{
    remove(FileName); 
    exit(1);
}















